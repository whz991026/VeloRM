##' Estimate RNA velocity using gene-relative slopes
##'
##' @param emat - spliced (exonic) count matrix
##' @param nmat - unspliced (nascent) count matrix
##' @param deltaT - amount of time to project the cell forward
##' @param smat - optional spanning read matrix (used in offset calculations)
##' @param steady.state.cells - optional set of steady-state cells on which the gamma should be estimated (defaults to all cells)
##' @param kCells - number of k nearest neighbors (NN) to use in slope calculation smoothing
##' @param cellKNN - optional pre-calculated cell KNN matrix
##' @param kGenes - number of genes (k) to use in gene kNN pooling
##' @param geneKNN - optional pre-calculated gene KNN matrix
##' @param old.fit - optional old result (in this case the slopes and offsets won't be recalculated, and the same kNN graphs will be used)
##' @param mult - library scaling factor (1e6 in case of FPM)
##' @param min.nmat.smat.correlation - minimum required Spearman rank correlation between n and s counts of a gene
##' @param min.nmat.emat.correlation - minimum required Spearman rank correlation between n and e counts of a gene
##' @param min.nmat.emat.slope - minimum sloope of n~e regression
##' @param zero.offset - should offset be set to zero, or determined (through smat regression or using near-0 e cases)
##' @param deltaT2 - scaling of the projected difference vector (normally should be set to 1)
##' @param delta_model estimate the delat use model 1 or model 2
##' @param fit.quantile perform gamma fit on a top/bottom quantiles of expression magnitudes
##' @param diagonal.quantiles whether extreme quantiles should be computed diagonally
##' @param show.gene an optional name of a gene for which the velocity estimation details should be shown (instead of estimating all velocities)
##' @param do.par whether the graphical device parameters should be reset as part of show.gene (default=TRUE)
##' @param cell.dist - cell distance to use in cell kNN pooling calculations
##' @param emat.size - pre-calculated cell sizes for the emat (spliced) matrix
##' @param nmat.size - pre-calculated cell sizes for the nmat (unspliced) matrix
##' @param cell.emb - cell embedding to be used in show.gene function
##' @param cell.colors - cell colors to be used in show.gene function
##' @param expression.gradient - color palette used to show the expression magnitudes in show.gene function
##' @param residual.gradient - color palette used to show the u residuals in show.gene function
##' @param n.cores - number of cores to use
##' @param verbose - output messages about progress
##'
##' @import MASS
##' @import Matrix
##' @import rjags
##' @importFrom methods as
##' @importFrom stats as.dist coef dist lm median na.omit predict quantile resid sd
##' @importFrom graphics abline box lines par
##'
##' @return a list with velocity results, including the current normalized expression state ($current), projected ($projected) over a certain time ($deltaT), unscaled transcriptional change ($deltaE), fit results ($gamma, $ko, $sfit if spanning reads were used), optional cell pooling parameters ($cellKNN, $kCells), kNN-convolved normalized matrices (conv.nmat.norm and conv.emat.norm), library scale ($mult)
##' @examples
##' \dontrun{
##'  # use min/max quantile gamma fit (recommended option when one can afford to do cell kNN smoothing)
##'  # The example below uses k=5 cell kNN pooling, and top/bottom 2% exprssion quantiles
##'  # emat and nmat are spliced (exonic) and unspliced (intronic) molecule/read count matirces
##' (preferably filtered for informative genes)
##'  rvel <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells = 5,fit.quantile = 0.02)
##'
##'  # alternativly, the function can be used to visualize gamma fit and regression for a
##' particular gene. here we pass embedding (a matrix/data frame with rows named with cell names,
##' and columns corresponding to the x/y coordinates)
##'
##'  # and cell colors. old.fit is used to save calculation time.
##'  gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells = 5,fit.quantile = 0.02,
##'     old.fit=rvel,show.gene='Chga',cell.emb=emb,cell.colors=cell.colors)
##' }
##' @export
gene.relative.velocity.estimates <- function (
        emat, nmat, deltaT = 1, smat = NULL, steady.state.cells = colnames(emat),
        kCells = 10, cellKNN = NULL, kGenes = 1, geneKNN = NULL, old.fit = NULL,
        mult = 1000, min.nmat.smat.correlation = 0.05, min.nmat.emat.correlation = 0.05,
        min.nmat.emat.slope = 0.05, zero.offset = FALSE, deltaT2 = 1,delta_model="model 2",
        fit.quantile = NULL, diagonal.quantiles = FALSE, show.gene = NULL,
        do.par = TRUE, cell.dist = NULL, emat.size = NULL, nmat.size = NULL,
        cell.emb = NULL, cell.colors = NULL, expression.gradient = NULL,
        residual.gradient = NULL, n.cores = 1, verbose = TRUE)
{

  if (!all(colnames(emat) == colnames(nmat)))
    stop("emat and nmat must have the same columns (cells)")
  if (!is.null(smat)) {
    if (!all(colnames(emat) == colnames(smat)))
      stop("smat must have the same columns (cells) as emat")
  }
  resl <- list()
  vg <- intersect(rownames(emat), rownames(nmat))
  if (is.null(smat)) {
    emat <- emat[vg, ]
    nmat <- nmat[vg, ]
  }
  else {
    vg <- intersect(vg, rownames(smat))
    emat <- emat[vg, ]
    nmat <- nmat[vg, ]
    smat <- smat[vg, ]
  }
  if (!is.null(show.gene)) {
    if (!show.gene %in% rownames(emat)) {
      stop(paste("gene", show.gene, "is not present in the filtered expression matrices"))
    }
  }
  pcount <- 1
  if (!is.null(cell.dist)) {
    if (!inherits(cell.dist, "dist")) {
      stop("cell.dist must be of a class dist")
    }
    if (!all(labels(cell.dist) == colnames(emat))) {
      cat("matching cells between cell.dist and emat/nmat ... ")
      cell.dist <- as.matrix(cell.dist)
      cn <- intersect(colnames(emat), colnames(cell.dist))
      cell.dist <- as.dist(cell.dist[cn, cn])
      emat <- emat[, cn]
      nmat <- nmat[, cn]
      if (!is.null(smat)) {
        smat <- smat[, cn]
      }
      cat("done\n")
    }
  }
  if (is.null(emat.size)) {
    emat.size <- Matrix::colSums(emat)
  }
  if (is.null(nmat.size)) {
    nmat.size <- Matrix::colSums(nmat)
  }
  emat.cs <- emat.size[colnames(emat)]/mult
  nmat.cs <- nmat.size[colnames(nmat)]/mult
  emat.log.norm <- log(as.matrix(t(t(emat)/emat.cs)) + pcount)
  if (!is.null(old.fit)) {
    cellKNN <- old.fit[["cellKNN"]]
  }
  knn.maxl <- 100
  if (kCells > 1) {
    if (is.null(cellKNN)) {
      cat("calculating cell knn ... ")
      if (is.null(cell.dist)) {
        cellKNN <- balanced_knn(as.matrix(cor(emat.log.norm)),
                                kCells,kCells * knn.maxl,method="cor",n.cores= n.cores)
      }
      else {
        cellKNN <- balanced_knn(as.matrix(cell.dist), kCells,
                               kCells * knn.maxl,method=attributes(cell.dist)$method,
                               n.cores= n.cores)
      }
      diag(cellKNN) <- 1
      colnames(cellKNN) <- rownames(cellKNN) <-colnames(emat)
      resl$cellKNN <- cellKNN
      cat("done\n")
    }
    rm(emat.log.norm)
    cat("calculating convolved matrices ... ")
    conv.emat <- emat %*% cellKNN[colnames(emat), colnames(emat)]
    conv.nmat <- nmat %*% cellKNN[colnames(nmat), colnames(nmat)]
    conv.emat.cs <- (emat.cs %*% cellKNN[colnames(emat),
                                         colnames(emat)])[1, ]
    conv.nmat.cs <- (nmat.cs %*% cellKNN[colnames(nmat),
                                         colnames(nmat)])[1, ]
    cat("done\n")
  }
  else {
    conv.emat <- emat
    conv.nmat <- nmat
    cellKNN <- NULL
    conv.emat.cs <- emat.cs
    conv.nmat.cs <- nmat.cs
  }
  conv.emat.norm <- t(t(conv.emat)/conv.emat.cs)
  conv.nmat.norm <- t(t(conv.nmat)/conv.nmat.cs)
  emat.norm <- t(t(emat)/emat.cs)
  nmat.norm <- t(t(nmat)/nmat.cs)
  if (kGenes > 1) {
    if (is.null(geneKNN)){
      if (!is.null(old.fit) && !is.null(old.fit$geneKNN)) {
        geneKNN <- old.fit$geneKNN
      }
      else {
        cat("gene kNN ... ")
        geneKNN <- balanced_knn(as.matrix(dist((log(as.matrix(conv.emat.norm) +
                                    pcount)))), kGenes, kGenes * 1200, n.cores =n.cores)
        diag(geneKNN) <- 1
      }
    }
    
    resl$geneKNN <- geneKNN
    cat("Bayesian linear regression for each gene ... ")
    
    scaledGeneKNN <- geneKNN/rowSums(geneKNN)
    cat("calculate the prior information of the bayesian regression ... ")
    conv.emat.norm_mean <- scaledGeneKNN %*% conv.emat.norm
    conv.nmat.norm_mean <- scaledGeneKNN %*% conv.nmat.norm
    gama_prior_mean <- apply(conv.nmat.norm_mean/conv.emat.norm_mean, 1, function(x) mean(x,na.rm=TRUE))
    gama_prior_sd <- apply(conv.nmat.norm_mean/conv.emat.norm_mean, 1, function(x) sd(x,na.rm=TRUE))
    names(gama_prior_mean) <- names(gama_prior_sd) <- rownames(conv.emat.norm)
    cat("done\n")
  }
  if (!is.null(smat)) {
    if (kCells > 1) {
      conv.smat <- smat %*% cellKNN[colnames(smat), colnames(smat)]
    }
    else {
      conv.smat <- smat
    }
    conv.smat.cs <- Matrix::colSums(conv.smat)/mult
    conv.smat.norm <- t(t(conv.smat)/conv.smat.cs)
    colnames(conv.smat.norm) <- colnames(conv.emat.norm)
    rownames(conv.smat.norm) <- rownames(conv.emat.norm)
    if (kGenes > 1) {
      conv.smat.norm <- scaledGeneKNN %*% conv.smat.norm
      colnames(conv.smat.norm) <- colnames(conv.emat.norm)
      rownames(conv.smat.norm) <- rownames(conv.emat.norm)
    }
    if (is.null(old.fit)) {
      cat("fitting smat-based offsets ... ")
      sfit <- data.frame(do.call(rbind, parallel::mclapply(sn(rownames(conv.emat.norm)),
               function(gn) {
                 df <- data.frame(n = (conv.nmat.norm[gn, steady.state.cells]),
                                  e = (conv.emat.norm[gn, steady.state.cells]),
                                  s = conv.smat.norm[gn, steady.state.cells])
                 sd <- lm(n ~ s, data = df)
                 r <- with(df[df$s > 0, ], cor(n, s, method = "spearman"),
                           3)
                 return(c(o = pmax(0, as.numeric(sd$coef[1])),
                          s = as.numeric(sd$coef[2]), r = r))
               }, mc.cores = n.cores, mc.preschedule = T)))
      cat("done\n")
    }
    else {
      sfit <- old.fit$sfit
    }
  }
  resl$conv.nmat.norm <- conv.nmat.norm
  resl$conv.emat.norm <- conv.emat.norm
  if (!is.null(show.gene)) {
    gn <- show.gene
    if (!is.null(cell.emb)) {
      cc <- intersect(rownames(cell.emb), colnames(conv.emat.norm))
      if (do.par) {
        par(mfrow = c(1, 4), mar = c(2.5, 2.5, 2.5, 0.5),
            mgp = c(1.5, 0.65, 0), cex = 0.85)
      }
      plot(cell.emb[cc, ], pch = 21, col = ac(1, alpha = 0.2),
           bg = val2col(conv.emat.norm[gn, cc], gradientPalette = expression.gradient),
           cex = 0.8, xlab = "", ylab = "", main = paste(gn,
                                                         "s"), axes = F)
      box()
      plot(cell.emb[cc, ], pch = 21, col = ac(1, alpha = 0.2),
           bg = val2col(conv.nmat.norm[gn, cc], gradientPalette = expression.gradient),
           cex = 0.8, xlab = "", ylab = "", main = paste(gn,
                                                         "u"), axes = F)
      box()
    }
    do <- NULL
    if (!is.null(smat)) {
      df <- data.frame(n = (conv.nmat.norm[gn, steady.state.cells]),
                       e = (conv.emat.norm[gn, steady.state.cells]),
                       o = sfit[gn, "o"])
      if (zero.offset)
        df$o <- 0
    }
    else {
      df <- data.frame(n = (conv.nmat.norm[gn, steady.state.cells]),
                       e = (conv.emat.norm[gn, steady.state.cells]))
      o <- 0
      df$o <- o
      if (!zero.offset) {
        zi <- df$e < 1/conv.emat.cs[steady.state.cells]
        if (any(zi)) {
          o <- sum(df$n[zi])/(sum(zi) + 1)
        }
      }
      df$o <- o
    }
    d <- lm(n ~ e + offset(o) + 0, data = df, weights = df$e^4 +
              df$n^4)
    cell.col <- ac(rep(1, nrow(df)), alpha = 0.1)
    names(cell.col) <- rownames(df)
    if (!is.null(cell.colors)) {
      cc <- intersect(names(cell.colors), rownames(df))
      cell.col[cc] <- cell.colors[cc]
    }
    plot(df$e, df$n, pch = 21, bg = ac(cell.col, alpha = 0.3),
         col = ac(1, alpha = 0.1), cex = 0.8, xlab = "s",
         ylab = "u", main = paste(gn, "fit"))
    if (!is.null(do)) {
      abline(do, lty = 2, col = 8)
    }
    if (!is.null(fit.quantile)) {
      if (diagonal.quantiles) {
        emax <- quantile(df$e, p = 0.99)
        nmax <- quantile(df$n, p = 0.99)
        if (emax == 0)
          emax <- max(max(df$e), 0.001)
        if (nmax == 0)
          nmax <- max(max(df$n), 0.001)
        x <- df$e/emax + df$n/nmax
        eq <- quantile(x, p = c(fit.quantile, 1 - fit.quantile))
        if (!is.null(smat)) {
          pw <- as.numeric(x >= eq[2])
        }
        else {
          pw <- as.numeric(x >= eq[2] | x <= eq[1])
        }
      }
      else {
        eq <- quantile(df$e, p = c(fit.quantile, 1 -
                                     fit.quantile))
        if (!is.null(smat) || zero.offset) {
          pw <- as.numeric(df$e >= eq[2])
        }
        else {
          pw <- as.numeric(df$e >= eq[2] | df$e <= eq[1])
        }
      }
      if (!is.null(smat) || zero.offset) {
        d <- lm(n ~ e + offset(o) + 0, data = df, weights = pw)
      }
      else {
        d <- lm(n ~ e, data = df, weights = pw)
      }
    }
    df <- df[order(df$e, decreasing = T), ]
    lines(df$e, predict(d, newdata = df), lty = 2, col = 2)
    if (!is.null(cell.emb)) {
      plot(cell.emb[cc, ], pch = 21, col = ac(1, alpha = 0.2),
           bg = val2col(resid(d)[cc], gradientPalette = residual.gradient),
           cex = 0.8, xlab = "", ylab = "", main = paste(gn,
                                                         "resid"), axes = F)
      box()
    }
    if (kGenes > 1) {
      return(invisible(geneKNN))
    }
    else {
      return(1)
    }
  }
  cat("fitting gamma coefficients ... ")
  if (is.null(old.fit)) {
    ko <- data.frame(do.call(rbind, parallel::mclapply(sn(rownames(conv.emat.norm)),
           function(gn) {
             if (!is.null(smat)) {
               df <- data.frame(n = (conv.nmat.norm[gn, steady.state.cells]),
                                e = (conv.emat.norm[gn, steady.state.cells]),
                                o = sfit[gn, "o"])
               if (zero.offset)
                 df$o <- 0
             }
             else {
               df <- data.frame(n = (conv.nmat.norm[gn, steady.state.cells]),
                                e = (conv.emat.norm[gn, steady.state.cells]))
               o <- 0
               if (!zero.offset) {
                 zi <- df$e < 1/conv.emat.cs[steady.state.cells]
                 if (any(zi)) {
                   o <- sum(df$n[zi])/(sum(zi) + 1)
                 }
               }
               df$o <- o
             }
             if (is.null(fit.quantile)) {
                 d <- lm(n ~ e + offset(o) + 0, data = df, weights = df$e^4 +
                           df$n^4)
                 return(c(o = df$o[1], g = as.numeric(coef(d)[1]),
                          r = cor(df$e, df$n, method = "spearman")))
             }
             else {
               if (diagonal.quantiles) {
                 emax <- quantile(df$e, p = 0.99)
                 nmax <- quantile(df$n, p = 0.99)
                 if (emax == 0)
                   emax <- max(max(df$e), 0.001)
                 if (nmax == 0)
                   nmax <- max(max(df$n), 0.001)
                 x <- df$e/emax + df$n/nmax
                 eq <- quantile(x, p = c(fit.quantile, 1 -
                                           fit.quantile))
                 if (!is.null(smat)) {
                   pw <- as.numeric(x >= eq[2])
                 }
                 else {
                   pw <- as.numeric(x >= eq[2] | x <= eq[1])
                 }
               }
               else {
                 eq <- quantile(df$e, p = c(fit.quantile,
                                            1 - fit.quantile))
                 if (!is.null(smat) || zero.offset) {
                   pw <- as.numeric(df$e >= eq[2])
                 }
                 else {
                   pw <- as.numeric(df$e >= eq[2] | df$e <=
                                      eq[1])
                 }
               }
               if (!is.null(smat) || zero.offset) {
                 d <- lm(n ~ e + offset(o) + 0, data = df,
                         weights = pw)
                 return(c(o = df$o[1], g = as.numeric(coef(d)[1]),
                          r = cor(df$e, df$n, method = "spearman")))
               }
               else {
                 d <- lm(n ~ e, data = df, weights = pw)
                 return(c(o = as.numeric(coef(d)[1]), g = as.numeric(coef(d)[2]),
                          r = cor(df$e, df$n, method = "spearman")))
               }
             }
             
             if (kGenes>1){
               mean_prior <- gama_prior_mean[gn]
               sd_prior <- gama_prior_sd[gn]
               data_list <- list(y = df$n-df$o, x = df$e)
               model_string <- "model {
                   for (i in 1:N) {
                     y[i] ~ dnorm(mu[i], tau)
                     mu[i] <- beta0 + beta1 * x[i]
                   }
                   beta0 <-0
                   beta1 ~ dnorm(mean_prior, sd_prior)
                   tau ~ dgamma(0.001, 0.001)
                 }"
               data_list <- list(y = df$n-df$o, x = df$e, N = nrow(df),
                                 mean_prior=mean_prior,sd_prior=sd_prior)
               model <- rjags::jags.model(textConnection(model_string), data = data_list, n.chains = 3)
               update(model, 1000)
               samples <- rjags::coda.samples(model, variable.names =  "beta1", n.iter = 5000)
               samples_s <- summary(samples)
               return(c(o = df$o[1], g = as.numeric(samples_s$statistics[1]),
                        r = cor(df$e, df$n, method = "spearman")))
             }
             
           }, mc.cores = n.cores, mc.preschedule = T)))
    ko <- na.omit(ko)
    cat("done. succesfful fit for", nrow(ko), "genes\n")
  }
  else {
    full.ko <- ko <- na.omit(old.fit$ko)
  }
  if (!is.null(smat)) {
    sfit <- na.omit(sfit)
    ko <- ko[rownames(ko) %in% rownames(sfit), ]
    vi <- sfit$r > min.nmat.smat.correlation
    ko <- ko[vi, ]
    if (!all(vi))
      cat("filtered out", sum(!vi), "out of", length(vi),
          "genes due to low nmat-smat correlation\n")
  }
  full.ko <- ko
  vi <- ko$r > min.nmat.emat.correlation
  if (!all(vi))
    cat("filtered out", sum(!vi), "out of", length(vi), "genes due to low nmat-emat correlation\n")
  ko <- ko[vi, ]
  vi <- ko$g > min.nmat.emat.slope
  if (!all(vi))
    cat("filtered out", sum(!vi), "out of", length(vi), "genes due to low nmat-emat slope\n")
  ko <- ko[vi, ]
  gamma <- ko$g
  offset <- ko$o
  names(gamma) <- names(offset) <- rownames(ko)
  cat("calculating RNA velocity shift ... ")
  
  if (delta_model=="model 1"){
    cat("calculating delta use model 1 ... ")
    deltaE <- t.get.projected.delta.model1(conv.emat.norm, conv.nmat.norm,
                                    gamma, offset = offset, delta = deltaT)
  }
  if (delta_model=="model 2"){
    cat("calculating delta use model 2 ... ")
    deltaE <- t.get.projected.delta.model2(conv.emat.norm, conv.nmat.norm,
                                    gamma, offset = offset, delta = deltaT)
  }

  resl$gamma <- gamma
  cat("done\n")
  cat("calculating extrapolated cell state ... ")
  emat.norm <- emat[rownames(emat) %in% rownames(deltaE), ]
  emat.sz <- emat.cs
  emat.norm <- t(t(emat.norm)/(emat.sz))
  emn <- t.get.projected.cell2(emat.norm, emat.sz, as.matrix(deltaE),
                               mult = mult, delta = deltaT2)
  cat("done\n")
  full.ko$valid <- rownames(full.ko) %in% rownames(ko)
  resl <- c(resl, list(projected = emn, current = emat.norm,
                       deltaE = deltaE, deltaT = deltaT, ko = full.ko, mult = mult,
                       kCells = kCells,kGenes = kGenes))
  if (!is.null(smat)) {
    resl$sfit <- sfit
  }
  return(resl)
}
