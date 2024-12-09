##' Estimate RNA velocity using site-relative slopes
##'
##' @param emat - spliced (exonic) count matrix
##' @param nmat - unspliced (nascent) count matrix
##' @param deltaT - amount of time to project the cell forward
##' @param steady.state.cells - optional set of steady-state cells on which the gamma should be estimated (defaults to all cells)
##' @param kCells - number of k nearest neighbors (NN) to use in slope calculation smoothing
##' @param cellKNN - optional pre-calculated cell KNN matrix
##' @param kSites - number of sites (k) to use in site kNN pooling
##' @param siteKNN - optional pre-calculated site KNN matrix
##' @param old.fit - optional old result (in this case the slopes and offsets won't be recalculated, and the same kNN graphs will be used)
##' @param mult - library scaling factor (1e6 in case of FPM)
##' @param min.nmat.emat.correlation - minimum required Spearman rank correlation between n and e counts of a site
##' @param min.nmat.emat.slope - minimum sloope of n~e regression
##' @param zero.offset - should offset be set to zero
##' @param deltaT2 - scaling of the projected difference vector (normally should be set to 1)
##' @param delta_model estimate the delat use model 1 or model 2
##' @param fit.quantile perform gamma fit on a top/bottom quantiles of expression magnitudes
##' @param diagonal.quantiles whether extreme quantiles should be computed diagonally
##' @param show.site an optional name of a site for which the velocity estimation details should be shown (instead of estimating all velocities)
##' @param cell.dist - cell distance to use in cell kNN pooling calculations
##' @param emat.size - pre-calculated cell sizes for the emat (spliced) matrix
##' @param nmat.size - pre-calculated cell sizes for the nmat (unspliced) matrix
##' @param cell.emb - cell embedding to be used in show.site function
##' @param cell.colors - cell colors to be used in show.site function
##' @param expression.gradient - color palette used to show the expression magnitudes in show.site function
##' @param residual.gradient - color palette used to show the u residuals in show.site function
##' @param n.cores - number of cores to use
##' @param verbose - output messages about progress
##' @param nrows the hyper-parameter for show site plot the number of rows
##' @param low_color the hyper-parameter for show site low color
##' @param high_color the hyper-parameter for show site high color
##' @param point.size size of the point
##'
##' @import MASS
##' @import Matrix
##' @import rjags
##' @import cowplot
##' @importFrom methods as
##' @importFrom stats as.dist coef dist lm median na.omit predict quantile resid sd
##' @importFrom graphics abline box lines par
##'
##' @return a list with velocity results, including the current normalized expression state ($current), projected ($projected) over a certain time ($deltaT), unscaled transcriptional change ($deltaE), fit results ($gamma, $ko, $sfit if spanning reads were used), optional cell pooling parameters ($cellKNN, $kCells), kNN-convolved normalized matrices (conv.nmat.norm and conv.emat.norm), library scale ($mult)
##' 
##' @export
site.relative.velocity.estimates <- function (
        emat, nmat, deltaT = 1, steady.state.cells = colnames(emat),
        kCells = 10, cellKNN = NULL, kSites = 1, siteKNN = NULL, old.fit = NULL,
        mult = 1000, min.nmat.emat.correlation = 0.05,
        min.nmat.emat.slope = 0.05, zero.offset = FALSE, deltaT2 = 1,delta_model="model 2",
        fit.quantile = NULL, diagonal.quantiles = FALSE, show.site = NULL,
        cell.dist = NULL, emat.size = NULL, nmat.size = NULL,
        cell.emb = NULL, cell.colors = NULL, expression.gradient = NULL,
        residual.gradient = NULL, n.cores = 1, verbose = TRUE,nrows=2,
        low_color="#FFCC15",high_color="#9933FF",point.size=3)
{
  point_shape=21
  if (!all(colnames(emat) == colnames(nmat)))
    stop("emat and nmat must have the same columns (cells)")
 
  resl <- list()
  vg <- intersect(rownames(emat), rownames(nmat))
  
    emat <- emat[vg, ]
    nmat <- nmat[vg, ]
  
  if (!is.null(show.site)) {
    if (!show.site %in% rownames(emat)) {
      stop(paste("site", show.site, "is not present in the filtered expression matrices"))
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
  if (kSites > 1) {
    if (is.null(siteKNN)){
      if (!is.null(old.fit) && !is.null(old.fit$siteKNN)) {
        siteKNN <- old.fit$siteKNN
      }
      else {
        cat("site kNN ... ")
        siteKNN <- balanced_knn(as.matrix(dist((log(as.matrix(conv.emat.norm) +
                                    pcount)))), kSites, kSites * 1200, n.cores =n.cores)
        diag(siteKNN) <- 1
      }
    }
    
    resl$siteKNN <- siteKNN
    cat("Bayesian linear regression for each site ... ")
    
    scaledSiteKNN <- siteKNN/rowSums(siteKNN)
    cat("calculate the prior information of the bayesian regression ... ")
    conv.emat.norm_mean <- scaledSiteKNN %*% conv.emat.norm
    conv.nmat.norm_mean <- scaledSiteKNN %*% conv.nmat.norm
    gama_prior_mean <- apply(conv.nmat.norm_mean/conv.emat.norm_mean, 1, function(x) mean(x,na.rm=TRUE))
    gama_prior_sd <- apply(conv.nmat.norm_mean/conv.emat.norm_mean, 1, function(x) sd(x,na.rm=TRUE))
    names(gama_prior_mean) <- names(gama_prior_sd) <- rownames(conv.emat.norm)
    cat("done\n")
  }

  resl$conv.nmat.norm <- conv.nmat.norm
  resl$conv.emat.norm <- conv.emat.norm
  if (!is.null(show.site)) {
    gn <- show.site
    if (!is.null(cell.emb)) {
      cc <- intersect(rownames(cell.emb), colnames(conv.emat.norm))
      
      data.frame.s.ggplot2 <- as.data.frame(cell.emb[cc, 1])
      data.frame.s.ggplot2$PC1 <- cell.emb[cc, 1]
      data.frame.s.ggplot2$PC2 <- cell.emb[cc, 2]
      data.frame.s.ggplot2$color <- conv.emat.norm[gn, cc]
      data.frame.s.ggplot2 <- data.frame.s.ggplot2[,c("PC1","PC2","color")]
      
      colnames(data.frame.s.ggplot2) <- c("PC1","PC2","color")
      spliced.plot <- ggplot(data.frame.s.ggplot2)+
        aes(x = .data$PC1, y = .data$PC2, fill = .data$color) + 
        geom_point(alpha=1,shape=point_shape,size=point.size) + labs(x = "PC1",y = "PC2",title = paste0(show.site," s"),colour="s")+ 
        theme(plot.title = element_text(size=12,hjust=0.5))+
        scale_fill_gradient(low = low_color, high = high_color)
      
      data.frame.u.ggplot2 <- as.data.frame(cell.emb[cc, 1])
      data.frame.u.ggplot2$PC1 <- cell.emb[cc, 1]
      data.frame.u.ggplot2$PC2 <- cell.emb[cc, 2]
      data.frame.u.ggplot2$color <- conv.nmat.norm[gn, cc]
      data.frame.u.ggplot2 <- data.frame.u.ggplot2[,c("PC1","PC2","color")]
      colnames(data.frame.u.ggplot2) <- c("PC1","PC2","color")
      unspliced.plot <- ggplot(data.frame.s.ggplot2)+
        aes(x = .data$PC1, y = .data$PC2, fill = .data$color) + 
        geom_point(alpha=1,shape=point_shape,size=point.size) + labs(x = "PC1",y = "PC2",title = paste0(show.site," u"),colour="u")+ 
        theme(plot.title = element_text(size=12,hjust=0.5))+
        scale_fill_gradient(low = low_color, high = high_color)
      
    }

    
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
    
    
    d <- lm(n ~ e + offset(o) + 0, data = df, weights = df$e^4 +
              df$n^4)
    cell.col <- ac(rep(1, nrow(df)), alpha = 0.1)
    names(cell.col) <- rownames(df)
    if (!is.null(cell.colors)) {
      cc <- intersect(names(cell.colors), rownames(df))
      cell.col[cc] <- cell.colors[cc]
    }
    
    data.frame.fit.ggplot2 <- df
    data.frame.fit.ggplot2$PC1 <-as.numeric(df$e)
    data.frame.fit.ggplot2$PC2 <-as.numeric(df$n)
    data.frame.fit.ggplot2$color <- ac(cell.col, alpha = 0.3)
    data.frame.fit.ggplot2 <- data.frame.fit.ggplot2[,c("PC1","PC2","color")]
    
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

        pw <- as.numeric(x >= eq[2] | x <= eq[1])
        
      }
      else {
        eq <- quantile(df$e, p = c(fit.quantile, 1 -
                                     fit.quantile))
        
          pw <- as.numeric(df$e >= eq[2] | df$e <= eq[1])
        
      }
      
        d <- lm(n ~ e, data = df, weights = pw)
    }
    df <- df[order(df$e, decreasing = T), ]
    data.frame.fit.line.ggplot2 <- df
    data.frame.fit.line.ggplot2$pred <-  predict(d, newdata = df)
    data.frame.fit.line.ggplot2 <- data.frame.fit.line.ggplot2[,c("e","pred")]
    
    colnames(data.frame.fit.ggplot2) <- c("PC1","PC2","color")
    colnames(data.frame.fit.line.ggplot2) <- c("e","pred")
    fit.plot <- ggplot()+
        geom_point(data=data.frame.fit.ggplot2,aes(x = .data$PC1, y = .data$PC2,
                  fill = .data$color), alpha = 0.5,shape=point_shape,size=point.size) +
      geom_line(data=data.frame.fit.line.ggplot2,aes(x = .data$e, y = .data$pred),
                color = "red",linetype= 2)+
      labs(x = "s",y = "u",title = paste0(show.site," fit"))+ 
      theme( plot.title = element_text(size=12,hjust=0.5))#legend.position =" none ",
    
    if (!is.null(cell.emb)) {
      data.frame.resid.ggplot2 <- as.data.frame(cell.emb[cc, 1])
      data.frame.resid.ggplot2$PC1 <- cell.emb[cc, 1]
      data.frame.resid.ggplot2$PC2 <- cell.emb[cc, 2]
      data.frame.resid.ggplot2$color <- resid(d)[cc]
      data.frame.resid.ggplot2 <- data.frame.resid.ggplot2[, c("PC1","PC2","color")]
      colnames(data.frame.resid.ggplot2) <- c("PC1","PC2","color")
      resid.plot <- ggplot(data.frame.resid.ggplot2)+
        aes(x = .data$PC1, y = .data$PC2, fill = .data$color) + 
        geom_point(alpha=1,shape=point_shape,size=point.size) + labs(x = "PC1",y = "PC2",title = paste0(show.site," resid"),colour="resid")+ 
        theme(plot.title = element_text(size=12,hjust=0.5))+ 
        scale_fill_gradient(low = low_color, high = high_color)
      
      combine_plot <- cowplot::plot_grid(spliced.plot, unspliced.plot, fit.plot, 
                                         resid.plot, nrow=nrows,labels = LETTERS[1:4])#spliced.plot, unspliced.plot,
      return(combine_plot)
      } else{return(fit.plot)}
  }
  cat("fitting gamma coefficients ... ")
  if (is.null(old.fit)) {
    ko <- data.frame(do.call(rbind, parallel::mclapply(sn(rownames(conv.emat.norm)),
           function(gn) {
             
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
             
             if (is.null(fit.quantile)) {
                 d <- lm(n ~ e + offset(o) + 0, data = df, weights = df$e^4 +
                           df$n^4)
                 pw <- rep(df$n,length(df$n))
                 return(c(o = df$o[1], g = as.numeric(coef(d)[1]),
                          r = cor(df$e, df$n, method = "spearman"),pw))
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
                 
                   pw <- as.numeric(x >= eq[2] | x <= eq[1])
                 
               }
               else {
                 eq <- quantile(df$e, p = c(fit.quantile,
                                            1 - fit.quantile))
                 
                   pw <- as.numeric(df$e >= eq[2] | df$e <=
                                      eq[1])
                 
               }
               
                 d <- lm(n ~ e, data = df, weights = pw)
                 return(c(o = as.numeric(coef(d)[1]), g = as.numeric(coef(d)[2]),
                          r = cor(df$e, df$n, method = "spearman"),pw))
               
             }
             
             if (kSites>1){
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
               pw <- rep(df$n,length(df$n))
               return(c(o = df$o[1], g = as.numeric(samples_s$statistics[1]),
                        r = cor(df$e, df$n, method = "spearman"),pw))
             }
             
           }, mc.cores = n.cores, mc.preschedule = T)))
    
    ko <- na.omit(ko)
    cat("done. succesfful fit for", nrow(ko), "sites\n")
  }
  else {
    full.ko <- ko <- na.omit(old.fit$ko)
  }
  
  full.ko <- ko
  vi <- ko$r > min.nmat.emat.correlation
  if (!all(vi))
    cat("filtered out", sum(!vi), "out of", length(vi), "sites due to low nmat-emat correlation\n")
  ko <- ko[vi, ]
  vi <- ko$g > min.nmat.emat.slope
  if (!all(vi))
    cat("filtered out", sum(!vi), "out of", length(vi), "sites due to low nmat-emat slope\n")
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
  colnames(full.ko)[-c(1,2,3,length(colnames(full.ko)))] <- paste0("cell_",1:(length(colnames(full.ko))-4))
  pw.data.frame <- full.ko[,paste0("cell_",1:(length(colnames(full.ko))-4))]
  cell.use.list <- apply(pw.data.frame,1,function(x) which(x==1))
  full.ko <-full.ko[,c("o","g","r","valid")]
  resl <- c(resl, list(projected = emn, current = emat.norm,
                       deltaE = deltaE, deltaT = deltaT, ko = full.ko, mult = mult,
                       kCells = kCells,kSites = kSites,cell.use.list=cell.use.list))
  
  return(resl)
}
