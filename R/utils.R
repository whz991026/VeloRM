# calculates the difference in the number of counts based on the library size, renormalizes
# note: also introduces centripetal velocity
#
# @param em normalized expression matrix
# @param cellSize the cell size of the spliced read
# @param deltae the delta spliced read count matrix
# @param delta  time to project forward
# @param mult library scaling factor (1e6 in case of FPM)
#
# @return the projected cell spliced read counts
t.get.projected.cell2 <- function(em,cellSize,deltae,mult=1e3,delta=1) {
  rz <- matrix(0,nrow=nrow(em),ncol=ncol(em)); colnames(rz) <- colnames(em);
  rownames(rz) <- rownames(em)
  gn <- intersect(rownames(deltae),rownames(rz))
  rz[match(gn,rownames(rz)),colnames(deltae)] <- deltae[gn,];
  # translate fpm delta into the number of molecules based on the current cell size
  rz <- t(t(rz)*cellSize)
  emm <- t(t(em)*cellSize)
  emn <- emm + rz*delta;
  emn[emn<0] <- 0;
  newCellSize <- (cellSize+Matrix::colSums(emn-emm)/mult)
  emn <- t(t(emn)/newCellSize)

  #emn <- t(t(emn)/Matrix::colSums(emn)*Matrix::colSums(em))
  emn
}

sizeFactor <- function(data, narm=FALSE) {
  # make the elements no smaller than 0
  data <- as.matrix(data)
  
  
  # first log
  log_data <- log(data)
  log_data [is.infinite(log_data)] <- NA
  log_mean <- rowMeans(log_data, na.rm=narm)
  log_s <- log_data-log_mean
  
  # then exp
  s_size <- exp(apply(log_s,2,function(x)median(x,na.rm=TRUE)))
  return(s_size)
}



# estimate projected delta under model 1
#
# @param em  normalized expression matrix
# @param nm normalized nascent matrix
# @param gamma inferred degradation coefficients
# @param offset inferred offset (assumed to be zero by default)
# @param delta time to project forward
#
# @return the delta spliced read count matrix
t.get.projected.delta.model1 <- function(em,nm,gamma,offset=rep(0,length(gamma)),delta=0.5) {
  # adjust rownames
  gn <- intersect(names(gamma),rownames(em));
  if(is.null(names(offset))) { names(offset) <- names(gamma); }
  em <- em[gn,]; nm <- nm[gn,]; gamma <- gamma[gn]; offset <- offset[gn];
  
  y <- (nm-offset - gamma*em)*delta
}


# estimate projected delta under model 2
#
# @param em  normalized expression matrix
# @param nm normalized nascent matrix
# @param gamma inferred degradation coefficients
# @param offset inferred offset (assumed to be zero by default)
# @param delta time to project forward
#
# @return the delta spliced read count matrix
t.get.projected.delta.model2 <- function(em,nm,gamma,offset=rep(0,length(gamma)),delta=0.5) {
  # adjust rownames
  gn <- intersect(names(gamma),rownames(em));
  if(is.null(names(offset))) { names(offset) <- names(gamma); }
  em <- em[gn,]; nm <- nm[gn,]; gamma <- gamma[gn]; offset <- offset[gn];
  # time effect constant
  egt <- exp(-gamma*delta);
  y <- nm-offset; y[y<0] <- 0; # zero out entries with a negative n levels after offset adjustment
  em*egt + (1-egt)*y/gamma  - em
}









# balanced_knn calculate the balanced k neighborhood
#
# @param d the distance matrix or the correlation matrix (should be square matrix)
# @param k The number of neighborhood
# @param maxl The maximal number of times that a cell be chosen as neighborhood
# @param method whether the matrix is correlation matrix or distance matrix c("dist","cor")
# @param n.cores The number of cores to use for parallel computation, default is 1
# @importFrom parallel mclapply
# @importFrom Matrix sparseMatrix
#
# @return knn matrix
balanced_knn <- function(d, k, maxl, method="dist", n.cores=1) {
  l <- rep(0, ncol(d)) # reciprocal neighbor count
  dsi <- matrix(0, nrow = nrow(d), ncol = ncol(d)) # will store column sort indices of the d matrix

  if(ncol(d) != nrow(d)){
    stop("the input d should be a square matrix")
  }

  if(ncol(d) < k){
    stop("the k should be smaller or equal to the number of the columns in the d matrix")
  }

  sort_indices <- function(i) {
    if(method == "dist") {
      order(d[, i])
    } else if(method == "cor") {
      order(d[, i], decreasing = TRUE)
    } else {
      stop("method parameter needs to be set as dist or cor")
    }
  }

  sorted_indices <- parallel::mclapply(1:ncol(d), sort_indices, mc.cores = n.cores)
  for (i in 1:ncol(d)) {
    si <- sorted_indices[[i]]
    dsi[, i] <- si
    l[si[1:k]] <- l[si[1:k]] + 1 # record the number of the calling
  }

  lsi <- order(l, decreasing = TRUE) # greedy order of column considerations
  l <- rep(0, ncol(d)) # reset so that l can be used to keep track of reciprocal counts in the kNN being constructed

  rowind <- numeric(k * ncol(d))
  vals <- rep(1, k * ncol(d)) # values

  for (i in 1:ncol(d)) {
    el <- lsi[i] # The cell that is called the most times needs to be seen first.
    si <- dsi[, el]
    p <- 1
    j <- 1
    while (j <= nrow(dsi) && p <= k) {
      m <- si[j]
      if (el == m) {
        j <- j + 1
        next
      } # dont record or count self-relationships
      if (l[m] >= maxl) {
        j <- j + 1
        next
      } # element has already maxed out its neighbors
      rowind[(el - 1) * k + p] <- m  # Adjusted indices to start from 1
      l[m] <- l[m] + 1
      p <- p + 1 # record neighbor

      j <- j + 1
    }
    if (j > nrow(dsi) && p < k) {
      rowind[(el - 1) * k + p] <- el
      p <- p + 1 # fill in last element(s) with self-identity if there were not enough spare elements
      while (p <= k) {
        warning(paste("unable to find unfilled neighbors for i =", i, ", el =", el, ", p =", p))
        rowind[(el - 1) * k + p] <- el
        p <- p + 1
      }
    }
  }


  colptr <- seq(0, ncol(d) * k, by = k)
  knn <- Matrix::sparseMatrix(i = rowind, p = colptr, x = rep(1, length(rowind)),
                              dims = c(nrow(d), ncol(d))) # kNN matrix
  # the knn only contains the number of cells that are similar to the cell, without itself.

  return(knn)
}



# Calculate the correlation matrix between the delta and the difference of the
# current spliced read counts
#
# @param current The current spliced read counts
# @param delta The delta calculated from the velocity model
# @param method_delta_scale scale ("log","sqrt","rank","linear")
# @param n.cores number of cores
# @importFrom stats cor
#
# @return  pearson correlation matrix
colDeltaCor <- function(current, delta, method_delta_scale, n.cores = 1) {
  options(warn = -1)
  if (method_delta_scale == "log") {
    transfor_function <- function(df) {
      log10(abs(df) + 1) * sign(df)
    }
  } else if (method_delta_scale == "sqrt") {
    transfor_function <- function(df) {
      (sqrt(abs(df)) * sign(df))
    }
  } else if (method_delta_scale == "rank") {
    transfor_function <- function(df) {
      apply(df, 2, rank)
    }
  } else if (method_delta_scale == "linear") {
    transfor_function <- function(df) {
      df
    }
  } else {
    stop("the method_delta_scale parameter should be c(log,sqrt,rank,linear)")
  }

  delta_transfor <- transfor_function(delta)

  if (method_delta_scale == "rank" || method_delta_scale == "linear") {
    current <- transfor_function(current)
    data.frame <- as.data.frame(parallel::mclapply(1:dim(delta_transfor)[2], function(ii) {
      current_difference <- current - current[, ii]
      unlist(parallel::mclapply(1:dim(delta_transfor)[2], function(jj) {
        cor(delta_transfor[, ii], current_difference[, jj])
      }))
    }, mc.cores = n.cores))
  } else {
    data.frame <- as.data.frame(parallel::mclapply(1:dim(delta_transfor)[2], function(ii) {
      current_difference <- transfor_function(current - current[, ii])
      unlist(parallel::mclapply(1:dim(delta_transfor)[2], function(jj) {
        cor(delta_transfor[, ii], current_difference[, jj])
      }))
    }, mc.cores = n.cores))
  }

  data.frame[is.na(data.frame)] <- 0
  data.frame <- as.matrix(data.frame)
  rownames(data.frame) <- colnames(data.frame) <- colnames(delta)
  data.frame
}

# Calculate the correlation matrix between the delta and the difference of the
# current spliced read counts
#
# @param e a matrix
# @param p a matrix
# @param n.cores number of cores
#
# @return  Euclid distance matrix
colEuclid <- function(e, p, n.cores = 1) {
  options(warn = -1)

  euclidean_distances <- function(x, y) {
    sqrt(sum((x - y)^2))
  }

  data.frame <- parallel::mclapply(seq_len(ncol(p)), function(j) {
    unlist(parallel::mclapply(seq_len(ncol(e)), function(k) {
      euclidean_distances(p[, j], e[, k])
    }, mc.cores = n.cores))
  }, mc.cores = n.cores)

  data.frame <- as.data.frame(data.frame)
  rownames(data.frame) <- colnames(data.frame) <- colnames(e)
  data.frame
}







# calculate the delta embedding values
#
# @param emb the embedding matrix
# @param tp transition probabilities matrix
# @param rm_uniform mode1: tp minus the uniform, mode2: tp minus the uniform and discard the negative,
# mode3: do nothing for the tp
# @param arrow.scale the scale of the arrow
# @param n.cores number of cores
#
# @return the delta embedding values
embArrows <- function(emb, tp,rm_uniform="mode1", arrow.scale = 1.0,n.cores=1) {
  dm <- matrix(0, nrow = ncol(emb), ncol = nrow(emb))
  tpb <- as.matrix(tp)  # Convert sparse matrix to matrix
  tpb[tpb != 0] <- 1
  tprs <- colSums(tpb != 0)
  tpb <- t(t(tpb) / tprs)  # Normalize

  temb <- t(emb)

  # Define a function to calculate dm[, i]
  calc_dm <- function(i) {
    di <- (temb - temb[, i])  # Coordinate difference
    di <- apply(di, 2, function(df) df / sum(abs(df)))  # Normalize, scale
    di[, i] <- 0  # No distance to itself
    di <- di * arrow.scale
    if(rm_uniform=="mode1"){
      ds <- di %*% tp[, i] - di %*% tpb[, i]
    } else if (rm_uniform=="mode2"){
      tp[, i] <- tp[, i] - tpb[, i]
      tp[, i] <- pmax(tp[, i],0)
      ds <- di %*% tp[, i] 
    } else if (rm_uniform=="mode3"){
      ds <- di %*% tp[, i] 
    } else{
      stop("rm_uniform should be mode1, mode2, and mode3")
    }
    return(ds)
  }

  # Use mclapply to parallelize the loop
  dm_list <- parallel::mclapply(1:nrow(emb), calc_dm, mc.cores = n.cores)

  # Combine results into matrix
  dm <- do.call(cbind, dm_list)

  return(dm)
}
#' Filter sites by requiring minimum average expression within at least one of the provided cell clusters
#'
#' @param emat spliced (exonic) count matrix
#' @param clusters named cell factor defining clusters
#' @param min.max.cluster.average required minimum average expression count (no normalization is perfomed)
#' @return filtered emat matrix
#' @export
filter_sites_by_cluster_expression <- function(emat,clusters,min.max.cluster.average=0.1) {
  if(!any(colnames(emat) %in% names(clusters))) stop("provided clusters do not cover any of the emat cells!")
  vc <- intersect(colnames(emat),names(clusters))
  cl.emax <- apply(do.call(cbind,tapply(vc,as.factor(clusters[vc]),function(ii) Matrix::rowMeans(emat[,ii]))),1,max)
  vi <- cl.emax>min.max.cluster.average;
  emat[vi,]
}

# quick self-naming vector routine
sn <- function(x) { names(x) <- x; x}

#' adjust colors, while keeping the vector names
#'
#' @param x color vector
#' @param alpha transparenscy value (passed to adjustcolors as alpha.f)
#' @param ... parameters passsed to adjustcolor
#' @importFrom grDevices adjustcolor colorRampPalette
#' @importFrom stats cor
#' @export
ac <- function(x, alpha=1, ...) { y <- adjustcolor(x, alpha.f=alpha, ...);
          names(y) <- names(x); return(y)}

# quick function to map value vector to colors
val2col <- function(x,gradientPalette=NULL,zlim=NULL,gradient.range.quantile=0.95) {
  if(all(sign(x)>=0)) {
    if(is.null(gradientPalette)) {
      gradientPalette <- colorRampPalette(c('gray90','red'), space = "Lab")(1024)
    }
    if(is.null(zlim)) {
      zlim <- as.numeric(quantile(na.omit(x),p=c(1-gradient.range.quantile,gradient.range.quantile)))
      if(diff(zlim)==0) {
        zlim <- as.numeric(range(na.omit(x)))
      }
    }
    x[x<zlim[1]] <- zlim[1]; x[x>zlim[2]] <- zlim[2];
    x <- (x-zlim[1])/(zlim[2]-zlim[1])

  } else {
    if(is.null(gradientPalette)) {
      gradientPalette <- colorRampPalette(c("blue", "grey90", "red"), space = "Lab")(1024)
    }
    if(is.null(zlim)) {
      zlim <- c(-1,1)*as.numeric(quantile(na.omit(abs(x)),p=gradient.range.quantile))
      if(diff(zlim)==0) {
        zlim <- c(-1,1)*as.numeric(na.omit(max(abs(x))))
      }
    }
    x[x<zlim[1]] <- zlim[1]; x[x>zlim[2]] <- zlim[2];
    x <- (x-zlim[1])/(zlim[2]-zlim[1])

  }

  gp <- gradientPalette[x*(length(gradientPalette)-1)+1]
  if(!is.null(names(x))) { names(gp) <- names(x) }
  gp
}



#' the diffusion density plot
#'
#' @param emb embedding to be used for projection
#' @param tp transition probability 
#' @param forward forward or backward
#' @param iteration the number of step
#' @param scale scale the density or not "linear" or "log"
#' @param log_unit the hyper parameter of the log scale
#' @param low_color the hyper-parameter for show site low color
#' @param high_color the hyper-parameter for show site high color
#' @param point.size size of the point
#' @param cell.border.alpha transparency for the cell border
#' 
#' @import ggplot2
#' @importFrom stats density 
#' @importFrom rlang .data
#' @return a data frame of the plot
#' @export
#'
diffusion_density_plot <- function(emb,tp,forward=TRUE,iteration=1000,scale="linear",
                                   log_unit=1e-5,low_color="#FFCC15",high_color="#9933FF",
                                   point.size=2,cell.border.alpha=1){
  ccells <- intersect(rownames(emb),rownames(tp));
  emb <- emb[ccells,];tp <- tp[ccells,ccells]
  
  cp <- Matrix::Diagonal(ncol(tp)); # cell position probabilities
  rownames(cp) <- colnames(cp) <- rownames(tp);
  if(forward==TRUE){
    ttp <- t(tp);
  } else{
    ttp <- tp/rowSums(tp)
  }
  
  
  # run diffusion steps to figure out end positions
  cat("simulating diffusion ... ")
  for(i in 1:iteration) {
    cp <- cp %*% ttp;
    #cp[cp<1e-5] <- 0;
  }
  cat("done\n");
  
  
  if (scale=="log"){
    cp_summary <- colSums(cp)/sum(colSums(cp))
    cp_summary <- log(cp_summary+log_unit)
  } 
  if (scale=="linear"){
    cp_summary <- colSums(cp)/sum(colSums(cp))
    cp_summary <- cp_summary
  } 
  
  
  
  data.frame.plot <- as.data.frame(emb)
  data.frame.plot$density <- cp_summary
  colnames(data.frame.plot)<- c("dim1","dim2","density")
  data.frame.plot$dim1 <- as.numeric(data.frame.plot$dim1)
  data.frame.plot$dim2 <- as.numeric(data.frame.plot$dim2)
  data.frame.plot$density <- as.numeric(data.frame.plot$density)
  
  polot_ggplot2 <- ggplot(data.frame.plot)+
    aes(x = .data$dim1, y = .data$dim2, color = .data$density) + 
    geom_point(size=point.size,alpha=cell.border.alpha) +
    scale_color_gradient(low = low_color, high = high_color)
  
  if (forward==TRUE){
    polot_ggplot2 <- polot_ggplot2 +
      labs(title="forward")+ 
      theme(plot.title=element_text(size=12,hjust=0.5))
  } else{
    polot_ggplot2 <- polot_ggplot2 +
      labs(title="backward")+ 
      theme(plot.title=element_text(size=12,hjust=0.5))
  }
  
  
  polot_ggplot2
}


