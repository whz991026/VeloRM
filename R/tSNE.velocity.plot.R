##' Joint t-SNE visualization of the velocities by joint t-SNE embedding of both current and extraploated cell positions
##'
##' @param vel velocity result
##' @param cell.colors named color vector for the cells
##' @param scale whether to rescale current/projected
##' @param do.par whether to reset par (default=T)
##' @param delta.norm whether to renormalize velocities following PCA projection
##' @param nPcs number of PCs onto which project the velocities
##' @param norm.nPcs number of PCs to use for velocity normalization
##' @param perplexity perplexity parameter to use in joint t-SNE calculation
##' @param show.grid.flow whether grid flow pattern should be drawn
##' @param grid.n number of grid points along each axis
##' @param grid.sd standard deviation of each grid point (used to determine the averaging radius for each grid point)
##' @param min.grid.cell.mass minimal number of cells around a grid point required for the grid point to show up
##' @param pcount pseudocount
##' @param verbose whether to show messages
##' @param min.arrow.median.ratio minimal ratio of arrow length (to the median arrow length) below which the arrows are not drawn (default=1/10)
##' @param max.arrow.quantile max arrow quantile that's used for arrow size calculation (default=0.9)
##' @param arrow.scale scaling factor for the arrows
##' @param arrow.lwd arrow line width
##' @param xlab x axis label
##' @param ylab y axis label
##' @param n.cores number of cores to use
##' @param size.norm whether to re-normalize current and projected cell sizes
##' @param ... extra parameters are passed to plot() routine.
##'
##' @import Rtsne
##'
##' @return invisible list containing embedding positions of current state (current.emb) and extrapolated states (projected.emb)
##' @export
tSNE.velocity.plot <- function(vel,cell.colors=NULL,scale='log',do.par=T,
                               delta.norm=TRUE, nPcs=15, norm.nPcs=nPcs*10,
                               perplexity=ncol(vel$current)/3, show.grid.flow=FALSE,
                               grid.n=20, grid.sd=NULL, min.grid.cell.mass=1,
                               pcount=0.1, verbose=TRUE, min.arrow.median.ratio=1/10,
                               max.arrow.quantile=0.9, arrow.scale=1, arrow.lwd=1,
                               xlab="", ylab="", n.cores=1, size.norm=TRUE,
                               ...) {
  x0 <- vel$current;
  x1 <- vel$projected;
  if(is.null(cell.colors)) { cell.colors <- ac(rep(1,ncol(x0)),alpha=0.3); names(cell.colors) <- colnames(x0) }

  if(size.norm) {
    # rescale to the same size
    cat("rescaling ... ")
    sz <- Matrix::colSums(x0);
    x0 <- t(t(x0)/sz)*mean(sz)
    x1 <- t(t(x1)/Matrix::colSums(x1))*mean(sz);
  }
  # transform
  if(scale=='log') {
    cat("log ... ")
    x0.log <- log2(x0+pcount)
    x1.log <- log2(x1+pcount)
  } else if(scale=='sqrt') {
    cat("sqrt ... ")
    x0.log <- sqrt(x0)
    x1.log <- sqrt(x1)
  } else { # linear
    cat("linear ... ")
    x0.log <- x0
    x1.log <- x1
  }
  if(!is.null(nPcs)) { # reduce using PCA first
    cat("pca ... ")
    cent <- rowMeans(x0.log);
    epc <- pcaMethods::pca(t(x0.log-cent),center=F,nPcs=ifelse(is.na(norm.nPcs),nPcs,norm.nPcs))
    x0.log <- epc@scores;
    x1.log <- t(x1.log - cent) %*% epc@loadings
    if(delta.norm) {
      # normalize velocities ...?
      cat("delta norm ... ")
      delta.pcs <- x1.log-x0.log;
      if(!is.na(norm.nPcs)){
        delta.pcs <- delta.pcs/mean(sqrt(rowSums(delta.pcs^2))) # ? unsure about this (cell-wise L2 norm)
      }
      x1.log <- x0.log+delta.pcs;
    }
    # drop extra Pcs
    x0.log <- t(x0.log[,1:nPcs])
    x1.log <- t(x1.log[,1:nPcs])
  }
  cat("tSNE ...")
  emb <- Rtsne::Rtsne(t(cbind(as.matrix(x0.log),as.matrix(x1.log))), num_threads=n.cores, perplexity=perplexity, verbose=verbose)$Y;
  x0.emb <- emb[1:ncol(x0.log),]
  x1.emb <- emb[-(1:ncol(x0.log)),]
  rownames(x0.emb) <- rownames(x1.emb) <- colnames(x0.log);

  cat("delta norm ... ") # again, somewhat unsure about this part
  delta.emb <- x1.emb - x0.emb;
  asize <- rowSums(delta.emb^2);
  no.arrow <- asize<= median(asize)*min.arrow.median.ratio;
  # restrict size top top 90% quantile
  delta.emb <- delta.emb/asize * pmin(asize,2*quantile(asize,p=max.arrow.quantile))*arrow.scale;
  delta.emb[no.arrow,] <- 0;
  cat("done\n")
  if(do.par) par(mfrow=c(1,1), mar = c(3.5,3.5,2.5,1.5), mgp = c(2,0.65,0), cex = 0.85);
  plot(x0.emb,bg=cell.colors[rownames(x0.emb)],pch=21,col=ac(1,alpha=0.3),xlab=ylab,ylab=xlab, ... ); box();


  if(show.grid.flow) { # show grid summary of the arrows
    # arrow estimates for each cell
    ars <- data.frame(x0.emb[,1],x0.emb[,2],x0.emb[,1]+delta.emb[,1],x0.emb[,2]+delta.emb[,2])
    colnames(ars) <- c('x0','y0','x1','y1')
    arsd <- data.frame(xd=ars$x1-ars$x0,yd=ars$y1-ars$y0)

    # set up a grid
    cat("grid estimates ... ")
    rx <- range(c(range(ars$x0),range(ars$x1)))
    ry <- range(c(range(ars$y0),range(ars$y1)))
    gx <- seq(rx[1],rx[2],length.out=grid.n)
    gy <- seq(ry[1],ry[2],length.out=grid.n)

    # for each grid point calculate Gaussian-weighted delta average
    if(is.null(grid.sd)) {
      grid.sd <- sqrt((gx[2]-gx[1])^2 + (gy[2]-gy[1])^2)/2
    }
    ginfo <- lapply(gx,function(x) {
      # cell distances (rows-cells,columsn - grid points)
      cd <- sqrt(outer(x0.emb[,2],-gy,'+')^2 + (x-x0.emb[,1])^2)
      cw <- dnorm(cd,sd=grid.sd)
      # calculate x and y delta expectations
      gw <- Matrix::colSums(cw)
      cws <- pmax(1,Matrix::colSums(cw));
      gxd <- Matrix::colSums(cw*arsd$xd)/cws
      gyd <- Matrix::colSums(cw*arsd$yd)/cws

      vg <- gw>=min.grid.cell.mass
      if(any(vg)) {
        suppressWarnings(arrows(rep(x,sum(vg)),gy[vg],x+gxd[vg],gy[vg]+gyd[vg],length=0.05,lwd=arrow.lwd))
      }
      points(rep(x,length(gy)),gy,pch='.')
    })
    cat("done\n")
  } else {
    # draw individual arrows
    suppressWarnings(arrows(x0.emb[,1],x0.emb[,2],x0.emb[,1]+delta.emb[,1],x0.emb[,2]+delta.emb[,2],length=0.05,lwd=arrow.lwd))
  }

  return(invisible(list(current.emb=x0.emb,projected.emb=x1.emb+delta.emb)))
}
