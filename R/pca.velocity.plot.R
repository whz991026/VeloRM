##' PCA-based visualization of the velocities
##'
##' @param vel velocity estimation (gene-relative or global)
##' @param nPcs number of successive PCs to visualize
##' @param cell.colors a named vector of cell colors for visualization
##' @param scale scale to use for expression state transform (default: 'log', other possible values are 'sqrt','linear')
##' @param plot.cols number of columns into which to arrange the plots
##' @param norm.nPcs optional total number of PCs to use for velocity magnitude normalization
##' @param do.par whether to set up graphical parameters of a plot
##' @param pc.multipliers an optional vector multipliers for the cell PC scores (useful for reorienting the PCs)
##' @param show.grid.flow whether a grid flow should be shown
##' @param grid.n number of grid points (on each axis)
##' @param grid.sd standard deviation of the grid
##' @param arrow.scale scale multiplier for the velocity estimates
##' @param min.grid.cell.mass minimum cellular mass
##' @param min.arrow.size minimum size of an arrow to show
##' @param pcount pseudocount
##' @param arrow.lwd thickness of arrows to plot
##' @param size.norm whether to rescale current and projected states by cell size (default=FALSE)
##' @param return.details whether to return detailed output
##' @param plot.grid.points whether to show dots at every grid point
##' @param fixed.arrow.length whether to use fixed-size arrow
##' @param max.grid.arrow.length limit to the size of the arrows that could be shown (when fixed.arrow.length=FALSE)
##' @param n.cores number of cores to use in the calculations
##' @param ... extra parameters are passed to plot() function
##'
##' @importFrom pcaMethods pca
##' @importFrom graphics arrows grid points
##' @importFrom stats dnorm
##'
##' @return If return.details=F, returns invisible list containing PCA info (epc) and projection of velocities onto the PCs (delta.pcs). If return.details=T, returns an extended list that can be passed into p1 app for velocity visualization.
##' @export
pca.velocity.plot <- function(vel,nPcs=4,cell.colors=NULL,scale='log',
                              plot.cols=min(3,nPcs-1),norm.nPcs=NA,do.par=T,
                              pc.multipliers=NULL, show.grid.flow=FALSE, grid.n=20,
                              grid.sd=NULL, arrow.scale=1, min.grid.cell.mass=1,
                              min.arrow.size=NULL, pcount=1, arrow.lwd=1,
                              size.norm=FALSE, return.details=FALSE,
                              plot.grid.points=FALSE, fixed.arrow.length=FALSE,
                              max.grid.arrow.length=NULL, n.cores=1, ...) {
  x0 <- vel$current;
  x1 <- vel$projected;
  if(is.null(cell.colors)) { cell.colors <- ac(rep(1,ncol(x0)),alpha=0.3); names(cell.colors) <- colnames(x0) }
  # rescale to the same size
  if(size.norm) {
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

  cat("pca ... ")
  cent <- rowMeans(x0.log);
  epc <- pcaMethods::pca(t(x0.log-cent),center=F,nPcs=ifelse(is.na(norm.nPcs),nPcs,norm.nPcs))

  if(!is.null(pc.multipliers)) { # apply multipliers (used for flipping the direction of PCs in the plots)
    if(length(pc.multipliers)!=nPcs) stop("pc.multipliers must be a vector equal in length to the number of PCs")
    cat("pc multipliers ... ")
    epc@loadings <- t(t(epc@loadings)*pc.multipliers)
    epc@scores <- scale(epc@completeObs,scale=F,center=T) %*% epc@loadings;
  }

  x1.scores <- t(x1.log - cent) %*% epc@loadings

  # normalize velocities ...?
  cat("delta norm ... ")
  delta.pcs <- as.matrix(x1.scores-epc@scores)
  if(!is.na(norm.nPcs)) {
    delta.pcs <- delta.pcs/mean(sqrt(rowSums(delta.pcs^2))) # suggested by Gioele, unsure about this ....
  }

  # browser()
  # z <- as.matrix(t(x1.log-x0.log)) %*% epc@loadings
  #
  # hist(apply(vel$deltaE,2,mean))
  # summary(apply(vel$deltaE,2,mean))
  # hist(apply(as.matrix(x1.log-x0.log),2,mean))
  # summary(apply(as.matrix(x1.log-x0.log),2,mean))
  #
  # z <- t(as.matrix(vel$deltaE)[rownames(epc@loadings),]) %*%  epc@loadings
  # str(z)
  # str(delta.pcs)
  # cn <- 'L6'
  # cn <- 'O15'
  # z <- (x1.log-x0.log)[,cn] * epc@loadings[,3]
  # z2 <- vel$deltaE[names(z),cn] * epc@loadings[,3]
  # sort(z,d=T)[1:20]
  # sum(z)
  # delta.pcs[cn,3]
  # summary(delta.pcs[,3])
  # str(epc@loadings[,3])
  # sort(delta.pcs[,3],d=T)[1:10]
  # z <- rowMeans(x1.log-x0.log) * epc@loadings[,3]
  #summary(z)
  #sort(z,d=T)[1:10]

  delta.pcs <- delta.pcs *arrow.scale;
  cat("done\n")
  if(do.par) par(mfrow=c(ceiling((nPcs-1)/plot.cols),plot.cols), mar = c(3.5,3.5,2.5,1.5), mgp = c(2,0.65,0), cex = 0.85);
  vinfo <- lapply(1:(nPcs-1),function(i) {
    pos <- epc@scores[,c((i-1)+1,(i-1)+2)];
    #ppos <- x1.scores[,c((i-1)+1,(i-1)+2)];
    ppos <- pos+delta.pcs[,c((i-1)+1,(i-1)+2)];
    plot(pos,bg=cell.colors[rownames(pos)],pch=21,col=ac(1,alpha=0.3),lwd=0.5,xlab=paste("PC",(i-1)+1),ylab=paste("PC",(i-1)+2),axes=T,main=paste('PC',(i-1)+1,' vs. PC',(i-1)+2,sep=''),  ...); box();

    if(show.grid.flow) { # show grid summary of the arrows
      # arrow estimates for each cell
      ars <- data.frame(pos[,1],pos[,2],ppos[,1],ppos[,2])
      colnames(ars) <- c('x0','y0','x1','y1')
      arsd <- data.frame(xd=ars$x1-ars$x0,yd=ars$y1-ars$y0)
      rownames(ars) <- rownames(arsd) <- rownames(pos);

      # set up a grid
      rx <- range(c(range(ars$x0),range(ars$x1)))
      ry <- range(c(range(ars$y0),range(ars$y1)))
      gx <- seq(rx[1],rx[2],length.out=grid.n)
      gy <- seq(ry[1],ry[2],length.out=grid.n)

      # for each grid point calculate Gaussian-weighted delta average
      if(is.null(grid.sd)) {
        grid.sd <- sqrt((gx[2]-gx[1])^2 + (gy[2]-gy[1])^2)/2
        cat("grid.sd=",grid.sd," ")
      }

      if(is.null(min.arrow.size)) {
        min.arrow.size <- sqrt((gx[2]-gx[1])^2 + (gy[2]-gy[1])^2)*1e-2;
        cat("min.arrow.size=",min.arrow.size," ")
      }

      if(is.null(max.grid.arrow.length)) {
        max.grid.arrow.length <- sqrt(sum((par('pin')/c(length(gx),length(gy)))^2))*0.25
        cat("max.grid.arrow.length=",max.grid.arrow.length," ")
      }


      garrows <- do.call(rbind,lapply(gx,function(x) {
        # cell distances (rows:cells, columns: grid points)
        cd <- sqrt(outer(pos[,2],-gy,'+')^2 + (x-pos[,1])^2)
        cw <- dnorm(cd,sd=grid.sd)
        # calculate x and y delta expectations
        gw <- Matrix::colSums(cw)
        cws <- pmax(1,Matrix::colSums(cw));
        gxd <- Matrix::colSums(cw*arsd$xd)/cws
        gyd <- Matrix::colSums(cw*arsd$yd)/cws

        al <- sqrt(gxd^2+gyd^2);
        vg <- gw>=min.grid.cell.mass & al>=min.arrow.size

        cbind(rep(x,sum(vg)),gy[vg],x+gxd[vg],gy[vg]+gyd[vg])
      }))
      colnames(garrows) <- c('x0','y0','x1','y1')

      # plot
      if(fixed.arrow.length) {
        suppressWarnings(arrows(garrows[,1],garrows[,2],garrows[,3],garrows[,4],length=0.05,lwd=arrow.lwd))
      } else {
        alen <- pmin(max.grid.arrow.length,sqrt( ((garrows[,3]-garrows[,1]) * par('pin')[1] / diff(par('usr')[c(1,2)]) )^2 + ((garrows[,4]-garrows[,2])*par('pin')[2] / diff(par('usr')[c(3,4)]) )^2))
        # can't specify different arrow lengths in one shot :(
        suppressWarnings(lapply(1:nrow(garrows),function(i) arrows(garrows[i,1],garrows[i,2],garrows[i,3],garrows[i,4],length=alen[i],lwd=arrow.lwd)))
      }
      if(plot.grid.points) points(rep(gx,each=length(gy)),rep(gy,length(gx)),pch='.',cex=1e-1,col=ac(1,alpha=0.4))

      if(return.details) { # for the p1 app
        # calculate expression shift
        cat("expression shifts .")
        # for individual cells
        es <- as.matrix(epc@loadings[,c((i-1)+1,(i-1)+2)] %*% t(delta.pcs[,c((i-1)+1,(i-1)+2)]))

        cat(".");
        gs <- epc@loadings[,c((i-1)+1,(i-1)+2)] %*% rbind(garrows[,3]-garrows[,1],garrows[,4]-garrows[,2])

        # note: here we're using deltaE vector, which may be normalized a bit differently from the $current/$projectted that was used above
        nd <- as.matrix(vel$deltaE)
        if(scale=='log') {
          nd <- (log10(abs(nd)+1)*sign(nd))
        } else if(scale=='sqrt') {
          nd <- (sqrt(abs(nd))*sign(nd))
        }
        cat(".");
        # velocity for the grid (weight-averaged velocity vectors)


        gv <- do.call(cbind,parallel::mclapply(gx,function(x) {
          # cell distances (rows:cells, columns: grid points)
          cd <- sqrt(outer(pos[,2],-gy,'+')^2 + (x-pos[,1])^2)
          cw <- dnorm(cd,sd=grid.sd)
          # calculate x and y delta expectations
          gw <- Matrix::colSums(cw)
          cws <- pmax(1,Matrix::colSums(cw));
          cw <- t(t(cw)/cws)
          gxd <- Matrix::colSums(cw*arsd$xd)
          gyd <- Matrix::colSums(cw*arsd$yd)
          al <- sqrt(gxd^2+gyd^2);
          vg <- gw>=min.grid.cell.mass & al>=min.arrow.size
          if(any(vg)) {
            z <- nd %*% cw[,vg]
          } else { NULL }
        },mc.cores=n.cores,mc.preschedule=T))
        cat(". done\n")

        return(invisible(list(garrows=garrows,arrows=as.matrix(ars),vel=nd,eshifts=es,gvel=gv,geshifts=gs,scale=scale,emb=pos,epc=epc)))
      }

    } else {
      # draw individual arrows
      grid();
      suppressWarnings(arrows(pos[,1],pos[,2],ppos[,1],ppos[,2],length=0.05,lwd=arrow.lwd))
    }
  })
  cat("done\n")
  if(return.details) { return(vinfo) }
  return(invisible(list(epc=epc,delta.pcs=delta.pcs)))

}
