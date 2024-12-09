##' PCA-based visualization of the velocities
##'
##' @param current current matrix
##' @param projected projected matrix
##' @param deltaE deltaE matrix
##' @param nPcs number of successive PCs to visualize
##' @param cell.colors a named vector of cell colors for visualization
##' @param scale scale to use for expression state transform (default: 'log', other possible values are 'sqrt','linear')
##' @param plot.cols number of columns into which to arrange the plots
##' @param norm.nPcs optional total number of PCs to use for velocity magnitude normalization
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
##' @param point.size size of the point
##' @param cell.border.alpha transparency for the cell border
##' @param arrow_size size of arrow
##' @param ... extra parameters are passed to plot() function
##'
##' @importFrom pcaMethods pca
##' @importFrom graphics arrows grid points
##' @importFrom stats dnorm
##' @importFrom grid arrow
##'
##' @return If return.details=F, returns invisible list containing PCA info (epc) and projection of velocities onto the PCs (delta.pcs). If return.details=T, returns an extended list that can be passed into p1 app for velocity visualization.
##' @export
pca.velocity.plot <- function(current,projected,deltaE,nPcs=4,cell.colors=NULL,scale='log',
                              plot.cols=min(3,nPcs-1),norm.nPcs=NA,
                              pc.multipliers=NULL, show.grid.flow=FALSE, grid.n=20,
                              grid.sd=NULL, arrow.scale=1, min.grid.cell.mass=1,
                              min.arrow.size=NULL, pcount=1, arrow.lwd=1,
                              size.norm=FALSE, return.details=FALSE,
                              plot.grid.points=FALSE, fixed.arrow.length=FALSE,
                              max.grid.arrow.length=NULL, n.cores=1,point.size=3,
                              cell.border.alpha=0.5,arrow_size=0.3,...) {
  x0 <- current
  x1 <- projected
  if(is.null(cell.colors)) { cell.colors <- ac(rep(1,ncol(x0)),alpha=0.3); names(cell.colors) <- colnames(x0) }
  # rescale to the same size
  if(size.norm) {
    cat("rescaling ... ")
    sz <- Matrix::colSums(x0)
    x0 <- t(t(x0)/sz)*mean(sz)
    x1 <- t(t(x1)/Matrix::colSums(x1))*mean(sz)
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
  cent <- rowMeans(x0.log)
  epc <- pcaMethods::pca(t(x0.log-cent),center=F,nPcs=ifelse(is.na(norm.nPcs),nPcs,norm.nPcs))
  
  if(!is.null(pc.multipliers)) { # apply multipliers (used for flipping the direction of PCs in the plots)
    if(length(pc.multipliers)!=nPcs) stop("pc.multipliers must be a vector equal in length to the number of PCs")
    cat("pc multipliers ... ")
    epc@loadings <- t(t(epc@loadings)*pc.multipliers)
    epc@scores <- scale(epc@completeObs,scale=F,center=T) %*% epc@loadings
  }
  
  x1.scores <- t(x1.log - cent) %*% epc@loadings
  
  # normalize velocities ...?
  cat("delta norm ... ")
  delta.pcs <- as.matrix(x1.scores-epc@scores)
  if(!is.na(norm.nPcs)) {
    delta.pcs <- delta.pcs/mean(sqrt(rowSums(delta.pcs^2))) # suggested by Gioele, unsure about this ....
  }
  
  delta.pcs <- delta.pcs * arrow.scale
  cat("done\n")
  
  vinfo <- lapply(1:(nPcs-1),function(i) {
    pos <- epc@scores[,c((i-1)+1,(i-1)+2)]
    ppos <- pos+delta.pcs[,c((i-1)+1,(i-1)+2)]
    point_data.frame <- as.data.frame(pos[,1])
    point_data.frame$PC1 <- pos[,1]
    point_data.frame$PC2 <- pos[,2]
    point_data.frame$color <- cell.colors[rownames(pos)]
    point_data.frame <- point_data.frame[,c("PC1","PC2","color")]
    
    if(show.grid.flow) { # show grid summary of the arrows
      # arrow estimates for each cell
      ars <- data.frame(pos[,1],pos[,2],ppos[,1],ppos[,2])
      colnames(ars) <- c('x0','y0','x1','y1')
      arsd <- data.frame(xd=ars$x1-ars$x0,yd=ars$y1-ars$y0)
      rownames(ars) <- rownames(arsd) <- rownames(pos)
      
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
        min.arrow.size <- sqrt((gx[2]-gx[1])^2 + (gy[2]-gy[1])^2)*1e-2
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
        cws <- pmax(1,Matrix::colSums(cw))
        gxd <- Matrix::colSums(cw*arsd$xd)/cws
        gyd <- Matrix::colSums(cw*arsd$yd)/cws
        
        al <- sqrt(gxd^2+gyd^2)
        vg <- gw>=min.grid.cell.mass & al>=min.arrow.size
        
        cbind(rep(x,sum(vg)),gy[vg],x+gxd[vg],gy[vg]+gyd[vg])
      }))
      colnames(garrows) <- c('x0','y0','x1','y1')
      
      # plot
      if(fixed.arrow.length) {
        plot.ggplot <- ggplot() +
          geom_point(data=point_data.frame, aes(x=.data$PC1, y=.data$PC2, 
                                    fill=.data$color), alpha=cell.border.alpha,size=point.size, shape = 21) + 
          labs(x=paste("PC",(i-1)+1), y=paste("PC",(i-1)+2), 
               title=paste('PC',(i-1)+1,' vs PC',(i-1)+2,sep='')) + 
          theme(legend.position="none", plot.title=element_text(size=12,hjust=0.5)) + 
          geom_segment(data=garrows, aes(x=.data$x0, y=.data$y0, xend=.data$x1, yend=.data$y1), 
                       arrow=arrow(length=unit(arrow_size, "inches")), size=arrow.lwd)
      } else {
        alen <- pmin(max.grid.arrow.length,sqrt(((garrows[,3]-garrows[,1]) * par('pin')[1] / diff(par('usr')[c(1,2)]))^2 + 
                                                  ((garrows[,4]-garrows[,2]) * par('pin')[2] / diff(par('usr')[c(3,4)]))^2))
        garrows <- data.frame(garrows, arrow_length=alen)
        custom_arrow <- function(length) {
          grid::arrow(length = unit(length, "inches"))
        }
        plot.ggplot <- ggplot() +
          geom_point(data=point_data.frame, aes(x=.data$PC1, y=.data$PC2,
                         fill=.data$color), alpha=cell.border.alpha,size=point.size, shape = 21) + 
          labs(x=paste("PC",(i-1)+1), y=paste("PC",(i-1)+2), 
               title=paste('PC',(i-1)+1,' vs PC',(i-1)+2,sep='')) + 
          theme(legend.position="none", plot.title=element_text(size=12,hjust=0.5))
        
        for (j in 1:nrow(garrows)) {
          plot.ggplot <- plot.ggplot + annotate("segment",
                                                x=garrows$x0[j], y=garrows$y0[j],
                                                xend=garrows$x1[j], yend=garrows$y1[j],
                                                arrow=custom_arrow(garrows$arrow_length[j]),
                                                size=arrow.lwd)
          if(plot.grid.points) {
            plot.ggplot <- plot.ggplot + geom_point(data=garrows[j,], 
              aes(x=.data$x0, y=.data$y0),color="black", alpha=0.7, shape = 17,size=2)
          } 
        }
      }
      
      
      if(return.details) { # for the p1 app
        # calculate expression shift
        cat("expression shifts .")
        # for individual cells
        es <- as.matrix(epc@loadings[,c((i-1)+1,(i-1)+2)] %*% t(delta.pcs[,c((i-1)+1,(i-1)+2)]))
        
        cat(".")
        gs <- epc@loadings[,c((i-1)+1,(i-1)+2)] %*% rbind(garrows[,3]-garrows[,1],garrows[,4]-garrows[,2])
        
        # note: here we're using deltaE vector, which may be normalized a bit differently from the $current/$projectted that was used above
        nd <- as.matrix(deltaE)
        if(scale=='log') {
          nd <- (log10(abs(nd)+1)*sign(nd))
        } else if(scale=='sqrt') {
          nd <- (sqrt(abs(nd))*sign(nd))
        }
        cat(".")
        # velocity for the grid (weight-averaged velocity vectors)
        gv <- do.call(cbind,parallel::mclapply(gx,function(x) {
          # cell distances (rows:cells, columns: grid points)
          cd <- sqrt(outer(pos[,2],-gy,'+')^2 + (x-pos[,1])^2)
          cw <- dnorm(cd,sd=grid.sd)
          # calculate x and y delta expectations
          gw <- Matrix::colSums(cw)
          cws <- pmax(1,Matrix::colSums(cw))
          cw <- t(t(cw)/cws)
          gxd <- Matrix::colSums(cw*arsd$xd)
          gyd <- Matrix::colSums(cw*arsd$yd)
          al <- sqrt(gxd^2+gyd^2)
          vg <- gw>=min.grid.cell.mass & al>=min.arrow.size
          if(any(vg)) {
            z <- nd %*% cw[,vg]
          } else { NULL }
        },mc.cores=n.cores,mc.preschedule=T))
        cat(". done\n")
        
        return(list(plot.ggplot,invisible(list(garrows=garrows,arrows=as.matrix(ars),vel=nd,eshifts=es,
                                   gvel=gv,geshifts=gs,scale=scale,emb=pos,epc=epc))))
      }else{
        plot.ggplot
      }
      
    } else {
      garrows <- as.data.frame(pos[,1])
      garrows$x0 <- pos[,1]
      garrows$y0 <- pos[,2]
      garrows$x1 <- ppos[,1]
      garrows$y1 <- ppos[,2]
      garrows <- garrows[,c("x0","y0","x1","y1")]
      if (fixed.arrow.length==TRUE){
        # draw individual arrows
        plot.ggplot <- ggplot() +
          geom_point(data=point_data.frame, aes(x=.data$PC1, y=.data$PC2, fill=.data$color),
                     alpha=cell.border.alpha,size=point.size, shape = 21) + 
          labs(x=paste("PC",(i-1)+1), y=paste("PC",(i-1)+2), 
               title=paste('PC',(i-1)+1,' vs PC',(i-1)+2,sep='')) + 
          theme(legend.position="none", plot.title=element_text(size=12,hjust=0.5)) 
        plot.ggplot <- plot.ggplot  + 
          geom_segment(data=garrows, aes(x=.data$x0, y=.data$y0, xend=.data$x1, yend=.data$y1), 
                       arrow=arrow(length=unit(arrow_size, "inches")), size=arrow.lwd)
        plot.ggplot
      } else{
        alen <- pmin(0.05,sqrt(((garrows[,3]-garrows[,1]) * par('pin')[1] / diff(par('usr')[c(1,2)]))^2 + 
                                                  ((garrows[,4]-garrows[,2]) * par('pin')[2] / diff(par('usr')[c(3,4)]))^2))
        garrows <- data.frame(garrows, arrow_length=alen)
        custom_arrow <- function(length) {
          grid::arrow(length = unit(length, "inches"))
        }
        plot.ggplot <- ggplot() +
          geom_point(data=point_data.frame, aes(x=.data$PC1, y=.data$PC2,
                            fill=.data$color), alpha=cell.border.alpha,size=point.size, shape = 21) + 
          labs(x=paste("PC",(i-1)+1), y=paste("PC",(i-1)+2), 
               title=paste('PC',(i-1)+1,' vs PC',(i-1)+2,sep='')) + 
          theme(legend.position="none", plot.title=element_text(size=12,hjust=0.5))
        
        for (j in 1:nrow(garrows)) {
          plot.ggplot <- plot.ggplot + annotate("segment",
                                                x=garrows$x0[j], y=garrows$y0[j],
                                                xend=garrows$x1[j], yend=garrows$y1[j],
                                                arrow=custom_arrow(garrows$arrow_length[j]),
                                                size=arrow.lwd)
        }
        plot.ggplot
      }
      
    }
  })
  cat("done\n")
  if(return.details==TRUE){
    plot_list <- lapply(vinfo,function(df) df[1][[1]])
    details_list <-  lapply(vinfo,function(df) df[2])
    
    plot <- do.call(cowplot::plot_grid, plot_list)
    return(list(plot,details_list,invisible(list(epc=epc,delta.pcs=delta.pcs))))
  } else{
    plot <- do.call(cowplot::plot_grid, vinfo)
    return(list(plot,invisible(list(epc=epc,delta.pcs=delta.pcs))))
  }
  
}
