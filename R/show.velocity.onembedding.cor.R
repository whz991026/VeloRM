##' Visualize RNA velocities on an existing embedding using correlation-based transition probability matrix within the kNN graph
##'
##' @param emb embedding onto which to project the velocities; The dimensions of coordinates should be on the order of 10x10 for the default values to make sense.
##' @param current current matrix
##' @param projected projected matrix
##' @param deltaE deltaE matrix
##' @param n neighborhood size (default=100 cells)
##' @param cell.colors name vector of cell colors
##' @param corr.sigma sigma parameter used to translate velocity-(expression delta) correlation into a transition probability
##' @param show.grid.flow whether to show grid velocity summary
##' @param grid.n number of grid points along each axis
##' @param grid.sd standard deviation (in embedding coordinate space) used to determine the weighting of individual cells around each grid point
##' @param min.grid.cell.mass minimal cell "mass" (weighted number of cells) around each grid point required for it to show up
##' @param min.arrow.size minimal arrow size
##' @param arrow.scale arrow scale multiplier
##' @param max.grid.arrow.length minimal arrow size
##' @param fixed.arrow.length whether to use fixed arrow width (default=FALSE)
##' @param plot.grid.points whether to mark all grid points with dots (even if they don't have valid velocities)
##' @param scale velocity scale to use (default: 'log', other values: 'sqrt','rank','linear')
##' @param nPcs number of PCs to use for velocity regularization (default NA, turns off regularization)
##' @param arrow.lwd arrow width (under fixed.arrow.length=T)
##' @param xlab x axis label
##' @param ylab y axls label
##' @param n.cores number of cores to use
##' @param show.cell whether to show detailed velocity estimates for a specified cell
##' @param cell.border.alpha transparency for the cell border
##' @param cc velocity-(exprssion delta) correlation matrix (can be passed back from previous results, as $cc) to save calculation time when replotting the same velocity estimates on the same embedding with different parameters
##' @param return.details whether to return detailed output (which can be passed to p1 app for visualization)
##' @param expression.scaling whether to scale the velocity length by the projection of velocity onto the expected expression change (based on the transition probability matrix)
##' @param point.size size of the point
##' @param arrow_size size of arrow
##' @param ... extra parameters passed to plot() function
##'
##' @importFrom stats rbinom
##'
##' @return if return.details=F, returns invisible list containing transition probability matrix ($tp) and the velocity-(expression delta) correlation matrix ($cc). If return.details=T, returns a more extended list that can be passed as veloinfo to pagoda2::p2.make.pagoda1.app() for visualization
##' @export
show.velocity.on.embedding.cor <- function(emb,current,projected,deltaE,n=100,cell.colors=NULL,
                                           corr.sigma=0.05, show.grid.flow=FALSE,
                                           grid.n=20, grid.sd=NULL, min.grid.cell.mass=1,
                                           min.arrow.size=NULL, arrow.scale=1,
                                           max.grid.arrow.length=NULL,
                                           fixed.arrow.length=FALSE,
                                           plot.grid.points=FALSE, scale='log',
                                           nPcs=NA,  arrow.lwd=1, xlab="", ylab="",
                                           n.cores=1, 
                                           show.cell=NULL, cell.border.alpha=0.5,
                                           cc=NULL, return.details=FALSE,
                                           expression.scaling=FALSE,point.size=3,arrow_size=0.3,  ...) {
  options(warn = -1)
  randomize <- FALSE;
  celcol <- 'white'
  if(is.null(show.cell)) { celcol <- cell.colors[rownames(emb)]; names(celcol) <- rownames(emb) }
  
  
  em <- as.matrix(current);
  ccells <- intersect(rownames(emb),colnames(em));
  em <- em[,ccells]; emb <- emb[ccells,]
  nd <- as.matrix(deltaE[,ccells])
  csites <- intersect(rownames(em),rownames(nd));
  nd <- nd[csites,]; em <- em[csites,]
  
  point.data.frame <- as.data.frame(emb[ccells,1])
  point.data.frame$PC1 <- emb[ccells,1]
  point.data.frame$PC2 <- emb[ccells,2]
  point.data.frame$color <- celcol[ccells]
  point.data.frame <- point.data.frame[,c("PC1","PC2","color")]
  plot.ggplot <- ggplot() +
    geom_point(data=point.data.frame, aes(x=.data$PC1, y=.data$PC2, 
                    fill=.data$color), alpha=cell.border.alpha,size=point.size,shape=21) + 
    labs(x="PC1", y="PC2", 
         title="PC1 vs PC2") + 
    theme(plot.title=element_text(size=12,hjust=0.5),legend.position="none")
  
  #browser()
  if(randomize) {
    # randomize cell and sign for each site
    nd <- t(apply(nd,1,function(x) (rbinom(length(x),1,0.5)*2-1)*abs(sample(x))))
  }
  #vg <- rownames(em) %in% rownames(r)


  if(is.null(cc)) {
    # cosine projections
    cat("delta projections ... ")

    if(scale=='log') {
      cat("log ")
      cc <- colDeltaCor (em, nd, "log", n.cores);
    } else if(scale=='sqrt') {
      cat("sqrt ")
      cc <- colDeltaCor (em, nd, "sqrt", n.cores);
    } else if(scale=='rank') {
      cat("rank ")
      cc <- colDeltaCor (em, nd, "rank", n.cores);
    } else { # linear
      cat("linear ")
      cc <- colDeltaCor (em, nd, "linear", n.cores);
    }
    colnames(cc) <- rownames(cc) <- colnames(em)
    diag(cc) <- 0;
  }
  if(n<nrow(emb)){
    cat("knn ... ")
    if(n>nrow(cc)) { n <- nrow(cc) }
    # TODO: add kNN based on high-dimensional correlation or Euclidean distances
    # define kNNs based on the embedding (L2 distance)
    emb.knn <- balanced_knn(as.matrix(dist(emb)),k=n,maxl=nrow(emb),'dist',n.cores)
    diag(emb.knn) <- 1
    # caluclate transition probabilities (from col to row)
  } else{
    emb.knn <- 1
  }
  
  cat("transition probs ... ")
  tp <- exp(cc/corr.sigma)*emb.knn
  #diag(tp) <- 0; #  should we allow the self-corelation for scaling?
  tp <- t(t(tp)/Matrix::colSums(tp)); # tp shows transition from a given column cell to different row cells
  tp <- as(tp,'dgCMatrix')
  cat("done\n")
  
  # arrow estimates for each cell
  cat("calculating arrows ... ")
  arsd <- data.frame(t(embArrows(emb,tp,arrow.scale,n.cores)))
  rownames(arsd) <- rownames(emb);

  if(expression.scaling) {
    tpb <- tp>0; tpb <- t(t(tpb)/colSums(tpb));
    es <- as.matrix(em %*% tp) -as.matrix(em %*% as.matrix(tpb));
    # project velocity onto expression shift
    #pm <- as.matrix(t(deltaE)/sqrt(colSums(deltaE*deltaE)))[colnames(es),] * (t(es)/sqrt(colSums(es*es)))
    #pl <- pmax(0,apply(pm,1,sum))
    pl <- pmin(1,pmax(0,apply(as.matrix(deltaE[,colnames(es)]) * es, 2, sum)/sqrt(colSums(es*es))))


    arsd <- arsd * pl;
  }


  ars <- data.frame(cbind(emb,emb+arsd));
  colnames(ars) <- c('x0','y0','x1','y1')
  colnames(arsd) <- c('xd','yd')
  rownames(ars) <- rownames(emb);
  cat("done\n")

  if(!is.null(show.cell)) {
    i <- match(show.cell,rownames(emb));
    if(is.na(i)) stop(paste('specified cell',i,'is not in the embedding'))
    
    ars_point <- ars[show.cell,]
    plot.ggplot <- plot.ggplot + geom_point(data= ars_point, 
                      aes(x=.data$x0, y=.data$y0),color="black",size=4)+
    geom_segment(data=ars_point, aes(x=.data$x0, y=.data$y0, xend=.data$x1, yend=.data$y1), 
                    arrow=arrow(length=unit(arrow_size, "inches")), size=arrow.lwd)+
      theme(legend.position="none") 
  } else{
    if(show.grid.flow) { # show grid summary of the arrows
      
      # set up a grid
      cat("grid estimates ... ")
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
        cd <- sqrt(outer(emb[,2],-gy,'+')^2 + (x-emb[,1])^2)
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
        plot.ggplot <- plot.ggplot + 
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
        
        for (j in 1:nrow(garrows)) {
          plot.ggplot <- plot.ggplot + annotate("segment",
                                                x=garrows$x0[j], y=garrows$y0[j],
                                                xend=garrows$x1[j], yend=garrows$y1[j],
                                                arrow=custom_arrow(garrows$arrow_length[j]),
                                                size=arrow.lwd)
        }
        if(plot.grid.points) {
          plot.ggplot <- plot.ggplot + geom_point(data=garrows[j,], 
                                                  aes(x=.data$x0, y=.data$y0),color="black", alpha=0.7, shape = 17,size=2)
        } 
      }
      
      cat("done\n")
      
      if(return.details) { # for the p1 app
        # calculate expression shift
        cat("expression shifts .")
        # for individual cells
        
        scale.int <- switch(scale,'log'=2,'sqrt'=3,1)
        #es <- expectedExpressionShift(e=as.matrix(em),tp=tp,scale=scale.int,nthreads=n.cores); colnames(es) <- colnames(em); rownames(es) <- rownames(em);
        if(!expression.scaling) { #otherwise it has already been calculated
          tpb <- tp>0; tpb <- t(t(tpb)/colSums(tpb));
          #es <- expectedExpressionShift(e=as.matrix(em %*% as.matrix(tpb)),tp=tp,scale=scale.int,nthreads=n.cores); colnames(es) <- colnames(em); rownames(es) <- rownames(em);
          es <- as.matrix(em %*% tp) -as.matrix(em %*% as.matrix(tpb));
        }
        cat(".");
        # for the grid
        gs <- do.call(cbind,parallel::mclapply(gx,function(x) {
          # cell distances (rows:cells, columns: grid points)
          cd <- sqrt(outer(emb[,2],-gy,'+')^2 + (x-emb[,1])^2)
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
            z <- es %*% cw[,vg]
          } else { NULL }
        },mc.cores=n.cores,mc.preschedule=T))
        
        if(scale=='log') {
          nd <- (log10(abs(nd)+1)*sign(nd))
        } else if(scale=='sqrt') {
          nd <- (sqrt(abs(nd))*sign(nd))
        }
        cat(".");
        # velocity for the grid
        gv <- do.call(cbind,parallel::mclapply(gx,function(x) {
          # cell distances (rows:cells, columns: grid points)
          cd <- sqrt(outer(emb[,2],-gy,'+')^2 + (x-emb[,1])^2)
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
        
        
        return(list(plot.ggplot,invisible(list(tp=tp,cc=cc,garrows=garrows,arrows=as.matrix(ars),vel=nd,eshifts=es,gvel=gv,geshifts=gs,scale=scale))))
      }
      
    } else { 
      if(fixed.arrow.length) {
        plot.ggplot <- plot.ggplot +
          geom_segment(data=ars, aes(x=.data$x0, y=.data$y0, xend=.data$x1, yend=.data$y1), 
                       arrow=arrow(length=unit(arrow_size, "inches")), size=arrow.lwd)+
          theme(legend.position="none") 
      } else {
        alen <- pmin(0.05,sqrt(((ars[,3]-ars[,1]) * par('pin')[1] / diff(par('usr')[c(1,2)]))^2 + 
                                 ((ars[,4]-ars[,2]) * par('pin')[2] / diff(par('usr')[c(3,4)]))^2))
        ars <- data.frame(ars, arrow_length=alen)
        custom_arrow <- function(length) {
          grid::arrow(length = unit(length, "inches"))
        }
        
        
        for (j in 1:nrow(ars)) {
          plot.ggplot <- plot.ggplot + annotate("segment",
                                                x=ars$x0[j], y=ars$y0[j],
                                                xend=ars$x1[j], yend=ars$y1[j],
                                                arrow=custom_arrow(ars$arrow_length[j]),
                                                size=arrow.lwd)
        }
      }
      
      
    }
  }
  
 
  
  return(list(plot.ggplot,invisible(list(tp=tp,cc=cc))))
}

