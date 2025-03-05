##' Visualize RNA velocities on an existing embedding using Euclidean-based transition probability matrix within the kNN graph.
##'
##'  based on Euclidean distance of the extrapolated cell to others
##' The direction of the arrow is towards n closest neighbors. The magnitude of the arrow is determined by the cosine projection of the velocity on to the chosen direction
##' n=1 will only show arrows for cells that end up being projected closer to some other cell than to the original position
##' n=k (k>1) will show an average direction
##' Given an expression distance between cells d, and ratio of extrapolated to current expression distances between cells f, the transition probability is calculated as exp(- (d*(f^beta))^2/(2*sigma^2) )
##'
##' @param emb embedding to be used for projection
##' @param current current matrix
##' @param projected projected matrix
##' @param deltaE deltaE matrix
##' @param n neighborhood size (default=30)
##' @param embedding.knn pre-calculated kNN
##' @param cell.colors named color vector for cell plotting
##' @param sigma sigma to use in calculating transition probability from the eucledian distance (estimated automatically by default)
##' @param beta beta parameter used in calculation of transition probability (by default=1)
##' @param arrow.scale additional scaling factor for the arrows (default=1)
##' @param scale scale to use in calculating distances (default: 'log', also supported 'sqrt'
##' @param nPcs number of PCs to project the cells onto (to perform distance calculations in lower dimensions), default=NA which turns off PCA dimensional reduction
##' @param arrow.lwd arrow line width
##' @param xlab x axis label
##' @param ylab y axis label
##' @param control.for.neighborhood.density compensate for cell density variations in the embedding (default: TRUE)
##' @param cell.dist - optional custom distance (must include all of the cells that are intersecting between emb and vel)
##' @param cell.border.alpha trasparency parameter to apply when showing cell colors
##' @param n.cores number of cores to use in calculations
##' @param show.cell whether to show detailed velocity estimates for a specified cell
##' @param show.grid.flow whether to show grid velocity summary
##' @param grid.n number of grid points along each axis
##' @param grid.sd standard deviation (in embedding coordinate space) used to determine the weighting of individual cells around each grid point
##' @param min.grid.cell.mass minimal cell "mass" (weighted number of cells) around each grid point required for it to show up
##' @param min.arrow.size minimal arrow size
##' @param max.grid.arrow.length minimal arrow size
##' @param fixed.arrow.length whether to use fixed arrow width (default=FALSE)
##' @param plot.grid.points whether to mark all grid points with dots (even if they don't have valid velocities)
##' @param return.details whether to return detailed output (which can be passed to p1 app for visualization)
##' @param expression.scaling whether to scale the velocity length by the projection of velocity onto the expected expression change (based on the transition probability matrix)
##' @param point.size size of the point
##' @param arrow_size size of arrow
##' @param rm_uniform mode1: tp minus the uniform, mode2: tp minus the uniform and discard the negative,
##' mode3: do nothing for the tp
##' @param ... extra parameters are passed to the plot() function

##' @importFrom cluster pam
##' @importFrom graphics xspline
##' @importFrom abind abind
##' @import igraph
##'
##' @return transition probability matrix
##' @export
show.velocity.on.embedding.eu <- function(emb,current,projected,deltaE,n=30,embedding.knn=TRUE,
                                          cell.colors=NULL, sigma=NA, beta=1,
                                          arrow.scale=1, scale='log', nPcs=NA,
                                          arrow.lwd=1, xlab="", ylab="",
                                          show.grid.flow=FALSE,
                                          grid.n=20, grid.sd=NULL, min.grid.cell.mass=1,
                                          min.arrow.size=NULL, show.cell=NULL,
                                          max.grid.arrow.length=NULL,
                                          fixed.arrow.length=FALSE,
                                          plot.grid.points=FALSE,
                                          control.for.neighborhood.density=TRUE, 
                                          cell.dist=NULL, return.details=FALSE,
                                          expression.scaling=FALSE,
                                          cell.border.alpha=0.5, n.cores=1,
                                          point.size=3,arrow_size=0.3,
                                          rm_uniform="mode1",...) {
  options(warn = -1)
  em <- current; emn <- projected;
  if(is.null(cell.colors)) { cell.colors <- ac(rep(1,ncol(em)),alpha=0.3); names(cell.colors) <- colnames(em) }
  celcol <- 'white'
  if(is.null(show.cell)) { celcol <- cell.colors[rownames(emb)] }
  
  

  ccells <- intersect(rownames(emb),colnames(em));
  emn <- emn[,ccells]; em <- em[,ccells]; emb <- emb[ccells,]
  
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
    theme(legend.position="none", plot.title=element_text(size=12,hjust=0.5))

  if(scale=='log') {
    cat("log scale ... ")
    em <- log10(em+1); emn <- log10(emn+1);
  } else if(scale=='sqrt') {
    cat("sqrt scale ... ")
    em <- sqrt(em); emn <- sqrt(emn);
  }else { # linear
    cat("linear scale ... ")
  }

  if(!is.na(nPcs)) { # run PCA reduction on the em
    cat("reducing to",nPcs,"PCs ... ")
    epc.center <- rowMeans(em);
    epc <- pcaMethods::pca(t(em-epc.center),center=F,nPcs=nPcs);
    em <- t(epc@scores)
    emn <- t(t(emn - epc.center) %*% epc@loadings)
  }

  cat("distance ... ")
  cc <- colEuclid(as.matrix(em),as.matrix(emn),n.cores)
  cc0 <- colEuclid(as.matrix(em),as.matrix(em),n.cores)
  cd <- (cc0-cc); # reduction in the Euclidean distance

  if(n>nrow(cc)) { n <- nrow(cc) }

  # pick reasonable sigma and beta if they weren't provided
  if(is.na(sigma) | is.na(beta)) { mcd <- mean(as.matrix(abs(cd)/cc)) }
  # TODO: adaptive methods for signal
  if(is.na(sigma)) { sigma <- mcd/10 }
  if(is.na(beta)) { beta <- mcd/20 }
  cat("sigma=",round(sigma,3)," beta=",round(beta,3)," transition probs ... ")

  # exp(- (d*(f^beta))^2/(2*sigma^2) )
  f <- (cc/cc0)^beta; diag(f) <- 1;
  tp <- exp(- ((cc0*f)^2) / (2*sigma^2))
  np <- exp(- ((cc0)^2) / (2*sigma^2))




  if(n<nrow(emb)) {
    if(!is.null(cell.dist)) {
      cat("kNN on provided distance ... ")
      if(!all(labels(cell.dist)==colnames(em))) {
        cat("matching cells between cell.dist and emat/nmat ... ")
        cell.dist <- as.matrix(cell.dist)
        cn <- colnames(em)
        cell.dist <- as.dist(cell.dist[cn,cn]);
      }
      cell.knn <- balanced_knn(as.matrix(cell.dist),k=n,maxl=nrow(emb),
                               method=attributes(cell.dist)$method,n.cores=n.cores)
      diag(cell.knn) <- 1;
    } else {
      if(embedding.knn) {
        cat("embedding kNN ... ")
        # define kNNs based on the embedding (L2 distance)
        cell.knn <- balanced_knn(as.matrix(dist((emb))),k=n,maxl=nrow(emb),
                                 method='dist',n.cores=n.cores)
        #diag(cell.knn) <- 0; # disallow self-transitions?
        diag(cell.knn) <- 1;
      } else {
        cat("expression kNN ... ")
        # define kNN based on the correlation distance in high-d
        cell.knn <- balanced_knn((cor(as.matrix(em))),k=n,maxl=ncol(em),
                                 method='cor',n.cores=n.cores)
        diag(cell.knn) <- 1;
      }
    }
    tp <- as.matrix(tp)*cell.knn;
    np <- as.matrix(np)*cell.knn;
  }

  # estimate density of the neighborhood
  tp <- t(t(tp)/Matrix::colSums(tp))
  np <- t(t(np)/Matrix::colSums(np))
  #diag(tp) <- diag(np) <- 0;
  #tp.nd <- colSums(tp); np.nd <- colSums(np);
  #tp <- tp * np.nd;

  if(control.for.neighborhood.density) {
    np.f <- Matrix::diag(np);
    tp <- tp*(np.f)
    np <- np*(np.f)
  }

  #diag(tp) <- diag(np) <- 0;
  # normalize
  tp <- t(t(tp)/Matrix::colSums(tp))
  np <- t(t(np)/Matrix::colSums(np))


  # normalize transition probabilities
  rownames(tp) <- colnames(tp) <- rownames(np) <- colnames(np) <- colnames(em);
  cat("done\n")
 
  # arrow estimates for each cell
  cat("calculating arrows ... ")
  arsd <- data.frame(t(embArrows(emb,tp,rm_uniform,arrow.scale,n.cores)))
  rownames(arsd) <- rownames(emb);
  
  if(expression.scaling) {
    tpb <- tp>0; tpb <- t(t(tpb)/colSums(tpb));
    es <- as.matrix(em %*% tp) -as.matrix(em %*% as.matrix(tpb));
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

