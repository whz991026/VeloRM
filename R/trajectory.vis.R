#' trajectory visualization 
#' @param emb embedding to be used for projection
#' @param tp transition probability 
#' @param vel velocity result
#' @param emb_local_influence whether consider the emb influence
#' @param local_influence_prob the proportion of the influence
#' @param sigma_d the hyper-parameter of local influence
#' @param cell.colors name vector of cell colors
#' @param show.cell.arrows show detailed velocity projection for the specified cell
#' @param show.cell.trajectories show trajectories for a specified cell
#' @param show.trajectories show top median diffusion trajectories
#' @param show.all.trajectories show all diffusion paths (messy)
#' @param arrow.lwd arrow line width
#' @param xlab x axis label
#' @param ylab y axls label
#' @param do.par whether to reset plotting parameters
#' @param n.cores number of cores to use in calculations
#' @param cell.border.alpha transparency for the cell border
#' @param show.cell.diffusion.posterior show diffusion posterior of a given cell
#' @param diffusion.steps number of diffusion steps to take forward (default=10)
#' @param arrow.scale additional scaling factor for the arrows (default=1)
#' @param show.grid.flow whether to show grid velocity summary
#' @param ntop.trajectories number of top trajectories to trace back for a given cell (when show.trajectories=TRUE)
#' @param trajectory.spline.shape shape parameter for smoothing median cell trajectories (default=1)
#' @param n.trajectory.clusters number of trajectory clusters to show median paths for (when show.trajectories=TRUE)
#' @param ... extra parameters passed to plot() function
#'
#' @export
#'
trajectory.vis <- function( emb,tp,vel,emb_local_influence=FALSE,local_influence_prob=0.2,
                            sigma_d=NA,arrow.scale=1,show.cell.arrows=NULL,
                            show.cell.trajectories=NULL,cell.colors=NULL,
                            show.trajectories=FALSE,arrow.lwd=1,
                            show.all.trajectories=FALSE,n.cores=1,
                            show.cell.diffusion.posterior=NULL, xlab="", ylab="",
                            show.grid.flow=FALSE, diffusion.steps=10,
                            cell.border.alpha=0.5,do.par=TRUE,
                            ntop.trajectories=1, trajectory.spline.shape=1,
                            n.trajectory.clusters=10,...){
  if(do.par) par(mfrow=c(1,1), mar = c(3.5,3.5,2.5,1.5), mgp = c(2,0.65,0), cex = 0.85);
  if (is.null(show.cell.diffusion.posterior)){
    celcol <- cell.colors[rownames(emb)] 
  } else{
    celcol <- "white"
  }
  
  plot(emb,bg=celcol,pch=21,col=ac(1,alpha=cell.border.alpha), xlab=xlab, ylab=ylab, ...);
  
  
  em <- vel$current; 
  ccells <- intersect(intersect(rownames(emb),colnames(em)),rownames(tp));
  em <- em[,ccells]; emb <- emb[ccells,];tp <- tp[ccells,ccells]
  
  if(emb_local_influence==TRUE){
    cat("re-calculate the transition probability ...")
    # calculate the euclidean distance according to the emb
    cc <- colEuclid(t(as.matrix(emb)),t(as.matrix(emb)),n.cores)
    cc <- as.matrix(cc)
    if(is.na(sigma_d)){
      sigma_d <- sqrt(mean(cc,na.rm=TRUE))
    }
    sigma_w <- sigma_d/2
    tp_new <- local_influence_prob* exp(- ((cc)^2) / (2*sigma_w^2))+
      (1-local_influence_prob)*tp*exp(- ((cc)^2) / (2*sigma_d^2))
    rownames(tp_new) <- colnames(tp_new) <- rownames(tp)
    tp <- tp_new
    tp <- tp/rowSums(tp)
    cat("done\n")
  }
  
  
  
  if(!is.null(show.cell.diffusion.posterior)) {
    i <- match(show.cell.diffusion.posterior,rownames(emb));
    if(is.na(i)) stop(paste('specified cell',i,'is not in the embedding'))
    # run diffusion
    cp <- Matrix::Diagonal(ncol(tp)); # cell position probabilities
    rownames(cp) <- colnames(cp) <- rownames(tp);
    ttp <- t(tp);
    
    # run diffusion steps to figure out end positions
    cat("simulating diffusion ... ")
    for(i in 1:diffusion.steps) {
      cp <- cp %*% ttp;
      #cp[cp<1e-5] <- 0;
    }
    cat("done\n");
    # plot
    points(emb,pch=19,col=ac(val2col(cp[show.cell.diffusion.posterior,rownames(emb)],gradient.range.quantile=1),alpha=0.5))
    points(emb[show.cell.diffusion.posterior,1],emb[show.cell.diffusion.posterior,2],pch=3,cex=1,col=1)
  } else if(!is.null(show.cell.arrows)) {
    i <- match(show.cell.arrows,rownames(emb));
    if(is.na(i)) stop(paste('specified cell',i,'is not in the embedding'))
    # plot transition prob for a given cell
    points(emb,pch=19,col=ac(val2col(tp[rownames(emb),show.cell.arrows],gradient.range.quantile=1),alpha=0.5))
    points(emb[i,1],emb[i,2],pch=3,cex=1,col=1)
    di <- t(t(emb)-emb[i,])
    di <- di/sqrt(Matrix::rowSums(di^2))*arrow.scale; di[i,] <- 0;
    dir <- Matrix::colSums(di*tp[,i])
    dic <- Matrix::colSums(di*(tp[,i]>0)/sum(tp[,i]>0)); # relative to expected kNN center
    dia <- dir-dic;
    #browser()
    suppressWarnings(arrows(emb[colnames(em)[i],1],emb[colnames(em)[i],2],emb[colnames(em)[i],1]+dic[1],emb[colnames(em)[i],2]+dic[2],length=0.05,lwd=1,col='blue'))
    suppressWarnings(arrows(emb[colnames(em)[i],1],emb[colnames(em)[i],2],emb[colnames(em)[i],1]+dir[1],emb[colnames(em)[i],2]+dir[2],length=0.05,lwd=1,col='red'))
    suppressWarnings(arrows(emb[colnames(em)[i],1]+dic[1],emb[colnames(em)[i],2]+dic[2],emb[colnames(em)[i],1]+dir[1],emb[colnames(em)[i],2]+dir[2],length=0.05,lwd=1,lty=1,col='grey50'))
    suppressWarnings(arrows(emb[colnames(em)[i],1],emb[colnames(em)[i],2],emb[colnames(em)[i],1]+dia[1],emb[colnames(em)[i],2]+dia[2],length=0.05,lwd=1,col='black'))
  } else if(show.trajectories) { # show diffusion paths
    cp <- Matrix::Diagonal(ncol(tp)); # cell position probabilities
    rownames(cp) <- colnames(cp) <- rownames(tp);
    
    #cpt <- as.array(cp);
    cpl <- list(); cpl[[1]] <- cp;
    #cp <- as.matrix(cp); tp <- as.matrix(tp)
    
    #ep <- as.array(emb)
    ttp <- t(tp);
    
    # run diffusion steps to figure out end positions
    cat("simulating diffusion ... ")
    for(i in 1:diffusion.steps) {
      cp <- cp %*% ttp;
      #cp[cp<1e-5] <- 0;
      #cpt <- abind(cpt,cp,along=3)
      cpl[[i+1]] <- cp;
      # clean up to zero out all but top n cells
      #cp <- t(apply(cp,1,function(x) { x[x<sort(x,decreasing=TRUE)[10]] <- 0; x }))
      #diag(cp) <- 0; #  prohibit the cell from returning to itself
      #cp <- cp/Matrix::rowSums(cp)
    }
    
    
    #cpt <- abind(lapply(cpl,as.matrix),along=3)
    
    # calculate probabilistic trajectories to the final ntop points
    
    # rank final points by probability
    cpo <- t(apply(-cp,1,order))
    
    # graph-based walkback approach
    
    # construct a walkback graph
    trp <- as(ttp,'TsparseMatrix')
    cat("constructing path graph ... ")
    x <- do.call(rbind,lapply(1:(diffusion.steps+1),function(i) {
      cbind(i=trp@i+1 + (i-1)*nrow(cp), # current time step
            j=trp@j+1 + (i)*nrow(cp))
    }))
    x <- x[x[,2]<=nrow(cp)*(diffusion.steps+1),]
    x <- Matrix::spMatrix(nrow(cp)*(diffusion.steps+1),nrow(cp)*(diffusion.steps+1),i=x[,1],j=x[,2],x=rep(1,nrow(x)))
    g <- igraph::graph.adjacency(x,mode='directed')
    rm(x); gc();
    
    # find topn trajectories for each cell
    cat("tracing shortest trajectories ... ")
    sps <- parallel::mclapply(1:nrow(cp),function(celli) {
      top.desti <- order(cp[celli,],decreasing=TRUE)[1:ntop.trajectories]
      # calculate cell-specific weights
      cw <- unlist(lapply(cpl,function(d) as.numeric(trp@x*(d[celli,trp@i+1]))))
      cw <- cw[1:igraph::ecount(g)] # trim extra edges
      # convert into penalty scores
      cw <- -log(cw)
      sp <- igraph::shortest_paths(g,from=celli,to=nrow(cp)*(diffusion.steps-1)+top.desti,weights=cw,mode='out')
      # remove time offset on the path nodes
      sp <- lapply(sp$vpath,function(x) { y <- (as.integer(x) %% nrow(cp)); y[y==0] <- nrow(cp); y});
      names(sp) <- rownames(cp)[top.desti]
      sp
    },mc.cores=n.cores)
    
    
    
    # cluster paths
    cat("clustering ... ")
    all.cells <- 1:nrow(cp)
    #spuci <- do.call(cbind,lapply(sps,function(x) all.cells %in% x[[1]]))
    # filter out empty paths
    sps <- lapply(sps,function(y) y[unlist(lapply(y,function(x) length(unique(x))))>1])
    spuci <- do.call(cbind,lapply(sps,function(y) do.call(cbind,lapply(y,function(x) all.cells %in% x))))
    usps <- unlist(sps,recursive=F); # will be used in looking up median trajectories in plotting
    
    
    spuci.dist <- as.matrix(dist(t(spuci),method = 'manhattan'))
    spuci.pam <- cluster::pam(spuci.dist,n.trajectory.clusters)
    cat("done.\n")
    
    # bezier
    # determine common start/end points
    #plot(emb,bg='white',pch=21,col=ac(1,alpha=cell.color.alpha), xlab=xlab, ylab=ylab);
    lapply(1:length(spuci.pam$id.med),function(cn) {
      if(length(usps[[spuci.pam$id.med[cn]]])>0){
        mp <- usps[[spuci.pam$id.med[cn]]]; mp <- mp[!duplicated(mp)]
        bp <- data.frame(do.call(cbind,xspline(emb[mp,],shape=trajectory.spline.shape,draw=F)))
        lines(bp$x,bp$y,col=ac(1,alpha=0.6))
        bp <- bp[abs(diff(bp$x))+abs(diff(bp$y))>1e-5,]
        ai <- round(length(bp$x)*c(0.2,0.8,0.5))
        arrows(bp$x[ai],bp$y[ai],bp$x[ai+1],bp$y[ai+1],angle=30,length=0.1,col=ac(1,alpha=0.6))
      }
    })
    return(invisible(list(spuci.dist=spuci.dist,sps=sps,tp=tp,cpl=cpl)))
  } else if(!is.null(show.cell.trajectories)) {
    # show optimal path(s) for a particular cell
    celli <- match(show.cell.trajectories,rownames(emb));
    if(is.na(celli)) stop(paste('specified cell',show.cell.trajectories,'is not in the embedding'))
    
    cp <- Matrix::Diagonal(ncol(tp)); # cell position probabilities
    rownames(cp) <- colnames(cp) <- rownames(tp);
    
    #cpt <- as.array(cp);
    cpl <- list(); cpl[[1]] <- cp;
    #cp <- as.matrix(cp); tp <- as.matrix(tp)
    
    #ep <- as.array(emb)
    ttp <- t(tp);
    
    # run diffusion steps to figure out end positions
    cat("simulating diffusion ... ")
    for(i in 1:diffusion.steps) {
      cp <- cp %*% ttp;
      #cp[cp<1e-5] <- 0;
      cpl[[i+1]] <- cp;
    }
    
    # rank final points by probability
    cpo <- t(apply(-cp,1,order))
    
    # graph-based walkback approach
    
    # construct a walkback graph
    trp <- as(ttp,'TsparseMatrix')
    
    cat("constructing path graph ... ")
    x <- do.call(rbind,lapply(1:(diffusion.steps+1),function(i) {
      cbind(i=trp@i+1 + (i-1)*nrow(cp), # current time step
            j=trp@j+1 + (i)*nrow(cp))
    }))
    x <- x[x[,2]<=nrow(cp)*(diffusion.steps+1),]
    x <- Matrix::spMatrix(nrow(cp)*(diffusion.steps+1),nrow(cp)*(diffusion.steps+1),i=x[,1],j=x[,2],x=rep(1,nrow(x)))
    g <- igraph::graph.adjacency(x,mode='directed')
    rm(x); gc();
    
    # find topn trajectories for each cell
    cat("tracing shortest trajectories ... ")
    top.desti <- order(cp[celli,],decreasing=TRUE)[1:ntop.trajectories]
    
    
    # calculate cell-specific weights
    cw <- unlist(lapply(cpl,function(d) as.numeric(trp@x*(d[celli,trp@i+1]))))
    cw <- cw[1:igraph::ecount(g)] # trim extra edges
    # convert into penalty scores
    cw <- -log(cw)
    sp <- igraph::shortest_paths(g,from=celli,to=nrow(cp)*(diffusion.steps)+top.desti,weights=cw,mode='out')
    # remove time offset on the path nodes
    sp <- lapply(sp$vpath,function(x) { y <- (as.integer(x) %% nrow(cp)); y[y==0] <- nrow(cp); y});
    names(sp) <- rownames(cp)[top.desti]
    cat("done.\n")
    lapply(sp,function(mp) {
      if(!is.null(mp) && length(mp)>0) {
        mp <- mp[!duplicated(mp)]
        if(length(mp)>1)  {
          lines(emb[mp,1],emb[mp,2],col=8,lty=3)
          bp <- data.frame(do.call(cbind,xspline(emb[mp,],shape=trajectory.spline.shape,draw=F)))
          lines(bp$x,bp$y,col=ac(1,alpha=0.6))
          bp <- bp[abs(diff(bp$x))+abs(diff(bp$y))>1e-5,]
          ai <- round(length(bp$x)*c(0.2,0.8,0.5))
          arrows(bp$x[ai],bp$y[ai],bp$x[ai+1],bp$y[ai+1],angle=30,length=0.1,col=ac(1,alpha=0.6))
        }
      }
    })
    return(invisible(sp))
  } else if(show.all.trajectories) { # show diffusion paths
    cp <- Matrix::Diagonal(ncol(tp)); # cell position probabilities row-from col-to
    rownames(cp) <- colnames(cp) <- rownames(tp);
    ep <- as.array(emb)
    
    for(i in 1:diffusion.steps) {
      cp <- cp %*% t(tp);
      # expected position
      cpm <- t(apply(cp,1,function(x) { x[x<sort(x,decreasing=TRUE)[3]] <- 0; x/sum(x) }))
      epi <- as.matrix(cpm %*% emb);
      #epi <- as.matrix(cp %*% emb);
      ep <- abind::abind(ep,epi,along=3)
    }
    apply(ep,c(1),function(d) {
      lines(d[1,],d[2,],lwd=1,col=ac(1,alpha=0.05))
    })
  } else {
    # calculate arrows, draw
    lapply(1:nrow(emb),function(i) {
      # normalized directions to each point
      di <- t(t(emb)-emb[i,])
      di <- di/sqrt(Matrix::rowSums(di^2))*arrow.scale; di[i,] <- 0;
      di <- Matrix::colSums(di*tp[,i])
      suppressWarnings(arrows(emb[colnames(em)[i],1],emb[colnames(em)[i],2],emb[colnames(em)[i],1]+di[1],emb[colnames(em)[i],2]+di[2],length=0.05,lwd=arrow.lwd))
    })
  }
}



