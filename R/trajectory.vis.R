#' trajectory visualization 
#' @param emb embedding to be used for projection
#' @param tp transition probability 
#' @param current current matrix
#' @param emb_local_influence whether consider the emb influence
#' @param local_influence_prob the proportion of the influence
#' @param sigma_d the hyper-parameter of local influence
#' @param cell.colors name vector of cell colors
#' @param show.cell.trajectories show trajectories for a specified cell
#' @param show.trajectories show top median diffusion trajectories
#' @param show.all.trajectories show all diffusion paths (messy)
#' @param arrow.lwd arrow line width
#' @param xlab x axis label
#' @param ylab y axls label
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
#' @param point.size size of the point
#' @param low_color the hyper-parameter for show site low color
#' @param high_color the hyper-parameter for show site high color
#' @param arrow_size size of arrow
#' 
#' @importFrom graphics plot.new 
#' @importFrom grDevices dev.off
#' @export
#'
trajectory.vis <- function( emb,tp,current,emb_local_influence=FALSE,local_influence_prob=0.2,
                            sigma_d=NA,arrow.scale=1,
                            show.cell.trajectories=NULL,cell.colors=NULL,
                            show.trajectories=FALSE,arrow.lwd=1,
                            show.all.trajectories=FALSE,n.cores=1,
                            show.cell.diffusion.posterior=NULL, xlab="", ylab="",
                            show.grid.flow=FALSE, diffusion.steps=10,
                            cell.border.alpha=0.5,
                            ntop.trajectories=1, trajectory.spline.shape=1,
                            n.trajectory.clusters=10,point.size=3,
                            low_color="#FFCC15",high_color="#9933FF",
                            arrow_size=0.1,...){
  if (is.null(show.cell.diffusion.posterior)){
    celcol <- cell.colors[rownames(emb)] 
  } else{
    celcol <- "white"
  }
  
 
  
  
  em <- current; 
  ccells <- intersect(intersect(rownames(emb),colnames(em)),rownames(tp));
  em <- em[,ccells]; emb <- emb[ccells,];tp <- tp[ccells,ccells]
  
  point.data.frame <- as.data.frame(emb[ccells,1])
  point.data.frame$PC1 <- emb[ccells,1]
  point.data.frame$PC2 <- emb[ccells,2]
  point.data.frame$color <- celcol[ccells]
  point.data.frame <- point.data.frame[,c("PC1","PC2","color")]
  
  
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
   
    data.frame.poster <- as.data.frame(emb[,1])
    data.frame.poster$PC1 <- emb[,1]
    data.frame.poster$PC2 <- emb[,2]
    data.frame.poster$color <- as.numeric(cp[show.cell.diffusion.posterior,])
    
    plot.ggplot <- ggplot() + 
      geom_point(data=data.frame.poster,
                 aes(x=.data$PC1, y=.data$PC2, fill=.data$color),
                 alpha=cell.border.alpha,size=point.size,shape=21)+
      labs(x="PC1", y="PC2", 
           title="PC1 vs PC2") + 
      theme(plot.title=element_text(size=12,hjust=0.5)) +
      scale_fill_gradient(low = low_color, high = high_color)+
      geom_point(data=point.data.frame[show.cell.diffusion.posterior,],
                 aes(x=.data$PC1, y=.data$PC2),color="black",shape=3,size=3)
    
    
    return(plot.ggplot)
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
    if(n.trajectory.clusters>=dim(spuci.dist)[1]) {
      cat("n.trajectory.clusters is too big ... set to be ")
      n.trajectory.clusters = dim(spuci.dist)[1]-1
      cat(n.trajectory.clusters)
    }
    spuci.pam <- cluster::pam(spuci.dist,n.trajectory.clusters)
    cat("done.\n")
    
    # bezier
    # determine common start/end points
    #plot(emb,bg='white',pch=21,col=ac(1,alpha=cell.color.alpha), xlab=xlab, ylab=ylab);
    
    plot.ggplot <- ggplot() +
      geom_point(data=point.data.frame, aes(x=.data$PC1, y=.data$PC2,
         fill=.data$color), alpha=cell.border.alpha,size=point.size,shape=21) +
      labs(x="PC1", y="PC2",
           title="PC1 vs PC2") +
      theme(plot.title=element_text(size=12,hjust=0.5))
    
    
    for (cn in 1:length(spuci.pam$id.med)) {
      if (length(usps[[spuci.pam$id.med[cn]]]) > 0) {
        mp <- usps[[spuci.pam$id.med[cn]]]
        mp <- mp[!duplicated(mp)]
        emb_data <- data.frame(emb[mp, ])
        colnames(emb_data) <- c("x", "y")
        emb_data <- data.frame(emb[mp, ])
        colnames(emb_data) <- c("x", "y")
        
        
        # Initialize a new plot (necessary for xspline to work)
        plot.new()
        
        # Use xspline to generate spline points
        bp <- as.data.frame(xspline(emb_data$x, emb_data$y, 1, draw = FALSE))
        
        dev.off()
        # Close the invisible plotting device
        # bp <- as.data.frame(bp_$x)
        # 
        # bp$y <- bp_$y
        # colnames(bp) <- c("x", "y")
        # bp <- bp[abs(diff(emb_data$x)) + abs(diff(emb_data$y)) > 1e-5, ]
        plot.ggplot <- plot.ggplot +
          geom_path(data = bp, aes(x = .data$x, y = .data$y),
                     color = "black",  size = arrow.lwd)#, span = 0.3, method = "loess", alpha = cell.border.alpha
        
        bp <- bp[abs(diff(emb_data$x)) + abs(diff(emb_data$y)) > 1e-5, ]
        ai <- round(length(bp$x) * c(0.2, 0.8, 0.5))
        
        bp_arrow <- data.frame(x0 = bp$x[ai], y0 = bp$y[ai], x1 = bp$x[ai + 1], y1 = bp$y[ai + 1])

        plot.ggplot <- plot.ggplot +
          geom_segment(data = bp_arrow, aes(x = .data$x0, y = .data$y0, xend = .data$x1, yend = .data$y1),
                       arrow = arrow(length = unit(arrow_size, "inches")), size = arrow.lwd) +
          theme(legend.position = "none")
      }
    }
    return(list(plot.ggplot,invisible(list(spuci.dist=spuci.dist,sps=sps,tp=tp,cpl=cpl))))
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
    
    plot.ggplot <- ggplot() +
      geom_point(data=point.data.frame, aes(x=.data$PC1, y=.data$PC2,
                                            fill=.data$color), alpha=cell.border.alpha,size=point.size,shape=21) +
      labs(x="PC1", y="PC2",
           title="PC1 vs PC2") +
      theme(plot.title=element_text(size=12,hjust=0.5))
    
    
    for (i in seq_along(sp)) {
      mp <- sp[[i]]
      if (!is.null(mp) && length(mp) > 0) {
        mp <- mp[!duplicated(mp)]
        if(length(mp)>1)  {
          emb_data <- data.frame(emb[mp, ])
          colnames(emb_data) <- c("x", "y")
          plot.new()
          
          # Use xspline to generate spline points
          bp <- as.data.frame(xspline(emb_data$x, emb_data$y, 1, draw = FALSE))
          
          dev.off()
          
          colnames(bp) <- c("x", "y")
          plot.ggplot <- plot.ggplot +
            geom_path(data = bp, aes(x = .data$x, y = .data$y),
                        color = "black",  size = arrow.lwd, alpha = cell.border.alpha)
          plot.ggplot <- plot.ggplot +
            geom_path(data = emb_data, aes(x = .data$x, y = .data$y),
                      color = "black",  size = arrow.lwd, alpha = cell.border.alpha)
          bp <- bp[abs(diff(bp$x))+abs(diff(bp$y))>1e-5,]
          ai <- round(length(bp$x)*c(0.2,0.8,0.5))
          bp_arrow <- data.frame(x0 = bp$x[ai], y0 = bp$y[ai], x1 = bp$x[ai + 1], y1 = bp$y[ai + 1])
          
          plot.ggplot <- plot.ggplot +
            geom_segment(data = bp_arrow, aes(x = .data$x0, y = .data$y0, xend = .data$x1, yend = .data$y1),
                         arrow = arrow(length = unit(arrow_size, "inches")), size = arrow.lwd) +
            theme(legend.position = "none")
        }
      }
    }
    
    
    return(list(plot.ggplot,invisible(sp)))
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
    plot.ggplot <- ggplot() +
      geom_point(data=point.data.frame, aes(x=.data$PC1, y=.data$PC2,
                                            fill=.data$color), alpha=cell.border.alpha,size=point.size,shape=21) +
      labs(x="PC1", y="PC2",
           title="PC1 vs PC2") +
      theme(plot.title=element_text(size=12,hjust=0.5))
    
    
    # for (i in 1:length(ep)){
    #   plot.ggplot <- plot.ggplot +
    #     geom_path(data = ep, aes(x = .data$x, y = .data$y),
    #               color = "black",  size = arrow.lwd, alpha = cell.border.alpha)
    # }
    #  apply(ep,c(1),function(d) {
    #    lines(d[1,],d[2,],lwd=1,col=ac(1,alpha=0.05))
    #  })
    for (i in 1:dim(ep)[1]) {
      plot.ggplot <- plot.ggplot  +
        geom_path(data = data.frame(x = ep[i,1 , ], y = ep[i,2 , ]), aes(x = .data$x, y = .data$y),
                  alpha = cell.border.alpha, size = 1, color = "black")
    }
    return(plot.ggplot)
  } else {
    stop("please change the parameters to output the trajectory")
  }
}



