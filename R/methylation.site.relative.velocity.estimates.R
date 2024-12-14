##' Estimate RNA methylation velocity using site-relative slopes
##'
##' @param test.list  a list of the test cells, with four data.frame of the methylation.spliced,  unmethylation.spliced,methylation.unspliced, unmethylation.unspliced read counts
##' @param control.list   a list of the control cells, with four data.frame of the methylation.spliced,  unmethylation.spliced,methylation.unspliced, unmethylation.unspliced read counts
##' @param deltaT - amount of time to project the cell forward
##' @param steady.state.cells - optional set of steady-state cells on which the gamma should be estimated (defaults to all cells)
##' @param kCells - number of k nearest neighbors (NN) to use in slope calculation smoothing
##' @param cellKNN - optional pre-calculated cell KNN matrix
##' @param kSites - number of sites (k) to use in site kNN pooling
##' @param siteKNN - optional pre-calculated site KNN matrix
##' @param old.fit - optional old result (in this case the slopes and offsets won't be recalculated, and the same kNN graphs will be used)
##' @param mult - library scaling factor (1e6 in case of FPM)
##' @param min.nmat.emat.correlation - minimum required Spearman rank correlation between n and e counts of a site
##' @param p_value calculate the p-value or not
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
##' @param epsilon_Beta epsilon of Beta
##' @param epsilon_M epsilon of M
##' @param epsilon_OR epsilon of OR
##' @param epsilon_RR epsilon of RR
##' @param epsilon_TCR epsilon of TCR
##' @param epsilon_control_M epsilon of control_M
##' @param epsilon_control_Beta epsilon of control_Beta
##' @param narm size factor rm na or not
#'@importFrom stats anova
#' @return list(meth_result,unmeth_result,methylation_level_list,control_information)
#' @export
#'
#'
methylation.site.relative.velocity.estimates <- function (
    test.list, control.list, deltaT = 1, steady.state.cells = colnames(test.list[[1]]),
    kCells = 5, cellKNN = NULL, kSites = 1, siteKNN = NULL, old.fit = NULL,
    mult = 1000, min.nmat.emat.correlation = 0.05,p_value=TRUE,
    min.nmat.emat.slope = 0.05, zero.offset = FALSE, deltaT2 = 1,delta_model="model 2",
    fit.quantile = NULL, diagonal.quantiles = FALSE, show.site = NULL,
    cell.dist = NULL, emat.size = NULL, nmat.size = NULL,
    cell.emb = NULL, cell.colors = NULL, expression.gradient = NULL,
    residual.gradient = NULL, n.cores = 1, verbose = TRUE,nrows=2,
    low_color="#FFCC15",high_color="#9933FF",point.size=3,epsilon_Beta=1e-7,
    epsilon_M=1e-7,epsilon_OR=1e-7,epsilon_RR=1e-7,epsilon_TCR=1e-7,
    epsilon_control_M=1,epsilon_control_Beta=1,narm=FALSE)
{
  
  # check test.list 
  if (length(test.list)!=4)
    stop("the length of test list need to be 4")
  
  if (dim(test.list[[1]])[1]!=dim(test.list[[2]])[1])
    stop("the dimension of test list need to be same")
  if (dim(test.list[[1]])[1]!=dim(test.list[[3]])[1])
    stop("the dimension of test list need to be same")
  if (dim(test.list[[1]])[1]!=dim(test.list[[4]])[1])
    stop("the dimension of test list need to be same")
  if (dim(test.list[[1]])[2]!=dim(test.list[[2]])[2])
    stop("the dimension of test list need to be same")
  if (dim(test.list[[1]])[2]!=dim(test.list[[3]])[2])
    stop("the dimension of test list need to be same")
  if (dim(test.list[[1]])[2]!=dim(test.list[[4]])[2])
    stop("the dimension of test list need to be same")
  
  spliced.meth.test <- test.list[[1]]
  spliced.unmeth.test <- test.list[[2]]
  unspliced.meth.test <- test.list[[3]]
  unspliced.unmeth.test <- test.list[[4]]
  
  if (!is.null(control.list)){
    if (length(control.list)!=4)
      stop("the length of control list need to be 4")
    
    if (dim(control.list[[1]])[1]!=dim(control.list[[2]])[1])
      stop("the dimension of control list need to be same")
    if (dim(control.list[[1]])[1]!=dim(control.list[[3]])[1])
      stop("the dimension of control list need to be same")
    if (dim(control.list[[1]])[1]!=dim(control.list[[4]])[1])
      stop("the dimension of control list need to be same")
    if (dim(control.list[[1]])[2]!=dim(control.list[[2]])[2])
      stop("the dimension of control list need to be same")
    if (dim(control.list[[1]])[2]!=dim(control.list[[3]])[2])
      stop("the dimension of control list need to be same")
    if (dim(control.list[[1]])[2]!=dim(control.list[[4]])[2])
      stop("the dimension of control list need to be same")
    spliced.meth.control <- control.list[[1]]
    spliced.unmeth.control <- control.list[[2]]
    unspliced.meth.control <- control.list[[3]]
    unspliced.unmeth.control <- control.list[[4]]
  }
  
  
  
  meth_result <- site.relative.velocity.estimates(emat = spliced.meth.test, nmat = unspliced.meth.test, deltaT = deltaT, steady.state.cells = colnames(spliced.meth.test),
                                                  kCells = kCells, cellKNN = cellKNN, kSites = kSites, siteKNN = siteKNN, old.fit = old.fit,
                                                  mult = mult, min.nmat.emat.correlation = min.nmat.emat.correlation,
                                                  min.nmat.emat.slope = min.nmat.emat.slope, zero.offset = zero.offset, deltaT2 = deltaT2,delta_model=delta_model,
                                                  fit.quantile = fit.quantile, diagonal.quantiles = diagonal.quantiles, show.site = show.site,
                                                  cell.dist = cell.dist, emat.size = emat.size, nmat.size = nmat.size,
                                                  cell.emb = cell.emb, cell.colors = cell.colors, expression.gradient = expression.gradient,
                                                  residual.gradient = residual.gradient, n.cores = n.cores, verbose = verbose,nrows=nrows,
                                                  low_color=low_color,high_color=high_color,point.size=point.size)
  
  unmeth_result <- site.relative.velocity.estimates(emat = spliced.unmeth.test, nmat = unspliced.unmeth.test, deltaT = deltaT, steady.state.cells = colnames(spliced.unmeth.test),
                                                  kCells = kCells, cellKNN = cellKNN, kSites = kSites, siteKNN = siteKNN, old.fit = old.fit,
                                                  mult = mult, min.nmat.emat.correlation = min.nmat.emat.correlation,
                                                  min.nmat.emat.slope = min.nmat.emat.slope, zero.offset = zero.offset, deltaT2 = deltaT2,delta_model=delta_model,
                                                  fit.quantile = fit.quantile, diagonal.quantiles = diagonal.quantiles, show.site = show.site,
                                                  cell.dist = cell.dist, emat.size = emat.size, nmat.size = nmat.size,
                                                  cell.emb = cell.emb, cell.colors = cell.colors, expression.gradient = expression.gradient,
                                                  residual.gradient = residual.gradient, n.cores = n.cores, verbose = verbose,nrows=nrows,
                                                  low_color=low_color,high_color=high_color,point.size=point.size)
  
  
  index_meth <- rownames(meth_result[["deltaE"]])
  index_unmeth <- rownames(unmeth_result[["deltaE"]])
  
  if(length(intersect(index_meth,index_unmeth))==0){
    stop("the results of meth and unmeth without any intersection, please change the parameter")
  }
  index_intersect <- intersect(index_meth,index_unmeth)
  
  meth_result[["gamma_intersect"]]<- meth_result[["gamma"]][index_intersect]
  meth_result[["projected_intersect"]]<- meth_result[["projected"]][index_intersect,]
  meth_result[["current_intersect"]]<- meth_result[["current"]][index_intersect,]
  meth_result[["deltaE_intersect"]]<- meth_result[["deltaE"]][index_intersect,]
  
  unmeth_result[["gamma_intersect"]]<- unmeth_result[["gamma"]][index_intersect]
  unmeth_result[["projected_intersect"]]<- unmeth_result[["projected"]][index_intersect,]
  unmeth_result[["current_intersect"]]<- unmeth_result[["current"]][index_intersect,]
  unmeth_result[["deltaE_intersect"]]<- unmeth_result[["deltaE"]][index_intersect,]
  
  if (p_value==TRUE){
    cat("calculating P-value ... ")
    
    # Initialize vector to store p-values
    p_value <- c()
    
    # Loop over rows of the data
    for (i in 1:dim(meth_y_value)[1]) {
      # Extract the cells for meth and unmeth separately
      meth_cells <- meth_result[["cell.use.list"]][index_intersect][[i]]
      unmeth_cells <- unmeth_result[["cell.use.list"]][index_intersect][[i]]
      
      # Create the data frame with separate observations for meth and unmeth
      data.frame.pvalue <- as.data.frame(rbind(
        data.frame(y = meth_y_value[i, meth_cells], 
                   x = meth_x_value[i, meth_cells], 
                   var = "meth"),
        data.frame(y = unmeth_y_value[i, unmeth_cells], 
                   x = unmeth_x_value[i, unmeth_cells], 
                   var = "unmeth")
      ))
      
      # Ensure columns are correctly typed
      data.frame.pvalue$var <- as.character(data.frame.pvalue$var)
      data.frame.pvalue$x <- as.numeric(data.frame.pvalue$x)
      data.frame.pvalue$y <- as.numeric(data.frame.pvalue$y)
      
      # Perform linear modeling
      lm_res <- lm(y ~ 0 + x + var, data = data.frame.pvalue)
      
      # Extract p-value from ANOVA
      anova_result <- anova(lm_res)
      p_value[i] <- anova_result["var", "Pr(>F)"]
    }
    
    # Name the p-values by the intersected indices
    names(p_value) <- index_intersect
    
    cat("done\n")
  }
  
  
  
  # calculate the methylation level
  cat("calculating methylation level ... ")
  
  Beta_value_current <- meth_result[["current"]][index_intersect,]/
    (meth_result[["current"]][index_intersect,]+unmeth_result[["current"]][index_intersect,]+epsilon_Beta)
  M_value_current <- log2((meth_result[["current"]][index_intersect,]+epsilon_M)/
    (unmeth_result[["current"]][index_intersect,]+epsilon_M))
  
  Beta_value_projected <- meth_result[["projected"]][index_intersect,]/
    (meth_result[["current"]][index_intersect,]+unmeth_result[["projected"]][index_intersect,]+epsilon_Beta)
  M_value_projected <- log2((meth_result[["projected"]][index_intersect,]+epsilon_M)/
                            (unmeth_result[["projected"]][index_intersect,]+epsilon_M))
  Beta_value_delta <- Beta_value_projected - Beta_value_current
  M_value_delta <- M_value_projected - M_value_current
  
  
  if (is.null(control.list)){
    OR_current <- RR_current <- TCR_current <- control_information <- NULL
    OR_projected <- RR_projected <- TCR_projected <- NULL
    OR_delta <- RR_delta <- TCR_delta <- NULL
  } else{
    control_information <- list()
    
    size <- sizeFactor(spliced.meth.control + spliced.unmeth.control,narm = narm)
    if (anyNA(size)) {
      size <- sizeFactor(spliced.meth.control + spliced.unmeth.control,narm = TRUE)
    }
    control_information_M <- (rowSums( t(t(spliced.meth.control)/size) )+epsilon_control_M) /
      (rowSums( t(t(spliced.unmeth.control)/size) )+epsilon_control_M)
    
    control_information_Beta <- (rowSums( t(t(spliced.meth.control)/size) )+epsilon_control_Beta) /
      (rowSums( t(t(spliced.unmeth.control)/size) )+rowSums( t(t(spliced.meth.control)/size) )+epsilon_control_Beta)
    
    control_information[["control_information_M"]] <- control_information_M
    control_information[["control_information_Beta"]] <- control_information_Beta
    
    OR_current <- log2(((meth_result[["current"]][index_intersect,]+epsilon_OR)/
                         (unmeth_result[["current"]][index_intersect,]+epsilon_OR))  /
                         control_information_M[index_intersect])
    RR_current <- log2(((meth_result[["current"]][index_intersect,]+epsilon_RR)/
                          (unmeth_result[["current"]][index_intersect,]+
                             meth_result[["current"]][index_intersect,]+epsilon_RR))  /
                         control_information_Beta[index_intersect])
    TCR_current <- ((meth_result[["current"]][index_intersect,]+epsilon_TCR)/
                          (unmeth_result[["current"]][index_intersect,]
                           +meth_result[["current"]][index_intersect,]+epsilon_TCR))- 
                         control_information_Beta[index_intersect]
    TCR_current[TCR_current<=0] <- 0
    
    OR_projected <- log2(((meth_result[["projected"]][index_intersect,]+epsilon_OR)/
                          (unmeth_result[["projected"]][index_intersect,]+epsilon_OR))  /
                         control_information_M[index_intersect])
    RR_projected <- log2(((meth_result[["projected"]][index_intersect,]+epsilon_RR)/
                          (unmeth_result[["projected"]][index_intersect,]+
                             meth_result[["projected"]][index_intersect,]+epsilon_RR))  /
                         control_information_Beta[index_intersect])
    TCR_projected <- ((meth_result[["projected"]][index_intersect,]+epsilon_TCR)/
                          (unmeth_result[["projected"]][index_intersect,]
                           +meth_result[["projected"]][index_intersect,]+epsilon_TCR))- 
                           control_information_Beta[index_intersect]
    
    TCR_projected[TCR_projected<=0] <- 0
    OR_delta <- OR_projected - OR_current
    RR_delta <- RR_projected - RR_current
    TCR_delta <- TCR_projected - TCR_current
    
    
  }
  
  
  methylation_level_list <- list()
  methylation_level_list[["current"]] <- list()
  methylation_level_list[["current"]][["M"]] <- M_value_current
  methylation_level_list[["current"]][["Beta"]] <- Beta_value_current
  methylation_level_list[["current"]][["OR"]] <- OR_current
  methylation_level_list[["current"]][["RR"]] <- RR_current
  methylation_level_list[["current"]][["TCR"]] <- TCR_current
  
  methylation_level_list[["projected"]] <- list()
  methylation_level_list[["projected"]][["M"]] <- M_value_projected
  methylation_level_list[["projected"]][["Beta"]] <- Beta_value_projected
  methylation_level_list[["projected"]][["OR"]] <- OR_projected
  methylation_level_list[["projected"]][["RR"]] <- RR_projected
  methylation_level_list[["projected"]][["TCR"]] <- TCR_projected
  
  methylation_level_list[["delta"]] <- list()
  methylation_level_list[["delta"]][["M"]] <- M_value_delta
  methylation_level_list[["delta"]][["Beta"]] <- Beta_value_delta
  methylation_level_list[["delta"]][["OR"]] <- OR_delta
  methylation_level_list[["delta"]][["RR"]] <- RR_delta
  methylation_level_list[["delta"]][["TCR"]] <- TCR_delta
  
  methylation_level_list[["p_value"]] <- p_value
  
  cat("done\n")
  
  
  return(list(meth_result,unmeth_result,methylation_level_list,control_information))
}



