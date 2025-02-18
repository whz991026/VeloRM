#' Title
#'
#' @param test.list a list of the test cells, with four data.frame of the methylation.spliced,  unmethylation.spliced,methylation.unspliced, unmethylation.unspliced read counts
#' @param control.list   a list of the control cells, with four data.frame of the methylation.spliced,  unmethylation.spliced,methylation.unspliced, unmethylation.unspliced read counts
#' @param control_size size factor for the control cell
#' @param methylation.site.relative.velocity.result results from methylation.site.relative.velocity.estimates
#' @param exp.site.velocity.result results from site.relative.velocity.estimates using meth reads + unmeth reads
#' @param delta_exp c("delta_from_exp_velocity", "delta_from_meta_methylation_site_velocity")
#' @param methylation_level c("M","Beta","RR","OR","TCR","meta_M_sep","meta_Beta_sep","meta_RR_sep","meta_OR_sep","meta_TCR_sep")
#' @param methylation_level_scale scale of the methylation level c("linear","log2")
#' @param cor_method c("pearson", "kendall", "spearman")
#' @param epsilon_Beta epsilon of Beta
#' @param epsilon_M epsilon of M
#' @param epsilon_OR epsilon of OR
#' @param epsilon_RR epsilon of RR
#' @param epsilon_TCR epsilon of TCR
#' @param epsilon_control_M epsilon of control_M
#' @param epsilon_control_Beta epsilon of control_Beta
#' @param na.action.lm 	a function which indicates what should happen when the data contain NAs. The default is set by the na.action setting of options, and is na.fail if that is unset. The ‘factory-fresh’ default is na.omit. Another possible value is NULL, no action. Value na.exclude can be useful.
#' @param narm size factor rm na or not
#' 
#' @importFrom stats na.exclude
#' 
#' @return list()
#' @export
#'
transcriptional.impact.analysis <- function(test.list,control.list=NULL,control_size=NULL,
                                            methylation.site.relative.velocity.result,
                                            exp.site.velocity.result,delta_exp=NULL,
                                            methylation_level,methylation_level_scale,
                                            cor_method,epsilon_Beta=1e-7,
                                            epsilon_M=1e-7,epsilon_OR=1e-7,epsilon_RR=1e-7,
                                            epsilon_TCR=1e-7,
                                            epsilon_control_M=1,epsilon_control_Beta=1,
                                            na.action.lm= "na.exclude",narm=FALSE){
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
  
  if (!is.null(control.list)){
    
    ##
    control_information_spliced <- list()
    if(length(control_size)==0){
      control_size <- sizeFactor(spliced.meth.control + spliced.unmeth.control,narm = narm)
      if (anyNA(control_size)) {
        control_size <- sizeFactor(spliced.meth.control + spliced.unmeth.control,narm = TRUE)
      }
      if (anyNA(control_size)) {
        index <- which(is.na(control_size))
        print(paste0("there are ",length(index), "number of cell missing the control size"))
        control_information_spliced_M <- (rowSums( t(t(spliced.meth.control[,-index])/control_size[-index]) )+epsilon_control_M) /
          (rowSums( t(t(spliced.unmeth.control[,-index])/control_size[-index]) )+epsilon_control_M)
        
        control_information_spliced_Beta <- (rowSums( t(t(spliced.meth.control[,-index])/control_size[-index]) )+epsilon_control_Beta) /
          (rowSums( t(t(spliced.unmeth.control[,-index])/control_size[-index]) )+rowSums( t(t(spliced.meth.control[,-index])/control_size[-index]) )+epsilon_control_Beta)
        
        control_information_spliced[["control_information_M"]] <- control_information_spliced_M
        control_information_spliced[["control_information_Beta"]] <- control_information_spliced_Beta
      } else{
        control_information_spliced_M <- (rowSums( t(t(spliced.meth.control)/control_size) )+epsilon_control_M) /
          (rowSums( t(t(spliced.unmeth.control)/control_size) )+epsilon_control_M)
        
        control_information_spliced_Beta <- (rowSums( t(t(spliced.meth.control)/control_size) )+epsilon_control_Beta) /
          (rowSums( t(t(spliced.unmeth.control)/control_size) )+rowSums( t(t(spliced.meth.control)/control_size) )+epsilon_control_Beta)
        
        control_information_spliced[["control_information_M"]] <- control_information_spliced_M
        control_information_spliced[["control_information_Beta"]] <- control_information_spliced_Beta
      }
    }else{
      control_information_spliced_M <- (rowSums( t(t(spliced.meth.control)/control_size) )+epsilon_control_M) /
        (rowSums( t(t(spliced.unmeth.control)/control_size) )+epsilon_control_M)
      
      control_information_spliced_Beta <- (rowSums( t(t(spliced.meth.control)/control_size) )+epsilon_control_Beta) /
        (rowSums( t(t(spliced.unmeth.control)/control_size) )+rowSums( t(t(spliced.meth.control)/control_size) )+epsilon_control_Beta)
      
      control_information_spliced[["control_information_M"]] <- control_information_spliced_M
      control_information_spliced[["control_information_Beta"]] <- control_information_spliced_Beta
    }
    
    
    ##
    control_information_unspliced <- list()
    if(length(control_size)==0){
      control_size <- sizeFactor(unspliced.meth.control + unspliced.unmeth.control,narm = narm)
      if (anyNA(control_size)) {
        control_size <- sizeFactor(unspliced.meth.control + unspliced.unmeth.control,narm = TRUE)
      }
      if (anyNA(control_size)) {
        index <- which(is.na(control_size))
        print(paste0("there are ",length(index), "number of cell missing the control size"))
        control_information_unspliced_M <- (rowSums( t(t(unspliced.meth.control[,-index])/control_size[-index]) )+epsilon_control_M) /
          (rowSums( t(t(unspliced.unmeth.control[,-index])/control_size[-index]) )+epsilon_control_M)
        
        control_information_unspliced_Beta <- (rowSums( t(t(unspliced.meth.control[,-index])/control_size[-index]) )+epsilon_control_Beta) /
          (rowSums( t(t(unspliced.unmeth.control[,-index])/control_size[-index]) )+rowSums( t(t(unspliced.meth.control[,-index])/control_size[-index]) )+epsilon_control_Beta)
        
        control_information_unspliced[["control_information_M"]] <- control_information_unspliced_M
        control_information_unspliced[["control_information_Beta"]] <- control_information_unspliced_Beta
      } else{
        control_information_unspliced_M <- (rowSums( t(t(unspliced.meth.control)/control_size) )+epsilon_control_M) /
          (rowSums( t(t(unspliced.unmeth.control)/control_size) )+epsilon_control_M)
        
        control_information_unspliced_Beta <- (rowSums( t(t(unspliced.meth.control)/control_size) )+epsilon_control_Beta) /
          (rowSums( t(t(unspliced.unmeth.control)/control_size) )+rowSums( t(t(unspliced.meth.control)/control_size) )+epsilon_control_Beta)
        
        control_information_unspliced[["control_information_M"]] <- control_information_unspliced_M
        control_information_unspliced[["control_information_Beta"]] <- control_information_unspliced_Beta
      }
    }else{
      control_information_unspliced_M <- (rowSums( t(t(unspliced.meth.control)/control_size) )+epsilon_control_M) /
        (rowSums( t(t(unspliced.unmeth.control)/control_size) )+epsilon_control_M)
      
      control_information_unspliced_Beta <- (rowSums( t(t(unspliced.meth.control)/control_size) )+epsilon_control_Beta) /
        (rowSums( t(t(unspliced.unmeth.control)/control_size) )+rowSums( t(t(unspliced.meth.control)/control_size) )+epsilon_control_Beta)
      
      control_information_unspliced[["control_information_M"]] <- control_information_unspliced_M
      control_information_unspliced[["control_information_Beta"]] <- control_information_unspliced_Beta
    }
    
    
    ##
    control_information_all <- list()
    if(length(control_size)==0){
      control_size <- sizeFactor(unspliced.meth.control + unspliced.unmeth.control+spliced.meth.control + spliced.unmeth.control,narm = narm)
      if (anyNA(control_size)) {
        control_size <- sizeFactor(unspliced.meth.control + unspliced.unmeth.control+spliced.meth.control + spliced.unmeth.control,narm = TRUE)
      }
      if (anyNA(control_size)) {
        index <- which(is.na(control_size))
        print(paste0("there are ",length(index), "number of cell missing the control size"))
        control_information_all_M <- (rowSums( t(t(unspliced.meth.control[,-index]+spliced.meth.control[,-index])/control_size[-index]) )+epsilon_control_M) /
          (rowSums( t(t(unspliced.unmeth.control[,-index]+spliced.unmeth.control[,-index])/control_size[-index]) )+epsilon_control_M)
        
        control_information_all_Beta <- (rowSums( t(t(unspliced.meth.control[,-index]+spliced.meth.control[,-index])/control_size[-index]) )+epsilon_control_Beta) /
          (rowSums( t(t(unspliced.unmeth.control[,-index]+spliced.unmeth.control[,-index])/control_size[-index]) )+rowSums( t(t(unspliced.meth.control[,-index]+spliced.meth.control[,-index])/control_size[-index]) )+epsilon_control_Beta)
        
        control_information_all[["control_information_M"]] <- control_information_all_M
        control_information_all[["control_information_Beta"]] <- control_information_all_Beta
      } else{
        control_information_all_M <- (rowSums( t(t(unspliced.meth.control+spliced.meth.control)/control_size) )+epsilon_control_M) /
          (rowSums( t(t(unspliced.unmeth.control+spliced.unmeth.control)/control_size) )+epsilon_control_M)
        
        control_information_all_Beta <- (rowSums( t(t(unspliced.meth.control+spliced.meth.control)/control_size) )+epsilon_control_Beta) /
          (rowSums( t(t(unspliced.unmeth.control+spliced.meth.control)/control_size) )+rowSums( t(t(unspliced.meth.control+spliced.meth.control)/control_size) )+epsilon_control_Beta)
        
        control_information_all[["control_information_M"]] <- control_information_all_M
        control_information_all[["control_information_Beta"]] <- control_information_all_Beta
      }
    }else{
      control_information_all_M <- (rowSums( t(t(unspliced.meth.control+spliced.meth.control)/control_size) )+epsilon_control_M) /
        (rowSums( t(t(unspliced.unmeth.control+spliced.unmeth.control)/control_size) )+epsilon_control_M)
      
      control_information_all_Beta <- (rowSums( t(t(unspliced.meth.control+spliced.meth.control)/control_size) )+epsilon_control_Beta) /
        (rowSums( t(t(unspliced.unmeth.control+spliced.unmeth.control)/control_size) )+rowSums( t(t(unspliced.meth.control+spliced.meth.control)/control_size) )+epsilon_control_Beta)
      
      control_information_all[["control_information_M"]] <- control_information_all_M
      control_information_all[["control_information_Beta"]] <- control_information_all_Beta
    }
    
  }
  
  
  
  
  
  
  
  # check meth and unmeth velocity result
  if (is.null(methylation.site.relative.velocity.result)){
    stop("methylation.site.relative.velocity.result cannot be NULL")
  }else{
    
    # check control_information
    if (is.null(control.list)){
      if(methylation_level%in%c("RR","OR","TCR","meta_RR_sep","meta_OR_sep","meta_TCR_sep")){
        stop("methylation level can not be RR, OR, TCR without control information")
      }
    } 
    
    if(is.null(exp.site.velocity.result)){
      intersect_names <- rownames(methylation.site.relative.velocity.result[[3]][["current"]][["M"]])
      if (is.null(control.list)){
        control_information_intersect <- NULL
      } else{
        control_information_intersect_spliced <- list()
        control_information_intersect_spliced [["control_information_M"]] <- control_information_spliced[["control_information_M"]][intersect_names] 
        control_information_intersect_spliced [["control_information_Beta"]] <- control_information_spliced[["control_information_Beta"]][intersect_names] 
        control_information_intersect_unspliced <- list()
        control_information_intersect_unspliced [["control_information_M"]] <- control_information_unspliced[["control_information_M"]][intersect_names] 
        control_information_intersect_unspliced [["control_information_Beta"]] <- control_information_unspliced[["control_information_Beta"]][intersect_names]
        control_information_intersect_all <- list()
        control_information_intersect_all [["control_information_M"]] <- control_information_all[["control_information_M"]][intersect_names] 
        control_information_intersect_all [["control_information_Beta"]] <- control_information_all[["control_information_Beta"]][intersect_names]
      }
      if(delta_exp=="delta_from_exp_velocity"){
        stop("delta_exp can not be delta_from_exp_velocity without exp.site.velocity.result")
      }else if (delta_exp=="delta_from_meta_methylation_site_velocity"){
        delta_exp_matrix_1 <- methylation.site.relative.velocity.result[[1]][["deltaE"]]
        delta_exp_matrix_2 <- methylation.site.relative.velocity.result[[2]][["deltaE"]]
        delta_exp_matrix_1 <- as.matrix(delta_exp_matrix_1)
        delta_exp_matrix_2 <- as.matrix(delta_exp_matrix_2)
        delta_exp_matrix <- delta_exp_matrix_1[intersect_names,] + delta_exp_matrix_2[intersect_names,]
      } else {
        stop("delta_exp should be delta_from_exp_velocity or delta_from_meta_methylation_site_velocity")
      }
    }else{
      if(delta_exp=="delta_from_exp_velocity"){
        delta_exp_matrix <- exp.site.velocity.result[["deltaE"]]
        delta_exp_matrix <- as.matrix(delta_exp_matrix)
        intersect_names <- intersect(rownames(methylation.site.relative.velocity.result[[3]][["current"]][["M"]]),
                                     rownames(exp.site.velocity.result[["current"]]))
        delta_exp_matrix <- delta_exp_matrix[intersect_names,]
      }else if (delta_exp=="delta_from_meta_methylation_site_velocity"){
        delta_exp_matrix_1 <- methylation.site.relative.velocity.result[[1]][["deltaE"]]
        delta_exp_matrix_2 <- methylation.site.relative.velocity.result[[2]][["deltaE"]]
        delta_exp_matrix_1 <- as.matrix(delta_exp_matrix_1)
        delta_exp_matrix_2 <- as.matrix(delta_exp_matrix_2)
        intersect_names <- rownames(methylation.site.relative.velocity.result[[3]][["current"]][["M"]])
        delta_exp_matrix <- delta_exp_matrix_1[intersect_names,] + delta_exp_matrix_2[intersect_names,]
        
      } else {
        stop("delta_exp should be delta_from_exp_velocity or delta_from_meta_methylation_site_velocity")
      }
      
      if (is.null(control.list)){
        control_information_intersect <- NULL
      } else{
        control_information_intersect_spliced <- list()
        control_information_intersect_spliced [["control_information_M"]] <- control_information_spliced[["control_information_M"]][intersect_names] 
        control_information_intersect_spliced [["control_information_Beta"]] <- control_information_spliced[["control_information_Beta"]][intersect_names] 
        control_information_intersect_unspliced <- list()
        control_information_intersect_unspliced [["control_information_M"]] <- control_information_unspliced[["control_information_M"]][intersect_names] 
        control_information_intersect_unspliced [["control_information_Beta"]] <- control_information_unspliced[["control_information_Beta"]][intersect_names]
        control_information_intersect_all <- list()
        control_information_intersect_all [["control_information_M"]] <- control_information_all[["control_information_M"]][intersect_names] 
        control_information_intersect_all [["control_information_Beta"]] <- control_information_all[["control_information_Beta"]][intersect_names]
        }
    }
  }
  
  
  
  
  
  
  common_rownames <- intersect(rownames(spliced.meth.test), intersect_names)
  if(!setequal(common_rownames, intersect_names)){
    stop("rownames(test.list[[1]]) should contains all the names in the methylation.site.relative.velocity.result")
  }
  if (!is.null(control.list)){
    common_rownames <- intersect(rownames(spliced.meth.control), intersect_names)
    if(!setequal(common_rownames, intersect_names)){
      stop("rownames(control.list[[1]]) should contains all the names in the methylation.site.relative.velocity.result")
    }
  }
  
  
  test.list.intersect <- list()
  
  test.list.intersect[[1]] <- spliced.meth.test[intersect_names,]
  test.list.intersect[[2]] <- spliced.unmeth.test[intersect_names,]
  test.list.intersect[[3]] <- unspliced.meth.test[intersect_names,]
  test.list.intersect[[4]] <- unspliced.unmeth.test[intersect_names,]
  
  exp_matrix <- test.list.intersect[[1]]+test.list.intersect[[2]]
  
  
  
  
  
  if (methylation_level=="M"){
    if(methylation_level_scale=="log2"){
      methylation_level_matrix <- (test.list.intersect[[1]]+test.list.intersect[[3]]+epsilon_M)/(test.list.intersect[[2]]+test.list.intersect[[4]]+epsilon_M)
      methylation_level_matrix <- log2(methylation_level_matrix)
    } else if (methylation_level_scale=="linear") {
      methylation_level_matrix <- (test.list.intersect[[1]]+test.list.intersect[[3]]+epsilon_M)/(test.list.intersect[[2]]+test.list.intersect[[4]]+epsilon_M)
      print("do not recommend use linear scale for M value")
    } else {
      stop(" methylation level scale should be linear or log2")
    }
    
  } else if (methylation_level=="Beta"){
    
    if(methylation_level_scale=="log2"){
      methylation_level_matrix <- (test.list.intersect[[1]]+test.list.intersect[[3]])/
        (test.list.intersect[[2]]+test.list.intersect[[4]]+test.list.intersect[[1]]+test.list.intersect[[3]]+epsilon_Beta)
      methylation_level_matrix <- log2(methylation_level_matrix)
      print("do not recommend use log2 scale for Beta value")
    } else if (methylation_level_scale=="linear") {
      methylation_level_matrix <- (test.list.intersect[[1]]+test.list.intersect[[3]])/
        (test.list.intersect[[2]]+test.list.intersect[[4]]+test.list.intersect[[1]]+test.list.intersect[[3]]+epsilon_Beta)
    } else {
      stop(" methylation level scale should be linear or log2")
    }
  } else if (methylation_level=="meta_M_sep"){
    
    if(methylation_level_scale=="log2"){
      methylation_level_matrix <- log2((test.list.intersect[[1]]+epsilon_M)/(test.list.intersect[[2]]+epsilon_M))+log2((test.list.intersect[[3]]+epsilon_M)/(test.list.intersect[[4]]+epsilon_M))
    } else if (methylation_level_scale=="linear") {
      methylation_level_matrix <- (test.list.intersect[[1]]+epsilon_M)/(test.list.intersect[[2]]+epsilon_M)+(test.list.intersect[[3]]+epsilon_M)/(test.list.intersect[[4]]+epsilon_M)
      print("do not recommend use linear scale for meta_M_sep value")
    } else {
      stop(" methylation level scale should be linear or log2")
    }
    
  } else if (methylation_level=="meta_Beta_sep"){
    
    
    if(methylation_level_scale=="log2"){
      methylation_level_matrix <- log2(test.list.intersect[[1]]/(test.list.intersect[[1]]+test.list.intersect[[2]]+epsilon_Beta))+
        log2(test.list.intersect[[3]]/(test.list.intersect[[3]]+test.list.intersect[[4]]+epsilon_Beta))
      print("do not recommend use log2 scale for meta_Beta_sep value")
    } else if (methylation_level_scale=="linear") {
      methylation_level_matrix <- test.list.intersect[[1]]/(test.list.intersect[[1]]+test.list.intersect[[2]]+epsilon_Beta)+
        test.list.intersect[[3]]/(test.list.intersect[[3]]+test.list.intersect[[4]]+epsilon_Beta)
    } else {
      stop(" methylation level scale should be linear or log2")
    }
    
    
  } else if (methylation_level=="RR"){
    
    if(methylation_level_scale=="log2"){
      methylation_level_matrix <- (test.list.intersect[[1]]+test.list.intersect[[3]]+epsilon_RR)/
        (test.list.intersect[[2]]+test.list.intersect[[4]]+test.list.intersect[[1]]+test.list.intersect[[3]]+epsilon_RR)
      methylation_level_matrix <- log2(methylation_level_matrix/control_information_intersect_all [["control_information_Beta"]])
      
    } else if (methylation_level_scale=="linear") {
      methylation_level_matrix <- (test.list.intersect[[1]]+test.list.intersect[[3]]+epsilon_RR)/
        (test.list.intersect[[2]]+test.list.intersect[[4]]+test.list.intersect[[1]]+test.list.intersect[[3]]+epsilon_RR)
      methylation_level_matrix <- (methylation_level_matrix/control_information_intersect_all [["_Beta"]])
      print("do not recommend use linear scale for RR value")
    } else {
      stop(" methylation level scale should be linear or log2")
    }
    
  } else if (methylation_level=="meta_RR_sep"){
    
    
    if(methylation_level_scale=="log2"){
      methylation_level_matrix_1 <- (test.list.intersect[[1]]+epsilon_RR)/
        (test.list.intersect[[1]]+test.list.intersect[[2]]+epsilon_RR)
      methylation_level_matrix_2 <- (test.list.intersect[[3]]+epsilon_RR)/
        (test.list.intersect[[3]]+test.list.intersect[[4]]+epsilon_RR)
      methylation_level_matrix <- log2(methylation_level_matrix_1/control_information_intersect_spliced [["control_information_Beta"]])+
        log2(methylation_level_matrix_2/control_information_intersect_unspliced [["control_information_Beta"]])
      
    } else if (methylation_level_scale=="linear") {
      methylation_level_matrix_1 <- (test.list.intersect[[1]]+epsilon_RR)/
        (test.list.intersect[[1]]+test.list.intersect[[2]]+epsilon_RR)
      methylation_level_matrix_2 <- (test.list.intersect[[3]]+epsilon_RR)/
        (test.list.intersect[[3]]+test.list.intersect[[4]]+epsilon_RR)
      methylation_level_matrix <- (methylation_level_matrix_1/control_information_intersect_spliced [["control_information_Beta"]])+
        (methylation_level_matrix_2/control_information_intersect_unspliced [["control_information_Beta"]])
      print("do not recommend use linear scale for meta_RR_sep value")
    } else {
      stop(" methylation level scale should be linear or log2")
    }
    
    
  } else if (methylation_level=="OR"){
    if(methylation_level_scale=="log2"){
      methylation_level_matrix <- (test.list.intersect[[1]]+test.list.intersect[[3]]+epsilon_OR)/
        (test.list.intersect[[2]]+test.list.intersect[[4]]+epsilon_OR)
      methylation_level_matrix <- log2(methylation_level_matrix/control_information_intersect_all [["control_information_M"]])
      
    } else if (methylation_level_scale=="linear") {
      methylation_level_matrix <- (test.list.intersect[[1]]+test.list.intersect[[3]]+epsilon_OR)/
        (test.list.intersect[[2]]+test.list.intersect[[4]]+epsilon_OR)
      methylation_level_matrix <- (methylation_level_matrix/control_information_intersect_all [["control_information_M"]])
      print("do not recommend use linear scale for OR value")
    } else {
      stop(" methylation level scale should be linear or log2")
    }
  } else if (methylation_level=="meta_OR_sep"){
    
    if(methylation_level_scale=="log2"){
      methylation_level_matrix_1 <- (test.list.intersect[[1]]+epsilon_OR)/
        (test.list.intersect[[2]]+epsilon_OR)
      methylation_level_matrix_2 <- (test.list.intersect[[3]]+epsilon_OR)/
        (test.list.intersect[[4]]+epsilon_OR)
      methylation_level_matrix <- log2(methylation_level_matrix_1/control_information_intersect_spliced [["control_information_M"]])+
        log2(methylation_level_matrix_2/control_information_intersect_unspliced [["control_information_M"]])
      
    } else if (methylation_level_scale=="linear") {
      methylation_level_matrix_1 <- (test.list.intersect[[1]]+epsilon_OR)/
        (test.list.intersect[[2]]+epsilon_OR)
      methylation_level_matrix_2 <- (test.list.intersect[[3]]+epsilon_OR)/
        (test.list.intersect[[4]]+epsilon_OR)
      methylation_level_matrix <- (methylation_level_matrix_1/control_information_intersect_spliced [["control_information_M"]])+
        (methylation_level_matrix_2/control_information_intersect_unspliced [["control_information_M"]])
      print("do not recommend use linear scale for meta_OR_sep value")
    } else {
      stop(" methylation level scale should be linear or log2")
    }
    
  } else if (methylation_level=="TCR"){
    
    if(methylation_level_scale=="log2"){
      methylation_level_matrix <- (test.list.intersect[[1]]+test.list.intersect[[3]]+epsilon_TCR)/
        (test.list.intersect[[2]]+test.list.intersect[[4]]+test.list.intersect[[1]]+test.list.intersect[[3]]+epsilon_TCR)
      methylation_level_matrix <- methylation_level_matrix-control_information_intersect_all [["control_information_Beta"]]
      methylation_level_matrix [methylation_level_matrix<=0] <- 0
      methylation_level_matrix <- log2(methylation_level_matrix)
      print("do not recommend use log2 scale for TCR value")
    } else if (methylation_level_scale=="linear") {
      methylation_level_matrix <- (test.list.intersect[[1]]+test.list.intersect[[3]]+epsilon_TCR)/
        (test.list.intersect[[2]]+test.list.intersect[[4]]+test.list.intersect[[1]]+test.list.intersect[[3]]+epsilon_TCR)
      methylation_level_matrix <- methylation_level_matrix-control_information_intersect_all [["control_information_Beta"]]
      methylation_level_matrix [methylation_level_matrix<=0] <- 0
      methylation_level_matrix <- log2(methylation_level_matrix)
    } else {
      stop(" methylation level scale should be linear or log2")
    }
    
  } else if (methylation_level=="meta_TCR_sep"){
    
    
    if(methylation_level_scale=="log2"){
      methylation_level_matrix_1 <- (test.list.intersect[[1]]+epsilon_TCR)/
        (test.list.intersect[[2]]+test.list.intersect[[1]]+epsilon_TCR)
      methylation_level_matrix_2 <- (test.list.intersect[[3]]+epsilon_TCR)/
        (test.list.intersect[[3]]+test.list.intersect[[4]]+epsilon_TCR)
      methylation_level_matrix_1 <- methylation_level_matrix_1-control_information_intersect_spliced [["control_information_Beta"]]
      methylation_level_matrix_1 [methylation_level_matrix_1<=0] <- 0
      methylation_level_matrix_2 <- methylation_level_matrix_2-control_information_intersect_unspliced [["control_information_Beta"]]
      methylation_level_matrix_2 [methylation_level_matrix_2<=0] <- 0
      methylation_level_matrix <- log2(methylation_level_matrix_1)+log2(methylation_level_matrix_2)
      print("do not recommend use log2 scale for meta_TCR_sep value")
    } else if (methylation_level_scale=="linear") {
      methylation_level_matrix_1 <- (test.list.intersect[[1]]+epsilon_TCR)/
        (test.list.intersect[[2]]+test.list.intersect[[1]]+epsilon_TCR)
      methylation_level_matrix_2 <- (test.list.intersect[[3]]+epsilon_TCR)/
        (test.list.intersect[[3]]+test.list.intersect[[4]]+epsilon_TCR)
      methylation_level_matrix_1 <- methylation_level_matrix_1-control_information_intersect_spliced [["control_information_Beta"]]
      methylation_level_matrix_1 [methylation_level_matrix_1<=0] <- 0
      methylation_level_matrix_2 <- methylation_level_matrix_2-control_information_intersect_unspliced [["control_information_Beta"]]
      methylation_level_matrix_2 [methylation_level_matrix_2<=0] <- 0
      methylation_level_matrix <- (methylation_level_matrix_1)+(methylation_level_matrix_2)
    } else {
      stop(" methylation level scale should be linear or log2")
    }
    
    
  } else {
    stop("methylation_level should be c('M','Beta','RR','OR','TCR','meta_M_sep','meta_Beta_sep','meta_RR_sep','meta_OR_sep','meta_TCR_sep')")
  }
  if (methylation_level%in%c("M","meta_M_sep")){
    if(methylation_level_scale=="log2"){
      methylation_level_matrix_spliced <- log2((test.list.intersect[[1]]+epsilon_M)/(test.list.intersect[[2]]+epsilon_M))
      methylation_level_matrix_unspliced <- log2((test.list.intersect[[3]]+epsilon_M)/(test.list.intersect[[4]]+epsilon_M))
    } else if (methylation_level_scale=="linear") {
      methylation_level_matrix_spliced <- ((test.list.intersect[[1]]+epsilon_M)/(test.list.intersect[[2]]+epsilon_M))
      methylation_level_matrix_unspliced <- ((test.list.intersect[[3]]+epsilon_M)/(test.list.intersect[[4]]+epsilon_M))
      print("do not recommend use linear scale for M or meta_M_sep value")
    } else {
      stop(" methylation level scale should be linear or log2")
    }
  } else if (methylation_level%in%c("Beta","meta_Beta_sep")){
    if(methylation_level_scale=="log2"){
      methylation_level_matrix_spliced <- log2(test.list.intersect[[1]]/(test.list.intersect[[1]]+test.list.intersect[[2]]+epsilon_Beta))
      methylation_level_matrix_unspliced <- log2(test.list.intersect[[3]]/(test.list.intersect[[3]]+test.list.intersect[[4]]+epsilon_Beta))
      print("do not recommend use log2 scale for Beta or meta_Beta_sep value")
    } else if (methylation_level_scale=="linear") {
      methylation_level_matrix_spliced <- (test.list.intersect[[1]]/(test.list.intersect[[1]]+test.list.intersect[[2]]+epsilon_Beta))
      methylation_level_matrix_unspliced <- (test.list.intersect[[3]]/(test.list.intersect[[3]]+test.list.intersect[[4]]+epsilon_Beta))
    } else {
      stop(" methylation level scale should be linear or log2")
    }
  } else if (methylation_level%in%c("RR","meta_RR_sep")){
    if(methylation_level_scale=="log2"){
      methylation_level_matrix_1 <- (test.list.intersect[[1]]+epsilon_RR)/
        (test.list.intersect[[1]]+test.list.intersect[[2]]+epsilon_RR)
      methylation_level_matrix_2 <- (test.list.intersect[[3]]+epsilon_RR)/
        (test.list.intersect[[3]]+test.list.intersect[[4]]+epsilon_RR)
      methylation_level_matrix_spliced <- log2(methylation_level_matrix_1/control_information_intersect_spliced [["control_information_Beta"]])
      methylation_level_matrix_unspliced <- log2(methylation_level_matrix_2/control_information_intersect_unspliced [["control_information_Beta"]])
      
    } else if (methylation_level_scale=="linear") {
      methylation_level_matrix_1 <- (test.list.intersect[[1]]+epsilon_RR)/
        (test.list.intersect[[1]]+test.list.intersect[[2]]+epsilon_RR)
      methylation_level_matrix_2 <- (test.list.intersect[[3]]+epsilon_RR)/
        (test.list.intersect[[3]]+test.list.intersect[[4]]+epsilon_RR)
      methylation_level_matrix_spliced <-  (methylation_level_matrix_1/control_information_intersect_spliced [["control_information_Beta"]])
      methylation_level_matrix_unspliced <-  (methylation_level_matrix_2/control_information_intersect_unspliced [["control_information_Beta"]])
      print("do not recommend use linear scale for RR or meta_RR_sep value")
    } else {
      stop(" methylation level scale should be linear or log2")
    }
  } else if (methylation_level%in%c("OR","meta_OR_sep")){
    if(methylation_level_scale=="log2"){
      methylation_level_matrix_1 <- (test.list.intersect[[1]]+epsilon_OR)/
        (test.list.intersect[[2]]+epsilon_OR)
      methylation_level_matrix_2 <- (test.list.intersect[[3]]+epsilon_OR)/
        (test.list.intersect[[4]]+epsilon_OR)
      methylation_level_matrix_spliced <- log2(methylation_level_matrix_1/control_information_intersect_spliced [["control_information_M"]])
      methylation_level_matrix_unspliced <- log2(methylation_level_matrix_2/control_information_intersect_unspliced [["control_information_M"]])
      
    } else if (methylation_level_scale=="linear") {
      methylation_level_matrix_1 <- (test.list.intersect[[1]]+epsilon_OR)/
        (test.list.intersect[[2]]+epsilon_OR)
      methylation_level_matrix_2 <- (test.list.intersect[[3]]+epsilon_OR)/
        (test.list.intersect[[4]]+epsilon_OR)
      methylation_level_matrix_spliced <-  (methylation_level_matrix_1/control_information_intersect_spliced [["control_information_M"]])
      methylation_level_matrix_unspliced <-  (methylation_level_matrix_2/control_information_intersect_unspliced [["control_information_M"]])
      print("do not recommend use linear scale for OR or meta_OR_sep value")
    } else {
      stop(" methylation level scale should be linear or log2")
    }
  } else if (methylation_level%in%c("TCR","meta_TCR_sep")){
    if(methylation_level_scale=="log2"){
      methylation_level_matrix_1 <- (test.list.intersect[[1]]+epsilon_TCR)/
        (test.list.intersect[[2]]+test.list.intersect[[1]]+epsilon_TCR)
      methylation_level_matrix_2 <- (test.list.intersect[[3]]+epsilon_TCR)/
        (test.list.intersect[[3]]+test.list.intersect[[4]]+epsilon_TCR)
      methylation_level_matrix_1 <- methylation_level_matrix_1-control_information_intersect_spliced [["control_information_Beta"]]
      methylation_level_matrix_1 [methylation_level_matrix_1<=0] <- 0
      methylation_level_matrix_2 <- methylation_level_matrix_2-control_information_intersect_unspliced [["control_information_Beta"]]
      methylation_level_matrix_2 [methylation_level_matrix_2<=0] <- 0
      methylation_level_matrix_spliced <- log2(methylation_level_matrix_1)
      methylation_level_matrix_unspliced <- log2(methylation_level_matrix_2)
      print("do not recommend use log2 scale for TCR or meta_TCR_sep value")
    } else if (methylation_level_scale=="linear") {
      methylation_level_matrix_1 <- (test.list.intersect[[1]]+epsilon_TCR)/
        (test.list.intersect[[2]]+test.list.intersect[[1]]+epsilon_TCR)
      methylation_level_matrix_2 <- (test.list.intersect[[3]]+epsilon_TCR)/
        (test.list.intersect[[3]]+test.list.intersect[[4]]+epsilon_TCR)
      methylation_level_matrix_1 <- methylation_level_matrix_1-control_information_intersect_spliced [["control_information_Beta"]]
      methylation_level_matrix_1 [methylation_level_matrix_1<=0] <- 0
      methylation_level_matrix_2 <- methylation_level_matrix_2-control_information_intersect_unspliced [["control_information_Beta"]]
      methylation_level_matrix_2 [methylation_level_matrix_2<=0] <- 0
      methylation_level_matrix_spliced <-  (methylation_level_matrix_1)
      methylation_level_matrix_unspliced <-  (methylation_level_matrix_2)
    } else {
      stop(" methylation level scale should be linear or log2")
    }
  }
  
  # exp_and_methylation_cor <- list()
  # for (i in 1: length(intersect_names)){
  #   exp_and_methylation_cor[[i]] <- cor(methylation_level_matrix[i,],exp_matrix[i,],method=cor_method)
  # }
  exp_and_methylation_cor <- mapply(function(meth, exp) {
    cor(meth, exp, method = cor_method)
  }, as.data.frame(t(as.matrix(methylation_level_matrix[intersect_names, ]))), 
  as.data.frame(t(as.matrix(exp_matrix[intersect_names, ]))))
  
  delta_exp_and_methylation_cor <- mapply(function(meth, exp) {
    cor(meth, exp, method = cor_method)
  }, as.data.frame(t(as.matrix(methylation_level_matrix[intersect_names, ]))), 
  as.data.frame(t(as.matrix(delta_exp_matrix[intersect_names, ]))))
  
  
   # exp_and_methylation_cor_lm <- list()
   # for (i in 1: length(intersect_names)){
   #   data.frame_ <- as.data.frame(cbind(methylation_level_matrix_spliced[i,],
   #                                      methylation_level_matrix_unspliced[i,],exp_matrix[i,]))
   #   colnames(data.frame_) <- c("methylation_level_spliced","methylation_level_unspliced","expression_level")
   #   exp_and_methylation_cor_lm[[i]] <- lm(expression_level~methylation_level_spliced+methylation_level_unspliced,
   #                                         data=data.frame_,na.action = na.action.lm)
   # }
  exp_and_methylation_cor_lm <- lapply(seq_along(intersect_names), function(i) {
    data.frame_ <- data.frame(
      methylation_level_spliced = methylation_level_matrix_spliced[i, ],
      methylation_level_unspliced = methylation_level_matrix_unspliced[i, ],
      expression_level = exp_matrix[i, ]
    )
    
    lm(expression_level ~ methylation_level_spliced + methylation_level_unspliced, 
       data = data.frame_, na.action = na.action.lm)
  })
  
  delta_exp_and_methylation_cor_lm <- lapply(seq_along(intersect_names), function(i) {
    data.frame_ <- data.frame(
      methylation_level_spliced = methylation_level_matrix_spliced[i, ],
      methylation_level_unspliced = methylation_level_matrix_unspliced[i, ],
      delta_expression_level = delta_exp_matrix[i, ]
    )
    
    lm(delta_expression_level ~ methylation_level_spliced + methylation_level_unspliced, 
       data = data.frame_, na.action = na.action.lm)
  })
  list_return <- list()
  list_return[["exp_and_methylation_cor"]] <- exp_and_methylation_cor
  list_return[["delta_exp_and_methylation_cor"]] <- delta_exp_and_methylation_cor
  list_return[["exp_and_methylation_cor_lm"]] <- exp_and_methylation_cor_lm
  list_return[["delta_exp_and_methylation_cor_lm"]] <- delta_exp_and_methylation_cor_lm
  other_information <- list()
  other_information[["delta_exp"]] <- delta_exp
  other_information[["methylation_level"]] <- methylation_level
  other_information[["methylation_level_scale"]] <- methylation_level_scale
  other_information[["cor_method"]] <- cor_method
  other_information[["na.action.lm"]] <- na.action.lm
  other_information[["methylation_level_matrix"]] <- methylation_level_matrix
  other_information[["methylation_level_matrix_spliced"]] <- methylation_level_matrix_spliced
  other_information[["methylation_level_matrix_unspliced"]] <- methylation_level_matrix_unspliced
  other_information[["exp_matrix"]] <- exp_matrix
  other_information[["delta_exp_matrix"]] <- delta_exp_matrix
  if(is.null(control.list)){
    control_information_intersect_spliced<- control_information_intersect_unspliced<- control_information_intersect_all <- NA
  }
  
  other_information[["control_information_intersect_spliced"]] <- control_information_intersect_spliced
  other_information[["control_information_intersect_unspliced"]] <- control_information_intersect_unspliced
  other_information[["control_information_intersect_all"]] <- control_information_intersect_all
  list_return[["other_information"]] <- other_information
  return(list_return)
}

