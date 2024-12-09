#' Title
#'
#' @param test.list - a list of the test cells, with six data.frame of the methylation.spliced, unmethylation.spliced,
#'  methylation.unspliced, unmethylation.unspliced, methylation.ambiguous, unmethylation.ambiguous,read counts
#' @param control.list - a list of the test cells, with six data.frame of the methylation.spliced, unmethylation.spliced,
#'  methylation.unspliced, unmethylation.unspliced,methylation.ambiguous, unmethylation.ambiguous,read counts
#' @param ambiguous.as.spliced prob ambiguous from spliced 
#' @param epsilon_M epsilon to calculate the M
#' @param epsilon_Beta epsilon to calculate the Beta
#' @param epsilon_RR epsilon to calculate the RR
#' @param epsilon_OR epsilon to calculate the OR
#' 
#' @return a list(test.return,control.return,index.site1,index.site2)
#'  1. test.return: a list of the test cells, with four data.frame of the methylation.spliced, 
#'  unmethylation.spliced,methylation.unspliced, unmethylation.unspliced read counts
#'  2. control.return: a list of the control cells, with four data.frame of the methylation.spliced, 
#'  unmethylation.spliced,methylation.unspliced, unmethylation.unspliced read counts
#'  3. index.site: a list of information of the sites, the test.return and control return only contains the 
#'  methylation occurs at both spliced and unspliced RNA
#' @export
#'
methylation.sites.preprocess <- function(test.list, control.list=NULL,ambiguous.as.spliced=NA,
                                         epsilon_M=1e-7,epsilon_Beta=1,epsilon_RR=1e-7,epsilon_OR=1e-7){
  # check length of list
  if (length(test.list)!=6)
    stop("the length of test list need to be 6")
  if (!is.null(control.list)&length(test.list)!=6)
    stop("the length of control list need to be 6")
  
  for (i in 1:6) {
    if (!is.null(test.list[[i]])){
      if (is.null(colnames(test.list[[i]]))|is.null(rownames(test.list[[i]])))
        stop("it need have colnames and rownames")
    }
    if(!is.null(control.list)){
      if (!is.null(control.list[[i]])){
        if (is.null(colnames(control.list[[i]]))|is.null(rownames(control.list[[i]])))
          stop("it need have colnames and rownames")
      }
    }
    
  }
  
  # check dim
  if (!all(colnames(test.list[[1]]) == colnames(test.list[[2]])))
    stop("emat and nmat must have the same columns (cells) for test cells")
  if (!all(colnames(test.list[[3]]) == colnames(test.list[[4]])))
    stop("emat and nmat must have the same columns (cells) for test cells")
  if (!all(colnames(test.list[[1]]) == colnames(test.list[[4]])))
    stop("emat and nmat must have the same columns (cells) for test cells")
  if (is.null(test.list[[5]])&is.null(test.list[[6]])){
    test.list[[5]] <- test.list[[6]] <- Matrix(0, dim(test.list[[1]])[1], dim(test.list[[1]])[2], sparse = TRUE)
    rownames(test.list[[5]]) <- rownames(test.list[[6]]) <- rownames(test.list[[1]])
    colnames(test.list[[5]]) <- colnames(test.list[[6]]) <- colnames(test.list[[1]])
  } else if((!is.null(test.list[[5]]))&(!is.null(test.list[[6]]))){
    if (!all(colnames(test.list[[5]]) == colnames(test.list[[6]])))
      stop("emat and nmat must have the same columns (cells) for test cells")
    if (!all(colnames(test.list[[1]]) == colnames(test.list[[5]])))
      stop("emat and nmat must have the same columns (cells) for test cells")
  } else{ 
    if(is.null(test.list[[5]])){
      if (!all(colnames(test.list[[1]]) == colnames(test.list[[6]])))
        stop("emat and nmat must have the same columns (cells) for test cells")
      test.list[[5]] <- Matrix(0, dim(test.list[[1]])[1], dim(test.list[[1]])[2], sparse = TRUE)
      rownames(test.list[[5]])  <- rownames(test.list[[1]])
      colnames(test.list[[5]])  <- colnames(test.list[[1]])
    }
    if(is.null(test.list[[6]])){
      if (!all(colnames(test.list[[1]]) == colnames(test.list[[5]])))
        stop("emat and nmat must have the same columns (cells) for test cells")
      test.list[[6]] <- Matrix(0, dim(test.list[[1]])[1], dim(test.list[[1]])[2], sparse = TRUE)
      rownames(test.list[[6]])  <- rownames(test.list[[1]])
      colnames(test.list[[6]])  <- colnames(test.list[[1]])
    }
  }
  if (length(colnames(test.list[[1]]))<=1){
    stop("length of the colnames need > 1")
  }
  
  
  # check dim
  if(!is.null(control.list)){
    if (!all(colnames(control.list[[1]]) == colnames(control.list[[2]])))
      stop("emat and nmat must have the same columns (cells) for control cells")
    if (!all(colnames(control.list[[3]]) == colnames(control.list[[4]])))
      stop("emat and nmat must have the same columns (cells) for control cells")
    if (!all(colnames(control.list[[1]]) == colnames(control.list[[4]])))
      stop("emat and nmat must have the same columns (cells) for control cells")
    if (is.null(control.list[[5]])&is.null(control.list[[6]])){
      control.list[[5]] <- control.list[[6]] <- Matrix(0, dim(control.list[[1]])[1], dim(control.list[[1]])[2], sparse = TRUE)
      rownames(control.list[[5]]) <- rownames(control.list[[6]]) <- rownames(control.list[[1]])
      colnames(control.list[[5]]) <- colnames(control.list[[6]]) <- colnames(control.list[[1]])
    } else if((!is.null(control.list[[5]]))&(!is.null(control.list[[6]]))){
      if (!all(colnames(control.list[[5]]) == colnames(control.list[[6]])))
        stop("emat and nmat must have the same columns (cells) for control cells")
      if (!all(colnames(control.list[[1]]) == colnames(control.list[[5]])))
        stop("emat and nmat must have the same columns (cells) for control cells")
    } else{
        if(is.null(control.list[[5]])){
          if (!all(colnames(control.list[[1]]) == colnames(control.list[[6]])))
            stop("emat and nmat must have the same columns (cells) for control cells")
          control.list[[5]] <- Matrix(0, dim(control.list[[1]])[1], dim(control.list[[1]])[2], sparse = TRUE)
          rownames(control.list[[5]])  <- rownames(control.list[[1]])
          colnames(control.list[[5]])  <- colnames(control.list[[1]])
        }
        if(is.null(control.list[[6]])){
          if (!all(colnames(control.list[[1]]) == colnames(control.list[[5]])))
            stop("emat and nmat must have the same columns (cells) for control cells")
          control.list[[6]] <- Matrix(0, dim(control.list[[1]])[1], dim(control.list[[1]])[2], sparse = TRUE)
          rownames(control.list[[6]])  <- rownames(control.list[[1]])
          colnames(control.list[[6]])  <- colnames(control.list[[1]])
        }
    }
    if (length(colnames(control.list[[1]]))<=1){
      stop("length of the colnames need > 1")
    }
  }
  
  
  vg1 <- intersect(rownames(test.list[[1]]),rownames(test.list[[2]]))
  vg2 <- intersect(rownames(test.list[[3]]),rownames(test.list[[4]]))
  vg3 <- intersect(rownames(test.list[[5]]),rownames(test.list[[6]]))
  vg.test <- intersect(intersect(vg1,vg2),vg3)
  
  vg1 <- intersect(rownames(control.list[[1]]),rownames(control.list[[2]]))
  vg2 <- intersect(rownames(control.list[[3]]),rownames(control.list[[4]]))
  vg3 <- intersect(rownames(control.list[[5]]),rownames(control.list[[6]]))
  vg.control <- intersect(intersect(vg1,vg2),vg3)
  
  if (!is.null(control.list)){
    if (!all(vg.test %in% vg.control)){
      stop("the control.list site information must larger or equal to the test.list")
    }
    
  }
  
  

  
  # filter the site by four class
  
  if (!is.na(ambiguous.as.spliced)){
    if (ambiguous.as.spliced>1|ambiguous.as.spliced<0){
      stop("ambiguous.as.spliced need from 0 to 1")
    }  
    spliced.meth.test.sum <- rowSums(test.list[[1]][vg.test,]+round(ambiguous.as.spliced*test.list[[5]][vg.test,])) 
    unspliced.meth.test.sum <- rowSums(test.list[[3]][vg.test,]+round((1-ambiguous.as.spliced)*test.list[[5]][vg.test,]))
    spliced.unmeth.test.sum <- rowSums(test.list[[2]][vg.test,]+round(ambiguous.as.spliced*test.list[[6]][vg.test,])) 
    unspliced.unmeth.test.sum <- rowSums(test.list[[4]][vg.test,]+round((1-ambiguous.as.spliced)*test.list[[6]][vg.test,]))
    if(!is.null(control.list)){
      spliced.meth.control.sum <- rowSums(control.list[[1]][vg.control,]+round(ambiguous.as.spliced*control.list[[5]][vg.control,])) 
      unspliced.meth.control.sum <- rowSums(control.list[[3]][vg.control,]+round((1-ambiguous.as.spliced)*control.list[[5]][vg.control,]))
      spliced.unmeth.control.sum <- rowSums(control.list[[2]][vg.control,]+round(ambiguous.as.spliced*control.list[[6]][vg.control,])) 
      unspliced.unmeth.control.sum <- rowSums(control.list[[4]][vg.control,]+round((1-ambiguous.as.spliced)*control.list[[6]][vg.control,]))
    }
    
  } else{
    spliced.meth.test.sum <- rowSums(test.list[[1]][vg.test,])
    unspliced.meth.test.sum <- rowSums(test.list[[3]][vg.test,])
    spliced.unmeth.test.sum <- rowSums(test.list[[2]][vg.test,])
    unspliced.unmeth.test.sum <- rowSums(test.list[[4]][vg.test,])
    if(!is.null(control.list)){
      spliced.meth.control.sum <- rowSums(control.list[[1]][vg.control,])
      unspliced.meth.control.sum <- rowSums(control.list[[3]][vg.control,])
      spliced.unmeth.control.sum <- rowSums(control.list[[2]][vg.control,])
      unspliced.unmeth.control.sum <- rowSums(control.list[[4]][vg.control,])
    }
  }
  
  # the site methylation only occurs at unspliced RNA
  index.1 <- vg.test [which(spliced.meth.test.sum==0&unspliced.meth.test.sum!=0)]
  # the site methylation only occurs at spliced RNA
  index.2 <- vg.test [which(spliced.meth.test.sum!=0&unspliced.meth.test.sum==0)]
  # the site without methylation 
  index.3 <- vg.test [which(spliced.meth.test.sum==0&unspliced.meth.test.sum==0)]
  # the site methylation occurs at both spliced and unspliced RNA
  index.4 <- vg.test [which(spliced.meth.test.sum!=0&unspliced.meth.test.sum!=0)]
  
  
  # filter the site by sub-four class
  
  if (!is.na(ambiguous.as.spliced)){
    spliced.unmeth.test.sum.index.4 <- rowSums(test.list[[2]][index.4,]+round(ambiguous.as.spliced*test.list[[6]][index.4,])) 
    unspliced.unmeth.test.sum.index.4 <- rowSums(test.list[[4]][index.4,]+round((1-ambiguous.as.spliced)*test.list[[6]][index.4,]))
  }else{
    spliced.unmeth.test.sum.index.4 <- rowSums(test.list[[2]][index.4,])
    unspliced.unmeth.test.sum.index.4 <- rowSums(test.list[[4]][index.4,])
  }
  
  
  
  # the site unmethylation only occurs at unspliced RNA
  index.5 <- index.4 [which(spliced.unmeth.test.sum.index.4==0&unspliced.unmeth.test.sum.index.4!=0)]
  # the site unmethylation only occurs at spliced RNA
  index.6 <- index.4 [which(spliced.unmeth.test.sum.index.4!=0&unspliced.unmeth.test.sum.index.4==0)]
  # the site without unmethylation 
  index.7 <- index.4 [which(spliced.unmeth.test.sum.index.4==0&unspliced.unmeth.test.sum.index.4==0)]
  # the site unmethylation occurs at both spliced and unspliced RNA
  index.8 <- index.4 [which(spliced.unmeth.test.sum.index.4!=0&unspliced.unmeth.test.sum.index.4!=0)]
  
  spliced_methylation_Beta_value_meta <- (spliced.meth.test.sum)/(spliced.unmeth.test.sum+spliced.meth.test.sum+epsilon_Beta)
  unspliced_methylation_Beta_value_meta <- (unspliced.meth.test.sum)/(unspliced.unmeth.test.sum+unspliced.meth.test.sum+epsilon_Beta)
  
  spliced_methylation_M_value_meta <- log2((spliced.meth.test.sum+epsilon_M)/(spliced.unmeth.test.sum+epsilon_M))
  unspliced_methylation_M_value_meta <- log2((unspliced.meth.test.sum+epsilon_M)/(unspliced.unmeth.test.sum+epsilon_M))
  
  if (!is.null(control.list)){
    
    spliced_methylation_Beta_value_meta_control <- (spliced.meth.control.sum)/(spliced.unmeth.control.sum+spliced.meth.control.sum+epsilon_Beta)
    unspliced_methylation_Beta_value_meta_control <- (unspliced.meth.control.sum)/(unspliced.unmeth.control.sum+unspliced.meth.control.sum+epsilon_Beta)
    
    spliced_methylation_M_value_meta_control <- log2((spliced.meth.control.sum+epsilon_M)/(spliced.unmeth.control.sum+epsilon_M))
    unspliced_methylation_M_value_meta_control <- log2((unspliced.meth.control.sum+epsilon_M)/(unspliced.unmeth.control.sum+epsilon_M))
    
    spliced_methylation_RR_value_meta <- (spliced_methylation_Beta_value_meta+epsilon_RR)/
      (spliced_methylation_Beta_value_meta_control+epsilon_RR)
    spliced_methylation_OR_value_meta <- (2^spliced_methylation_M_value_meta+epsilon_RR)/
      (2^spliced_methylation_M_value_meta_control+epsilon_OR)
    spliced_methylation_TCR_value_meta <- spliced_methylation_Beta_value_meta-spliced_methylation_Beta_value_meta_control
    spliced_methylation_TCR_value_meta [spliced_methylation_TCR_value_meta<0] <- 0
    
    unspliced_methylation_RR_value_meta <- (unspliced_methylation_Beta_value_meta+epsilon_RR)/
      (unspliced_methylation_Beta_value_meta_control+epsilon_RR)
    unspliced_methylation_OR_value_meta <- (2^unspliced_methylation_M_value_meta+epsilon_RR)/
      (2^unspliced_methylation_M_value_meta_control+epsilon_OR)
    unspliced_methylation_TCR_value_meta <- unspliced_methylation_Beta_value_meta-unspliced_methylation_Beta_value_meta_control
    unspliced_methylation_TCR_value_meta [unspliced_methylation_TCR_value_meta<0] <- 0
    
  } else{
    spliced_methylation_RR_value_meta <- NULL
    spliced_methylation_OR_value_meta <- NULL
    spliced_methylation_TCR_value_meta <- NULL
    
    unspliced_methylation_RR_value_meta <- NULL
    unspliced_methylation_OR_value_meta <- NULL
    unspliced_methylation_TCR_value_meta <- NULL
  }
  
  
  
  
  if (!is.na(ambiguous.as.spliced)){
    spliced.meth.test <- test.list[[1]][index.8,]+round(ambiguous.as.spliced*test.list[[5]][index.8,])
    unspliced.meth.test <- test.list[[3]][index.8,]+round((1-ambiguous.as.spliced)*test.list[[5]][index.8,])
    spliced.unmeth.test <- test.list[[2]][index.8,]+round(ambiguous.as.spliced*test.list[[6]][index.8,])
    unspliced.unmeth.test <- test.list[[4]][index.8,]+round((1-ambiguous.as.spliced)*test.list[[6]][index.8,])
    if (is.null(control.list)){
      spliced.meth.control <- NULL
      unspliced.meth.control <- NULL
      spliced.unmeth.control <- NULL
      unspliced.unmeth.control <- NULL
    } else{
      spliced.meth.control <- control.list[[1]][index.8,]+round(ambiguous.as.spliced*control.list[[5]][index.8,])
      unspliced.meth.control <- control.list[[3]][index.8,]+round((1-ambiguous.as.spliced)*control.list[[5]][index.8,])
      spliced.unmeth.control <- control.list[[2]][index.8,]+round(ambiguous.as.spliced*control.list[[6]][index.8,])
      unspliced.unmeth.control <- control.list[[4]][index.8,]+round((1-ambiguous.as.spliced)*control.list[[6]][index.8,])
    }
    
  } else{
    spliced.meth.test <- test.list[[1]][index.8,]
    unspliced.meth.test <- test.list[[3]][index.8,]
    spliced.unmeth.test <- test.list[[2]][index.8,]
    unspliced.unmeth.test <- test.list[[4]][index.8,]
    if (is.null(control.list)){
      spliced.meth.control <- NULL
      unspliced.meth.control <- NULL
      spliced.unmeth.control <- NULL
      unspliced.unmeth.control <- NULL
    } else{
      spliced.meth.control <- control.list[[1]][index.8,]
      unspliced.meth.control <- control.list[[3]][index.8,]
      spliced.unmeth.control <- control.list[[2]][index.8,]
      unspliced.unmeth.control <- control.list[[4]][index.8,]
    }
  }
  
  
  test.return <- list()
  test.return[[1]] <-  spliced.meth.test
  test.return[[2]] <-  spliced.unmeth.test
  test.return[[3]] <-  unspliced.meth.test
  test.return[[4]] <-  unspliced.unmeth.test
  names(test.return) <- c("spliced.meth.test","spliced.unmeth.test","unspliced.meth.test","unspliced.unmeth.test")
  
  if (is.null(control.list)){
    control.return <- NULL
  } else{
    control.return <- list()
    control.return[[1]] <-  spliced.meth.control
    control.return[[2]] <-  spliced.unmeth.control
    control.return[[3]] <-  unspliced.meth.control
    control.return[[4]] <-  unspliced.unmeth.control
    names(control.return) <- c("spliced.meth.control","spliced.unmeth.control","unspliced.meth.control","unspliced.unmeth.control")
  }
  
  index.site <- list("methylation only occurs at unspliced RNA"=index.1,
                     "methylation only occurs at spliced RNA"=index.2,
                     "without methylation"=index.3,
                     "have methylation on both but unmethylation.spliced no read counts"=index.5,
                     "have methylation on both but unmethylation.unspliced no read counts"=index.6,
                     "have methylation on both but without unmethylation read counts"=index.7,
                     "have methylation on both and unmethylation on both"=index.8)
  
  
  
  res_list <- list(test.return, control.return, index.site,
                   spliced_methylation_Beta_value_meta,unspliced_methylation_Beta_value_meta,
                   spliced_methylation_M_value_meta,unspliced_methylation_M_value_meta,
                   spliced_methylation_RR_value_meta,unspliced_methylation_RR_value_meta,
                   spliced_methylation_OR_value_meta,unspliced_methylation_OR_value_meta,
                   spliced_methylation_TCR_value_meta,unspliced_methylation_TCR_value_meta,vg.test)
  names(res_list) <- c("test","control","index",
                       "spliced_methylation_Beta_value_meta","unspliced_methylation_Beta_value_meta",
                       "spliced_methylation_M_value_meta","unspliced_methylation_M_value_meta",
                       "spliced_methylation_RR_value_meta","unspliced_methylation_RR_value_meta",
                       "spliced_methylation_OR_value_meta","unspliced_methylation_OR_value_meta",
                       "spliced_methylation_TCR_value_meta","unspliced_methylation_TCR_value_meta",
                       "vg.test")
  return(res_list)
  
}