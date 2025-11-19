## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(VeloRM)

## ----eval=FALSE---------------------------------------------------------------
#  library(DESeq2)
#  
#  # Read CSV file containing read counts
#  read_counts <- read.csv("/gpfs/work/bio/haozhewang17/data/GSE29714/SNP_37199_read_counts.csv", row.names = NULL, check.names = FALSE)
#  # Extract pre-mRNA and spliced-mRNA columns
#  pre_counts     <- read_counts[, c("SRR494613_pre-mRNA", "SRR494614_pre-mRNA", "SRR494615_pre-mRNA", "SRR494616_pre-mRNA")]
#  spliced_counts <- read_counts[, c("SRR494613_mRNA_spliced", "SRR494614_mRNA_spliced", "SRR494615_mRNA_spliced", "SRR494616_mRNA_spliced")]
#  
#  # Define control/test and methylated/unmethylated groups
#  # Control samples: SRR494613, SRR494615; Test samples: SRR494614, SRR494616
#  meth_control   <- spliced_counts[, c("SRR494614_mRNA_spliced", "SRR494616_mRNA_spliced")]  # methylated counts for control?
#  meth_test      <- pre_counts[, c("SRR494614_pre-mRNA", "SRR494616_pre-mRNA")]               # methylated counts for test
#  unmeth_control <- spliced_counts[, c("SRR494613_mRNA_spliced", "SRR494615_mRNA_spliced")]  # unmethylated counts for control
#  unmeth_test    <- pre_counts[, c("SRR494613_pre-mRNA", "SRR494615_pre-mRNA")]               # unmethylated counts for test
#  
#  # Custom DESeq2 analysis function
#  DESeq2_test <- function(unmeth_control, unmeth_test, meth_control, meth_test) {
#  
#    # Convert all inputs to data.frame
#    meth_control   <- data.frame(meth_control)
#    meth_test      <- data.frame(meth_test)
#    unmeth_control <- data.frame(unmeth_control)
#    unmeth_test    <- data.frame(unmeth_test)
#  
#    # Function to calculate custom size factors
#    sizeFactor <- function(data) {
#      data <- as.matrix(data)
#      log_data <- log(data)
#      log_data[is.infinite(log_data)] <- NA
#      log_mean <- rowMeans(log_data, na.rm = TRUE)
#      log_s <- log_data - log_mean
#      s_size <- exp(apply(log_s, 2, function(x) median(x, na.rm = TRUE)))
#      return(s_size)
#    }
#  
#    # Calculate total counts for size factor computation
#    spliced_total   <- meth_control + unmeth_control
#    unspliced_total <- meth_test + unmeth_test
#    size_factor <- sizeFactor(cbind(spliced_total, unspliced_total))
#  
#    num_control <- ncol(unmeth_control)
#    num_test    <- ncol(unmeth_test)
#  
#    # Combine all four types of counts
#    counts <- cbind(unmeth_control, unmeth_test, meth_control, meth_test)
#  
#    # Create design matrix for DESeq2
#    design <- data.frame(
#      Trt = rep(c(rep("control", num_control), rep("test", num_test)), 2),
#      Meth = c(rep("methnon", num_control + num_test), rep("methylation", num_control + num_test))
#    )
#  
#    # Define model formula
#    model <- ~ Meth + Trt + Meth:Trt
#  
#    # Create DESeqDataSet
#    dds <- DESeqDataSetFromMatrix(countData = counts,
#                                  colData = design,
#                                  design = model)
#  
#    # Set custom size factors
#    sizeFactors(dds) <- rep(size_factor, 2)
#  
#    # Run DESeq2 using Wald test
#    dds <- DESeq(dds, test = "Wald")
#  
#    # Extract results for the interaction term
#    res <- DESeq2::results(dds, name = "Methmethylation.Trttest", cooksCutoff = FALSE)
#  
#    return(res)
#  }
#  
#  # Run the DESeq2 analysis
#  result <- DESeq2_test(unmeth_control, unmeth_test, meth_control, meth_test)
#  
#  # Save the results as an RDS file
#  saveRDS(result, "/gpfs/work/bio/haozhewang17/data/GSE29714/res_DESeq2.rds")
#  
#  

## ----eval=FALSE---------------------------------------------------------------
#  # -----------------------------
#  # 1. Read input CSV
#  # -----------------------------
#  read_counts <- read.csv(
#    "/gpfs/work/bio/haozhewang17/data/GSE29714/SNP_37199_read_counts.csv",
#    row.names = NULL,
#    check.names = FALSE
#  )
#  
#  # Extract pre-mRNA and spliced-mRNA columns
#  pre_counts     <- read_counts[, c("SRR494613_pre-mRNA", "SRR494614_pre-mRNA", "SRR494615_pre-mRNA", "SRR494616_pre-mRNA")]
#  spliced_counts <- read_counts[, c("SRR494613_mRNA_spliced", "SRR494614_mRNA_spliced", "SRR494615_mRNA_spliced", "SRR494616_mRNA_spliced")]
#  
#  # Define groups (exact same logic as your DESeq2 code)
#  meth_control   <- spliced_counts[, c("SRR494614_mRNA_spliced", "SRR494616_mRNA_spliced")]
#  meth_test      <- pre_counts[,     c("SRR494614_pre-mRNA",     "SRR494616_pre-mRNA")]
#  unmeth_control <- spliced_counts[, c("SRR494613_mRNA_spliced", "SRR494615_mRNA_spliced")]
#  unmeth_test    <- pre_counts[,     c("SRR494613_pre-mRNA",     "SRR494615_pre-mRNA")]
#  
#  # -----------------------------
#  # 2. Fisher test function
#  # -----------------------------
#  Fisher_test <- function(unmeth_control, unmeth_test, meth_control, meth_test) {
#  
#    # Sum replicates
#    Uc <- rowSums(unmeth_control)  # control unmethylated
#    Mc <- rowSums(meth_control)    # control methylated
#    Ut <- rowSums(unmeth_test)     # test unmethylated
#    Mt <- rowSums(meth_test)       # test methylated
#  
#    n <- length(Uc)
#  
#    pvals <- numeric(n)
#    odds  <- numeric(n)
#  
#    for (i in seq_len(n)) {
#  
#      tbl <- matrix(c(Uc[i], Mc[i],
#                      Ut[i], Mt[i]),
#                    nrow = 2, byrow = TRUE)
#  
#      ft <- fisher.test(tbl)
#  
#      pvals[i] <- ft$p.value
#      odds[i]  <- ft$estimate
#    }
#  
#    # Adjust p-values (BH-FDR)
#    padj <- p.adjust(pvals, method = "BH")
#  
#    return(data.frame(
#      Uc = Uc, Mc = Mc,
#      Ut = Ut, Mt = Mt,
#      odds_ratio = odds,
#      pvalue = pvals,
#      padj = padj
#    ))
#  }
#  
#  # -----------------------------
#  # 3. Run Fisher test
#  # -----------------------------
#  fisher_result <- Fisher_test(unmeth_control, unmeth_test, meth_control, meth_test)
#  
#  # -----------------------------
#  # 4. Save output
#  # -----------------------------
#  saveRDS(fisher_result, "/gpfs/work/bio/haozhewang17/data/GSE29714/res_Fisher.rds")

## ----warning=FALSE------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)

res_DESeq2 <- readRDS("D:/data/VeloRM_data/MeRIP-seq/res_DESeq2.rds")

log2FC <-res_DESeq2@listData[["log2FoldChange"]]
pvals <- res_DESeq2@listData[["pvalue"]]

plot_data <- data.frame(log2FC = log2FC, 
                        pvalue = pvals,
                        neg_log10_p = -log10(pvals))


volcano_plot_DESeq2 <- ggplot(plot_data, aes(x = log2FC, y = neg_log10_p)) +
  # Add light gray background for -1 < log2FC < 1 (non-significant FC region)
  annotate("rect", xmin = -1, xmax = 1, ymin = 0, ymax = 30, 
           fill = "gray90", alpha = 0.3) +
  
  # Color points based on significance and log2FC:
  # - Red: log2FC > 1 & p < 0.05 (significantly up-regulated)
  # - Blue: log2FC < -1 & p < 0.05 (significantly down-regulated)
  # - Gray: Not significant or |log2FC| â‰¤ 1
  geom_point(aes(color = case_when(
    log2FC > 1 & pvalue < 0.05 ~ "Unspliced",
    log2FC < -1 & pvalue < 0.05 ~ "Spliced",
    TRUE ~ "Not significant"
  )), alpha = 0.3) +
  
  # Manual color scale (red, blue, gray)
  scale_color_manual(
    values = c(
      "Spliced" = "red",
      "Unspliced" = "blue",
      "Not significant" = "gray50"
    ),
    name = "Regulation"  # Legend title
  ) +
  
  # Add dashed lines for thresholds:
  # - Vertical: log2FC = Â±1
  # - Horizontal: p-value = 0.05 (-log10(0.05) â‰ˆ 1.3)
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  
  # Axis labels and theme
  labs(
    x = "log2 Fold Change", 
    y = "-log10(p-value)",
    title = NULL
  ) +
  theme_grey() +
  # Set axis limits (zoomed-in view)
  coord_cartesian(ylim = c(0, 4))
volcano_plot_DESeq2

n_unspliced <- length(which(log2FC >= 1 & pvals <= 0.05))  
n_spliced <- length(which(log2FC <= -1 & pvals <= 0.05))  
n_nonsig <- 37199-n_unspliced-n_spliced

total <- n_unspliced + n_spliced + n_nonsig
Percentage <- round(c(n_unspliced, n_spliced, n_nonsig) / total * 100, 1)

df <- data.frame(
  Category = c("Unspliced", "Spliced", "Insignificant"),
  Count = c(n_unspliced, n_spliced, n_nonsig)
)



df$logCount <- log10(df$Count + 1)  

log_bar <- ggplot(df, aes(x = Category, y = logCount, fill = Category)) +
  geom_col(alpha = 0.7) +
  geom_text(aes(label = Count), vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("Unspliced"="blue","Spliced"="red","Insignificant"="gray50")) +
  labs(
    x = "",
    y = "log10(Count)",
    title = NULL
  ) +
  theme_grey()+coord_cartesian(ylim=c(0,5))

log_bar


## ----warning=FALSE------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)

res_Fisher <- readRDS("D:/data/VeloRM_data/MeRIP-seq/res_Fisher.rds")

log2FC <- log2(res_Fisher$odds_ratio)
pvals <- res_Fisher$pvalue

plot_data <- data.frame(log2FC = log2FC, 
                        pvalue = pvals,
                        neg_log10_p = -log10(pvals))
plot_data <- plot_data[is.finite(plot_data$log2FC) &
                         is.finite(plot_data$neg_log10_p), ]

volcano_plot_Fisher <- ggplot(plot_data, aes(x = log2FC, y = neg_log10_p)) +
  # Add light gray background for -1 < log2FC < 1 (non-significant FC region)
  annotate("rect", xmin = -1, xmax = 1, ymin = 0, ymax = 30, 
           fill = "gray90", alpha = 0.3) +
  
  # Color points based on significance and log2FC:
  # - Red: log2FC > 1 & p < 0.05 (significantly up-regulated)
  # - Blue: log2FC < -1 & p < 0.05 (significantly down-regulated)
  # - Gray: Not significant or |log2FC| ?? 1
  geom_point(aes(color = case_when(
    log2FC > 1 & pvalue < 0.05 ~ "Unspliced",
    log2FC < -1 & pvalue < 0.05 ~ "Spliced",
    TRUE ~ "Not significant"
  )), alpha = 0.3) +
  
  # Manual color scale (red, blue, gray)
  scale_color_manual(
    values = c(
      "Spliced" = "red",
      "Unspliced" = "blue",
      "Not significant" = "gray50"
    ),
    name = "Regulation"  # Legend title
  ) +
  
  # Add dashed lines for thresholds:
  # - Vertical: log2FC = ??1
  # - Horizontal: p-value = 0.05 (-log10(0.05) ?? 1.3)
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  
  # Axis labels and theme
  labs(
    x = "log2 Fold Change", 
    y = "-log10(p-value)",
    title = NULL
  ) +
  theme_grey() +
  # Set axis limits (zoomed-in view)
  coord_cartesian(ylim = c(0, 5),xlim = c(-5, 5))
volcano_plot_Fisher

n_unspliced <- length(which(log2FC >= 1 & pvals <= 0.05))  
n_spliced <- length(which(log2FC <= -1 & pvals <= 0.05))  
n_nonsig <- 37199-n_unspliced-n_spliced

total <- n_unspliced + n_spliced + n_nonsig
Percentage <- round(c(n_unspliced, n_spliced, n_nonsig) / total * 100, 1)

df <- data.frame(
  Category = c("Unspliced", "Spliced", "Insignificant"),
  Count = c(n_unspliced, n_spliced, n_nonsig)
)

df$logCount <- log10(df$Count + 1)  

log_bar <- ggplot(df, aes(x = Category, y = logCount, fill = Category)) +
  geom_col(alpha = 0.7) +
  geom_text(aes(label = Count), vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("Unspliced"="blue","Spliced"="red","Insignificant"="gray50")) +
  labs(
    x = "",
    y = "log10(Count)",
    title = NULL
  ) +
  theme_grey()+coord_cartesian(ylim=c(0,5))

log_bar

