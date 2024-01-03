## Analyze RNA data


# Import utils and Setup --------------------------------------------------
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse)) 
suppressPackageStartupMessages(library(edgeR)) 

# Load Data ---------------------------------------------------------------
df_rna_pr_cr <- read.csv(file = "rna_counts_pr_cr.csv", header = T)
df_rna_sd_pd <- read.csv(file = "rna_counts_sd_pd.csv", header = T)

# Data Wrangling ----------------------------------------------------------
## Info
str(df_rna_pr_cr);dim(df_rna_pr_cr)
str(df_rna_sd_pd);dim(df_rna_sd_pd)

## Create a function to check if all values in a matrix are >= 0
check_0 <- function(data){
  if (all(data >= 0)) {
    print("All values are >= 0")
  } else {
    print("RNA-seq data cannot be negative")}
}
## Apply
check_0(df_rna_pr_cr)
check_0(df_rna_sd_pd)
## Check Na
table(is.na(df_rna_pr_cr)) ## No NA
table(is.na(df_rna_sd_pd)) ## No NA
## View
View(df_rna_pr_cr[1:10, 1:10])
View(df_rna_sd_pd[1:10, 1:10])
## Set row names
rownames(df_rna_pr_cr) <- df_rna_pr_cr$X
df_rna_pr_cr$X         <- NULL
rownames(df_rna_sd_pd) <- df_rna_sd_pd$X
df_rna_sd_pd$X         <- NULL
## Check
View(df_rna_pr_cr[1:10, 1:10])
View(df_rna_sd_pd[1:10, 1:10])
## Check if some patient is duplicate
length(unique(rownames(df_rna_pr_cr))) == nrow(df_rna_pr_cr)
length(unique(rownames(df_rna_sd_pd))) == nrow(df_rna_sd_pd)
## All the patients are repeated one time

## Make a control about values of patients
## Control if there are some patients with a total expression of 0
## Sum for each row
is_zero_row_pr_cr <- rowSums(df_rna_pr_cr) == ncol(df_rna_pr_cr)
is_zero_row_sd_pd <- rowSums(df_rna_sd_pd) == ncol(df_rna_sd_pd)
## Check
table(is_zero_row_pr_cr) ## No
table(is_zero_row_sd_pd) ## No
## The same control also for the genes
## Sum for each columns
is_zero_col_pr_cr <- colSums(df_rna_pr_cr) == nrow(df_rna_pr_cr)
is_zero_col_sd_pd <- colSums(df_rna_sd_pd) == nrow(df_rna_sd_pd)
## Check
table(is_zero_col_pr_cr) ## Yes, 67
table(is_zero_col_sd_pd) ## Yes, 40
## In this case, there are some genes that are never expressed and
## can be removed
df_rna_pr_cr <- df_rna_pr_cr[, colSums(df_rna_pr_cr) != 0]
df_rna_sd_pd <- df_rna_sd_pd[, colSums(df_rna_sd_pd) != 0]
## Check
dim(df_rna_pr_cr);dim(df_rna_sd_pd)
## Extract common genes
common_genes <- intersect(colnames(df_rna_pr_cr), colnames(df_rna_sd_pd))
length(common_genes)
## Final
final_df_rna_pr_cr <- df_rna_pr_cr[, common_genes]
final_df_rna_sd_pd <- df_rna_sd_pd[, common_genes]
## Check
dim(final_df_rna_pr_cr);dim(final_df_rna_sd_pd)
## Create an unique dataset
final_expr <- rbind(final_df_rna_pr_cr, final_df_rna_sd_pd)
## Check
dim(final_expr)
View(final_expr[1:10, 1:10])

# DEG Analysis ------------------------------------------------------------
## Create a new data with information about the status of the patients
## Exclude description column
description_gene                <- final_expr$Description
final_expr$Description          <- NULL
## Check
dim(final_expr)
## Create information about patients
patients_condition <- data.frame(Patients = colnames(final_expr),
                                 Condition = c(rep("cancer",
                                                   ncol(final_expr))))
patients_condition$Condition    <- factor(patients_condition$Condition)
rownames(patients_condition)    <- patients_condition$Patients
patients_condition$Patients     <- NULL

# DEG Analysis ------------------------------------------------------------
## It is absolutely critical that the columns of the count matrix and
## the rows of the column data (information about samples) are in the
## same order.
ncol(final_expr) == nrow(patients_condition)
## Creating DESeq2 object
DESeq_dt <- DESeqDataSetFromMatrix(countData = final_expr,
                                   colData = patients_condition,
                                   design = ~ 1)
## It was used design = ~1, because the condition for patients is
## unique ---> Cancer.
## Check
dim(DESeq_dt)
## 2457 rows and 152 columns

## Pre-filtering
## We apply a pre-filtering operation, selecting only the genes
## that have at least 10 reads.
## 1.removing rows in which there are very few reads
## 2.we reduce the memory size of the DESeq_dt data object
keep <- rowSums(counts(DESeq_dt)) >= 10
DESeq_dt_filtering <- DESeq_dt[keep, ]
DESeq_dt_filtering$BOR_cat <- relevel(DESeq_dt_filtering$BOR_cat,
                                        ref = "resistance")

dim(DESeq_dt_filtering)
## The number of observation is not reduced

## Identify DEGs specifying the thresholds setting using p-value and
## Fold Change FC >= |1.2|
## We apply Benjamini-Hochberg correction for multiple testing (FDR)
DESeq_output <- DESeq(DESeq_dt_filtering)
res          <- results(DESeq_output, alpha = 0.05,
                        lfcThreshold = 1.2, pAdjustMethod = "BH")
## See the results
res;summary(res)
## We obtain a set of 2457  genes
## LFC >  1.20  ---> Up-regulated genes   (2457, 100%)
## LFC <= -1.20 ---> Down regulated genes (0, 0%)

## Volcano Plot
volcano_dat <- as.data.frame(res)
volcano_dat <- volcano_dat %>% 
  mutate(Expression = case_when(log2FoldChange >= 1.2 & padj <= 0.05 ~ "Up-regulated",
                                log2FoldChange <= -1.2 & padj <= 0.05 ~ "Down-regulated",
                                TRUE ~ "Unchanged"))
## Plot
ggplot(data = volcano_dat, aes(x = log2FoldChange,
                               y = -log10(padj),
                               col = Expression)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values = c("red", "grey", "green")) +
  geom_vline(xintercept = c(-1.2, 1.2), col = "black") +
  geom_hline(yintercept = -log10(0.05), col = "black") +
  ggtitle("Volcano Plot")

## MA Plot
## The function plotMA shows the log2 fold changes attributable
## to a given variable over the mean of normalized counts for
## all the samples in the DESeqDataSet.
## Points will be colored red if the adjusted p value is less than 0.1.
## Points which fall out of the window are plotted as open triangles
## pointing either up or down.
# plotMA(res, ylim = c(-2, 2))

##
