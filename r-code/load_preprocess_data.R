# Install CRAN packages
#install.packages(c("data.table", "ggplot2", "ggprism"), dependencies = TRUE)

# Install Bioconductor packages
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install(c("GEOquery", "limma", "edgeR", "sva", "lme4"), dependencies = TRUE)

# Load the libraries
library(data.table)
library(ggplot2)
library(ggprism)
library(limma)
library(edgeR)
library(sva)
library(lme4)
library(readr)
library(ggfortify)
library(ggprism)
library(pheatmap)


# Load data
exprmx <- read_csv("~/PycharmProjects/internship-m2-melanoma/main/dataset/created/gex_mat/gex_pre.csv")
meta <- read_csv("~/PycharmProjects/internship-m2-melanoma/main/dataset/created/sample_source/gex_pre.csv")


meta$batch <- as.factor(meta$source)
# meta$treatment ...

gene_ids <- exprmx$HGNC
exprmx_mat <- as.matrix(exprmx[, -1])
rownames(exprmx_mat) <- gene_ids

# Count samples where gene is expressed (> 0), ignoring NAs
n_expressed <- rowSums(exprmx_mat > 0, na.rm = TRUE)

# Count samples where gene was actually measured (not NA)
n_measured  <- rowSums(!is.na(exprmx_mat))

# Filter low-expressed genes
perc_keep <- 0.8    # keep genes expressed at least 80% of samples

# Option A — expressed in ≥80% of MEASURED samples (recommended)
# respects the fact that NAs = not measured, not absent
gene_keep <- (n_expressed / n_measured) >= perc_keep

# Create a filtered count matrix
count_tbl_low_rm <- exprmx_mat[gene_keep, ]
dim(count_tbl_low_rm)

# Filter genes with >25% NA
gene_na_rate <- rowMeans(is.na(count_tbl_low_rm))
count_tbl_low_rm2 <- count_tbl_low_rm[gene_na_rate <= 0.25, ]
cat("Genes remaining:", nrow(count_tbl_low_rm2), "\n")  # expect ~15,570

# Row mean imputation for remaining NAs
count_tbl_pca <- t(apply(count_tbl_low_rm2, 1, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
}))

cat("NAs after imputation:", sum(is.na(count_tbl_pca)), "\n")  # should be 0
