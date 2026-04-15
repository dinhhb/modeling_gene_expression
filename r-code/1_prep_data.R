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


exprmx <- read_csv("~/PycharmProjects/internship-m2-melanoma/main/dataset/created/gex_mat.csv")
meta <- read_csv("~/PycharmProjects/internship-m2-melanoma/main/dataset/created/sample_source.csv")


meta$batch <- as.factor(meta$source)
meta$pfs <- as.factor(meta$pfs_label)

gene_ids <- exprmx$HGNC
exprmx_mat <- as.matrix(exprmx[, -1])
rownames(exprmx_mat) <- gene_ids

# Keep genes measured in ≥80% of all samples
gene_keep <- rowMeans(!is.na(exprmx_mat)) >= 0.8
count_tbl_low_rm <- exprmx_mat[gene_keep, ]

# Impute remaining NAs with row mean
count_tbl_pca <- t(apply(count_tbl_low_rm, 1, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
}))

cat("Genes retained:", nrow(count_tbl_pca), "\n")
cat("NAs remaining:", sum(is.na(count_tbl_pca)), "\n")
