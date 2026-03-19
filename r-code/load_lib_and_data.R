# Install CRAN packages
install.packages(c("data.table", "ggplot2", "ggprism"), dependencies = TRUE)

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery", "limma", "edgeR", "sva", "lme4"), dependencies = TRUE)

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


# Load data
exprmx <- read_csv("~/PycharmProjects/internship-m2-melanoma/main/dataset/created/gex_mat_rna_seq.csv")
meta <- read_csv("~/PycharmProjects/internship-m2-melanoma/main/dataset/created/sample_source_rna_seq.csv")



