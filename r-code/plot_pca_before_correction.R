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

# Log transform=
shift <- abs(min(count_tbl_pca)) + 1
count_tbl_log <- log2(count_tbl_pca + shift)

# Transpose count matrix for PCA
count_tbl_low_rm_t <- as.data.frame(t(count_tbl_log))

# Perform PCA
pca_prep <- prcomp(count_tbl_low_rm_t, scale. = TRUE)

# Create PCA plot colored by batch
pca_plot_batch <- autoplot(pca_prep, label = F, data = meta,
                           colour = "batch") +
  theme_prism(base_size = 16) +
  ggtitle("PCA Plot Before Batch Effect Correction") +
  theme(plot.title = element_text(hjust = 0.5))

# Save the plot
ggsave("plots/pca_plot_batch_unadjusted.png", pca_plot_batch,
       device = "png", units = "cm", height = 12, width = 22)
