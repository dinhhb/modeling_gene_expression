# --- 1. Load Libraries ---
library(data.table)
library(ggplot2)
library(ggprism)
library(limma)
library(edgeR)
library(readr)
library(ggfortify)
library(pheatmap)

# --- 2. Load and Initial Prep ---
exprmx <- read_csv("~/PycharmProjects/internship-m2-melanoma/main/dataset/created/gex_mat.csv")
meta <- read_csv("~/PycharmProjects/internship-m2-melanoma/main/dataset/created/sample_source.csv")

# Ensure metadata factors are set
meta$batch <- as.factor(meta$source)
meta$pfs   <- as.factor(meta$pfs_label)

# Convert expression to matrix and set gene names
gene_ids <- exprmx$HGNC
exprmx_mat <- as.matrix(exprmx[, -1])
rownames(exprmx_mat) <- gene_ids

# --- 3. Filtering & Imputation ---
# Filter: Keep genes measured (not NA) in at least 80% of samples
n_measured <- rowSums(!is.na(exprmx_mat))
gene_keep <- (n_measured / ncol(exprmx_mat)) >= 0.8
count_tbl_filtered <- exprmx_mat[gene_keep, ]

# Impute NAs using row means (Must be done before scaling)
count_tbl_imputed <- t(apply(count_tbl_filtered, 1, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
}))

# --- 4. Transformation & Batch-wise Alignment ---
# Log2 transform: handle the -16 values by capping at 0 (pmax)
count_tbl_log <- log2(pmax(count_tbl_imputed, 0) + 1)

# CRITICAL: Scale within each cohort to align the boxes from your plot
# This prevents the "Yan" cohort (218k) from dominating the "Kwong" cohort
count_tbl_aligned <- do.call(cbind, lapply(unique(meta$source), function(s) {
  cols <- which(meta$source == s)
  sub_mat <- count_tbl_log[, cols]
  
  # Safety check: if a gene has 0 variance in a batch, return 0s to avoid NaN
  scaled_sub <- t(apply(sub_mat, 1, function(row) {
    if(sd(row) == 0) return(rep(0, length(row))) 
    return(scale(row))
  }))
  colnames(scaled_sub) <- colnames(sub_mat)
  return(scaled_sub)
}))

# Remove genes that became constant across ALL samples (Prevents SVD error)
count_tbl_aligned <- count_tbl_aligned[apply(count_tbl_aligned, 1, var) > 0, ]

# --- 5. PCA: Before Correction (Aligned Data) ---
# We use scale. = FALSE because we already Z-scored manually
pca_before <- prcomp(t(count_tbl_aligned), scale. = FALSE)
pca_plot_before <- autoplot(pca_before, data = meta, colour = "pfs") +
  theme_prism(base_size = 14) +
  ggtitle("PCA: Aligned Cohorts (Pre-Limma)") +
  theme(plot.title = element_text(hjust = 0.5))
print(pca_plot_before)

# 1. Identify the top 25 most variable genes in the aligned data
# We use variance to find genes that actually "move" across the 125 samples
top25_idx_pre <- order(apply(count_tbl_aligned, 1, var), decreasing = TRUE)[1:25]
expr_top25_pre <- count_tbl_aligned[top25_idx_pre, ]

# 2. Setup Annotation Dataframe
# Using the meta_final (or meta) that matches the columns of count_tbl_aligned
annot_df <- data.frame(
  Batch = meta$batch, 
  PFS = meta$pfs, 
  row.names = colnames(count_tbl_aligned)
)

# 3. Create the Heatmap
pheatmap(expr_top25_pre,
         annotation_col = annot_df,
         scale = "none",          # Data is already Z-scored per batch
         clustering_method = "ward.D2",
         clustering_distance_cols = "euclidean",
         show_colnames = FALSE,
         show_rownames = TRUE,
         main = "Heatmap: Aligned Cohorts (Pre-Limma Correction)",
         color = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
         fontsize_row = 10)


# --- 6. Limma Batch Correction ---
# Create the protection design for your target class (PFS)
# Ensure your meta$pfs has no NAs here!
design_mat <- model.matrix(~pfs, data = meta)

batch_corrected <- removeBatchEffect(
  count_tbl_aligned,
  batch = meta$batch,
  design = design_mat
)

# --- 7. PCA: After Correction ---
pca_after <- prcomp(t(batch_corrected), scale. = FALSE)
pca_plot_after <- autoplot(pca_after, data = meta, colour = "pfs") +
  theme_prism(base_size = 14) +
  ggtitle("PCA: After Limma Correction") +
  theme(plot.title = element_text(hjust = 0.5))
print(pca_plot_after)

# --- 8. Heatmap: Top 25 Variable Genes ---
top25_idx <- order(apply(batch_corrected, 1, var), decreasing = TRUE)[1:25]
expr_top25 <- batch_corrected[top25_idx, ]

# Setup heatmap annotation
annot_df <- data.frame(Batch = meta$batch, PFS = meta$pfs, row.names = colnames(batch_corrected))

pheatmap(expr_top25,
         annotation_col = annot_df,
         scale = "none", # Data is already Z-scored
         clustering_method = "ward.D2",
         color = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
         main = "Final Corrected Heatmap (Top 25 Genes)",
         show_colnames = FALSE)

# --- 9. Export Corrected Data ---
final_df <- as.data.frame(batch_corrected)
final_df <- cbind(HGNC = rownames(batch_corrected), final_df)
write_csv(final_df, "~/PycharmProjects/internship-m2-melanoma/main/dataset/created/gex_mat_corrected.csv")