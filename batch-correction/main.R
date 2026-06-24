# In RStudio, open Project at path ~/PycharmProjects/internship-m2-melanoma/main/batch-correction

library(data.table)
library(ggplot2)
library(ggprism)
library(limma)
library(edgeR)
library(readr)
library(ggfortify)
library(pheatmap)
library(RColorBrewer)

# 1. FUNCTION: PLOT PCA
plot_pca_func <- function(gex_mat, meta, title, filename) {
  # Filter out 0-variance genes
  gene_var <- apply(gex_mat, 1, var)
  gex_mat  <- gex_mat[gene_var > 0, ]

  pca_res <- prcomp(t(gex_mat), scale. = TRUE)
  
  p <- autoplot(pca_res, label = FALSE, data = meta, colour = "batch") +
    theme_prism(base_size = 16) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(filename, p, device = "png", units = "cm", height = 15, width = 20)
  return(p)
}

# 2. FUNCTION: PLOT HEATMAP
plot_heatmap_func <- function(gex_mat, meta, title, filename) {
  annot_df <- data.frame(
    PFS = meta$pfs,
    Batch = meta$batch,
    row.names = colnames(gex_mat)
  )
  
  n_batches <- nlevels(meta$batch)
  n_pfs     <- nlevels(meta$pfs)
  
  batch_cols <- setNames(RColorBrewer::brewer.pal(max(3, n_batches), "Set2")[seq_len(n_batches)],
                         levels(meta$batch))
  pfs_cols <- setNames(RColorBrewer::brewer.pal(max(3, n_pfs), "Dark2")[seq_len(n_pfs)],
                       levels(meta$pfs))
  
  annot_cols <- list(Batch = batch_cols, PFS = pfs_cols)
  
  # Clustering/Top Var Genes
  top25_idx <- order(apply(gex_mat, 1, var), decreasing = TRUE)[1:20]
  gex_mat_top25 <- gex_mat[top25_idx, ]
  
  pheatmap(gex_mat_top25,
           annotation_col = annot_df,
           annotation_colors = annot_cols,
           scale = "row",
           color = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
           show_rownames = TRUE,
           show_colnames = FALSE,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "ward.D2",
           main = title,
           fontsize = 11,
           filename = filename,
           width = 22 / 2.54,
           height = 20 / 2.54)
}
# --- Prepare data ---

gex_mat_raw <- read_csv("~/PycharmProjects/internship-m2-melanoma/main/dataset/created/binary/gex_mat.csv")
meta <- read_csv("~/PycharmProjects/internship-m2-melanoma/main/dataset/created/binary/sample_source.csv")

meta$batch <- as.factor(meta$source)
meta$pfs <- as.factor(meta$pfs_label)

sample_ids <- gex_mat_raw[[1]]
gex_mat    <- as.matrix(gex_mat_raw[, -1])
rownames(gex_mat) <- sample_ids 

# Keep genes measured in ≥80% of all samples
gene_keep <- rowMeans(!is.na(gex_mat)) >= 0.8
gex_mat_filtered <- gex_mat[gene_keep, ]

# Impute remaining NAs with row mean
gex_mat_imputed <- t(apply(gex_mat_filtered, 1, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
}))

cat("Genes retained:", nrow(gex_mat_imputed), "\n")
cat("NAs remaining:", sum(is.na(gex_mat_imputed)), "\n")

# Log transform
shift <- abs(min(gex_mat_imputed)) + 1
gex_mat_log <- log2(gex_mat_imputed + shift)


# --- Execution ---

# Visualize Before
plot_pca_func(gex_mat_log, meta, "Before Correction (PCA)", "plots/pca_before.png")
plot_heatmap_func(gex_mat_log, meta, "Before Correction", "plots/heatmap_var_before.png")

# Batch Correction
design_mat <- model.matrix(~pfs, data = meta)
# corrected_mat <- removeBatchEffect(gex_mat_log, design = design_mat, batch = meta$batch)
corrected_mat <- removeBatchEffect(gex_mat_log, batch = meta$batch)

gex_mat_df <- as.data.frame(corrected_mat)
gex_mat_df <- cbind(sample_id = rownames(gex_mat_log), gex_mat_df)

#write_csv(gex_mat_df, 
#          "~/PycharmProjects/internship-m2-melanoma/main/dataset/created/binary/gex_mat_corrected.csv")

# Visualize After
plot_pca_func(corrected_mat, meta, "After limma Correction (PCA)", "plots/pca_after.png")
plot_heatmap_func(corrected_mat, meta, "After Correction", "plots/heatmap_var_after.png")

