# Log transform
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

# ─────────────────────────────────────────────
# Heatmap: Gene × Sample (top 25 variable genes)
# ─────────────────────────────────────────────
library(pheatmap)

batch_vec <- meta$batch          # factor: as.factor(meta$source)

# Annotation for heatmap columns (samples)
annot_df <- data.frame(Batch = batch_vec,
                       row.names = colnames(count_tbl_log))

# Colour palette per batch level
n_batches  <- nlevels(batch_vec)
batch_cols <- setNames(RColorBrewer::brewer.pal(max(3, n_batches), "Set2")[seq_len(n_batches)],
                       levels(batch_vec))
annot_cols <- list(Batch = batch_cols)

# ── Before correction — Gene × Sample heatmap (top 25 variable genes) ──────
top25_idx  <- order(apply(count_tbl_log, 1, var), decreasing = TRUE)[1:25]
expr_top25 <- count_tbl_log[top25_idx, ]

p_before <- pheatmap(expr_top25,
                     annotation_col    = annot_df,        # batch label on columns (samples)
                     annotation_colors = annot_cols,
                     scale             = "row",           # z-score per gene — essential for gene×sample
                     color             = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
                     show_rownames     = TRUE,           # 25 gene names won't fit
                     show_colnames     = FALSE,
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean",
                     clustering_method        = "ward.D2",
                     main              = "Before Batch Correction — (top 25 variable genes)",
                     fontsize          = 11,
                     filename          = "plots/heatmap_before_correction.png",
                     width             = 22 / 2.54,      # cm → inches
                     height            = 20 / 2.54)