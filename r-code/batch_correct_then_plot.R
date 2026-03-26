batch_corrected_limma <- removeBatchEffect(count_tbl_log,
                                         batch = meta$batch)

# Visualize the results with PCA
pca_prep_batch_limma <- prcomp(t(batch_corrected_limma), scale. = TRUE)

limma_pca_plot <- autoplot(pca_prep_batch_limma,
                           label  = F,
                           data   = meta,
                           colour = "batch") +
  theme_prism(base_size = 16) +
  ggtitle("PCA Plot After limma Batch Correction") +
  theme(plot.title = element_text(hjust = 0.5))

# Display and save the plot
print(limma_pca_plot)
ggsave("plots/pca_plot_batch_limma_corrected.png", limma_pca_plot,
       device = "png", units = "cm", height = 12, width = 22)

# Convert to dataframe with gene names as first column
batch_corrected_df <- as.data.frame(batch_corrected_limma)
batch_corrected_df <- cbind(HGNC = rownames(batch_corrected_limma), batch_corrected_df)

# Save to csv
write_csv(batch_corrected_df, 
          "~/PycharmProjects/internship-m2-melanoma/main/dataset/created/gex_mat_corrected/gex_pre.csv")


# ── After limma correction — Gene × Sample heatmap (top 25 variable genes) ──
top25_idx_corrected  <- order(apply(batch_corrected_limma, 1, var), decreasing = TRUE)[1:25]
expr_top25_corrected <- batch_corrected_limma[top25_idx_corrected, ]

p_after <- pheatmap(expr_top25_corrected,
                    annotation_col    = annot_df,
                    annotation_colors = annot_cols,
                    scale             = "row",            # z-score per gene
                    color             = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
                    show_rownames     = TRUE,
                    show_colnames     = FALSE,
                    clustering_distance_rows = "euclidean",
                    clustering_distance_cols = "euclidean",
                    clustering_method        = "ward.D2",
                    main              = "After limma Batch Correction (top 25 variable genes)",
                    fontsize          = 11,
                    filename          = "plots/heatmap_after_correction.png",
                    width             = 22 / 2.54,
                    height            = 20 / 2.54)