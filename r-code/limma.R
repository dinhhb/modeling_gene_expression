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
          "~/PycharmProjects/internship-m2-melanoma/main/dataset/gene_expressions_matrix_corrected.csv")