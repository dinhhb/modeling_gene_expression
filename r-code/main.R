library(data.table)
library(ggplot2)
library(ggprism)
library(limma)
library(edgeR)
library(readr)
library(ggfortify)
library(pheatmap)

#--------------------------1.PREPARE DATA--------------------------
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

# Log transform
shift <- abs(min(count_tbl_pca)) + 1
count_tbl_log <- log2(count_tbl_pca + shift)

#--------------------------2.Genes selection using lmFit()--------------------------
# Clean batch level names before building design matrix
meta$batch_clean <- make.names(meta$batch)
meta$batch_clean <- as.factor(meta$batch_clean)

# Verify
levels(meta$batch_clean)  # should show valid R names e.g. "Hugo.et.al."

# Differential expression across your 3 PFS groups
design <- model.matrix(~0 + pfs + batch_clean, data = meta)
fit <- lmFit(count_tbl_log, design)

# Define contrasts between groups
contrast_mat <- makeContrasts(
  pfs0 - pfs2,
  pfs0 - pfs1,
  pfs2 - pfs1,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_mat)
fit2 <- eBayes(fit2)

# F-test across all contrasts — picks genes that differ in any comparison
top_genes <- topTable(fit2, number = Inf, sort.by = "F")

# Subset corrected matrix to only significant genes
sig_genes <- rownames(top_genes[top_genes$adj.P.Val < 0.2, ])

#--------------------------3.1.PLOTS BEFORE CORRECTION (PCA)--------------------------
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
ggsave("plots/pca_before.png", pca_plot_batch,
       device = "png", units = "cm", height = 12, width = 22)

#--------------------------3.2.PLOTS BEFORE CORRECTION (Heatmap)--------------------------
# Annotation for heatmap columns (samples)
annot_df <- data.frame(
  PFS = meta$pfs,
  Batch = meta$batch,
  row.names = colnames(count_tbl_log)
)

# Colour palette per batch level
n_batches  <- nlevels(meta$batch)
batch_cols <- setNames(RColorBrewer::brewer.pal(max(3, n_batches), "Set2")[seq_len(n_batches)],
                       levels(meta$batch))
n_pfs_levels <- nlevels(meta$pfs)
pfs_cols <- setNames(RColorBrewer::brewer.pal(max(3, n_pfs_levels), "Dark2")[seq_len(n_pfs_levels)],
                     levels(meta$pfs))

annot_cols <- list(
  Batch = batch_cols,
  PFS = pfs_cols
)

# top 25 variable genes
top25_idx  <- order(apply(count_tbl_log, 1, var), decreasing = TRUE)[1:25]
expr_top25 <- count_tbl_log[top25_idx, ]

p_before_var <- pheatmap(expr_top25,
                     annotation_col    = annot_df,        # batch label on columns (samples)
                     annotation_colors = annot_cols,
                     scale             = "row",           # z-score per gene — essential for gene×sample
                     color             = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
                     show_rownames     = TRUE,           # 25 gene names won't fit
                     show_colnames     = FALSE,
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean",
                     clustering_method        = "ward.D2",
                     main              = "Before Batch Correction (top 25 most variable genes)",
                     fontsize          = 11,
                     filename          = "plots/heatmap_var_before.png",
                     width             = 22 / 2.54,      # cm → inches
                     height            = 20 / 2.54)

# top 25 genes by F-statistic 
top25_f_genes    <- rownames(top_genes)[1:25]          # top_genes already sorted by F
expr_top25_f     <- count_tbl_log[top25_f_genes, ]

p_before_f <- pheatmap(expr_top25_f,
                     annotation_col    = annot_df,        # batch label on columns (samples)
                     annotation_colors = annot_cols,
                     scale             = "row",           # z-score per gene — essential for gene×sample
                     color             = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
                     show_rownames     = TRUE,           # 25 gene names won't fit
                     show_colnames     = FALSE,
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean",
                     clustering_method        = "ward.D2",
                     main              = "Before Batch Correction (top 25 limma F-test p-value)",
                     fontsize          = 11,
                     filename          = "plots/heatmap_f_before.png",
                     width             = 22 / 2.54,      # cm → inches
                     height            = 20 / 2.54)


#--------------------------4.BATCH CORRECTION (limma)--------------------------
design_mat <- model.matrix(~pfs, data = meta)

batch_corrected_limma <- removeBatchEffect(
  count_tbl_log,
  design = design_mat,
  batch = meta$batch
)

corrected_sig <- batch_corrected_limma[sig_genes, ]

# Convert to dataframe with gene names as first column
corrected_sig_df <- as.data.frame(corrected_sig)
corrected_sig_df <- cbind(HGNC = rownames(corrected_sig), corrected_sig_df)

# Save to CSV
write_csv(corrected_sig_df, 
          "~/PycharmProjects/internship-m2-melanoma/main/dataset/created/gex_mat_corrected_sig.csv")


#--------------------------5.1.PLOTS AFTER BATCH CORRECTION (PCA)--------------------------
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
#write_csv(batch_corrected_df, 
#          "~/PycharmProjects/internship-m2-melanoma/main/dataset/created/gex_mat_corrected.csv")

#--------------------------5.1.PLOTS AFTER BATCH CORRECTION (HEATMAP)--------------------------
# top 25 variable gene
top25_idx_corrected  <- order(apply(batch_corrected_limma, 1, var), decreasing = TRUE)[1:25]
expr_top25_corrected <- batch_corrected_limma[top25_idx_corrected, ]

p_before_var <- pheatmap(expr_top25_corrected,
                     annotation_col    = annot_df,        # batch label on columns (samples)
                     annotation_colors = annot_cols,
                     scale             = "row",           # z-score per gene — essential for gene×sample
                     color             = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
                     show_rownames     = TRUE,           # 25 gene names won't fit
                     show_colnames     = FALSE,
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean",
                     clustering_method        = "ward.D2",
                     main              = "After Batch Correction (top 25 most variable genes)",
                     fontsize          = 11,
                     filename          = "plots/heatmap_var_after.png",
                     width             = 22 / 2.54,      # cm → inches
                     height            = 20 / 2.54)

# top 25 genes by F-statistic 
expr_top25_f_corrected     <- batch_corrected_limma[top25_f_genes, ]

p_after_f <- pheatmap(expr_top25_f_corrected,
                    annotation_col    = annot_df,
                    annotation_colors = annot_cols,
                    scale             = "row",            # z-score per gene
                    color             = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
                    show_rownames     = TRUE,
                    show_colnames     = FALSE,
                    clustering_distance_rows = "euclidean",
                    clustering_distance_cols = "euclidean",
                    clustering_method        = "ward.D2",
                    main              = "After limma Batch Correction (top 25 limma F-test p-value)",
                    fontsize          = 11,
                    filename          = "plots/heatmap_f_after.png",
                    width             = 22 / 2.54,
                    height            = 20 / 2.54)



