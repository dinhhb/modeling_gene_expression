# Differential expression across your 3 PFS groups
design <- model.matrix(~0 + pfs, data = meta)
fit <- lmFit(batch_corrected_limma, design)

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
sig_genes <- rownames(top_genes[top_genes$adj.P.Val < 0.05, ])
corrected_sig <- batch_corrected_limma[sig_genes, ]

# Convert to dataframe with gene names as first column
corrected_sig_df <- as.data.frame(corrected_sig)
corrected_sig_df <- cbind(HGNC = rownames(corrected_sig), corrected_sig_df)

# Save to CSV
write_csv(corrected_sig_df, 
          "~/PycharmProjects/internship-m2-melanoma/main/dataset/created/gex_mat_corrected_sig.csv")