# Plot NA rate distribution across genes
library(ggplot2)

gene_na_rate <- rowMeans(is.na(count_tbl_low_rm))

na_df <- data.frame(
  bucket = c("0–10%", "10–25%", "25–50%", "50–75%", "75–100%"),
  count  = as.numeric(table(cut(gene_na_rate,
                                breaks = c(0, 0.1, 0.25, 0.5, 0.75, 1),
                                include.lowest = TRUE))),
  status = c("keep", "keep", "borderline", "drop", "drop")
)

na_df$bucket <- factor(na_df$bucket,
                       levels = c("0–10%", "10–25%", "25–50%", "50–75%", "75–100%"))

ggplot(na_df, aes(x = bucket, y = count, fill = status)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = scales::comma(count)),
            vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("keep"       = "#1D9E75",
                               "borderline" = "#EF9F27",
                               "drop"       = "#E24B4A"),
                    labels = c("keep"       = "Keep",
                               "borderline" = "Borderline",
                               "drop"       = "Drop")) +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Gene missingness distribution",
       x     = "NA rate per gene",
       y     = "Number of genes",
       fill  = NULL) +
  theme_bw(base_size = 14) +
  theme(plot.title   = element_text(hjust = 0.5),
        legend.position = "top")
