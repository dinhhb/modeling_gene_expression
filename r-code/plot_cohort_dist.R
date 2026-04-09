library(ggplot2)
library(tidyr)

# 1. Prepare data for plotting
# We convert the matrix to a "long" format that ggplot likes
plot_data <- data.frame(t(exprmx_mat))
plot_data$Source <- meta$source

# Reshape to long format
plot_long <- pivot_longer(plot_data, 
                          cols = -Source, 
                          names_to = "Gene", 
                          values_to = "Expression")

# 2. Density Plot (Shows the 'Shape' and 'Skew')
# We use log10 on the x-axis to handle the 218k vs -16 range
ggplot(plot_long, aes(x = Expression, fill = Source)) +
  geom_density(alpha = 0.5) +
  scale_x_log10() + 
  theme_prism() +
  labs(title = "Gene Expression Density by Cohort",
       x = "Expression (Log10 Scale)", 
       y = "Density") +
  theme(legend.position = "bottom")

# 3. Boxplot (Shows the Medians and Outliers)
ggplot(plot_long, aes(x = Source, y = Expression, fill = Source)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_y_log10() +
  theme_prism() +
  labs(title = "Range Comparison Across Cohorts",
       y = "Expression (Log10 Scale)") +
  coord_flip() # Flip for easier reading of names