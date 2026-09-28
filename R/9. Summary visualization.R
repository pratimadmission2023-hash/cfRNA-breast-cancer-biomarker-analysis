# Summary visualization

# Load packages
installed.packages("ggplot2")
library(ggplot2)

# Load original dataset
data <- read.csv(
  "C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech/synthetic_cfRNA_expression_clinical_data.csv",
  check.names = FALSE
)

#  Prepare biomarker data

biomarkers <- data[, c(
  "miR451a",
  "miR486",
  "miR21",
  "miR155",
  "miR205"
)]

# Convert to long format
biomarker_long <- data.frame(
  Sample_ID = rep(data$Sample_ID, times = 5),
  Biomarker = rep(
    c("miR451a", "miR486", "miR21", "miR155", "miR205"),
    each = nrow(data)
  ),
  Expression = c(
    data$miR451a,
    data$miR486,
    data$miR21,
    data$miR155,
    data$miR205
  )
)

# Create summary plot
summary_plot <- ggplot(
  biomarker_long,
  aes(
    x = Biomarker,
    y = Expression
  )
) +
  geom_boxplot(
    na.rm = TRUE
  ) +
  geom_jitter(
    width = 0.15,
    size = 1.8,
    alpha = 0.7,
    na.rm = TRUE
  ) +
  labs(
    title = "Summary of Selected cfRNA Biomarker Expression",
    x = "Biomarker",
    y = "Expression"
  ) +
  theme_minimal()

print(summary_plot)

# Save summary figure
ggsave(
  "C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech/biomarker_summary.png",
  plot = summary_plot,
  width = 9,
  height = 6,
  units = "in",
  dpi = 300
)

# Calculate summary statistics
summary_statistics <- data.frame(
  Biomarker = c(
    "miR451a",
    "miR486",
    "miR21",
    "miR155",
    "miR205"
  ),
  Mean = c(
    mean(data$miR451a, na.rm = TRUE),
    mean(data$miR486, na.rm = TRUE),
    mean(data$miR21, na.rm = TRUE),
    mean(data$miR155, na.rm = TRUE),
    mean(data$miR205, na.rm = TRUE)
  ),
  Median = c(
    median(data$miR451a, na.rm = TRUE),
    median(data$miR486, na.rm = TRUE),
    median(data$miR21, na.rm = TRUE),
    median(data$miR155, na.rm = TRUE),
    median(data$miR205, na.rm = TRUE)
  ),
  SD = c(
    sd(data$miR451a, na.rm = TRUE),
    sd(data$miR486, na.rm = TRUE),
    sd(data$miR21, na.rm = TRUE),
    sd(data$miR155, na.rm = TRUE),
    sd(data$miR205, na.rm = TRUE)
  )
)

# Display results
print(summary_statistics)

# Save summary statistics
write.csv(
  summary_statistics,
  "C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech/biomarker_summary_statistics.csv",
  row.names = FALSE
)

cat("Step 9 summary visualisation completed successfully.\n")
dev.off()
