# Correlation analysis of cfRNA biomarkers

# Load packages
library(ggplot2)

# load data
data <- read.csv(
  "C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech/synthetic_cfRNA_expression_clinical_data.csv",
  check.names = FALSE
)

# Check data
head(data)

# Select biomarkers
correlation_data <- data[, c(
  "miR451a",
  "miR486"
)]

# Remove samples with missing values
correlation_data <- na.omit(correlation_data)

# Check number of samples
nrow(correlation_data)

#  Calculate Pearson correlation

correlation_test <- cor.test(
  correlation_data$miR451a,
  correlation_data$miR486,
  method = "pearson"
)

print(correlation_test)

# Extract correlation coefficient
correlation_coefficient <- correlation_test$estimate
p_value <- correlation_test$p.value

cat("Pearson correlation coefficient:", correlation_coefficient, "\n")
cat("P-value:", p_value, "\n")

# Create correlation plot
correlation_plot <- ggplot(
  correlation_data,
  aes(
    x = miR451a,
    y = miR486
  )
) +
  geom_point(
    size = 3
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  labs(
    title = "Correlation between miR451a and miR486",
    x = "miR451a expression",
    y = "miR486 expression"
  ) +
  theme_minimal()

print(correlation_plot)

# Save correlation plot
ggsave(
  "C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech/miR451a_miR486_correlation.png",
  plot = correlation_plot,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)

# Save statistical results
correlation_results <- data.frame(
  Biomarker_1 = "miR451a",
  Biomarker_2 = "miR486",
  Method = "Pearson correlation",
  Correlation = as.numeric(correlation_coefficient),
  P_value = p_value
)

write.csv(
  correlation_results,
  "C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech/correlation_results.csv",
  row.names = FALSE
)

cat("Correlation analysis completed successfully.\n")

read.csv("C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech/correlation_results.csv")
