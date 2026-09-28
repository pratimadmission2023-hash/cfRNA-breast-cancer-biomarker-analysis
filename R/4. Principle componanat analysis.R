# Principal Component Analysis of cfRNA biomarker expression

# Load transformed expression data
expression_data <- read.csv(
  "C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech//log2_transformed_expression.csv",
  check.names = FALSE
)

# Biomarkers to include in PCA
biomarkers <- c(
  "miR451a",
  "miR486",
  "miR21",
  "miR155",
  "miR205",
  "SCGB2A2",
  "AGR2",
  "TFF1",
  "MALAT1",
  "RPLP0"
)

# Extract expression data
pca_data <- expression_data[, biomarkers]

# Convert to matrix
pca_matrix <- as.matrix(pca_data)

# Add sample IDs as row names
rownames(pca_matrix) <- expression_data$Sample_ID

# Run PCA
pca_result <- prcomp(
  pca_matrix,
  scale. = TRUE
)

# View PCA summary
summary(pca_result)

# Calculate percentage variance explained
variance_explained <- pca_result$sdev^2 /
  sum(pca_result$sdev^2) * 100

# Create PCA plot
plot(
  pca_result$x[, 1],
  pca_result$x[, 2],
  xlab = paste0(
    "PC1 (",
    round(variance_explained[1], 1),
    "% variance)"
  ),
  ylab = paste0(
    "PC2 (",
    round(variance_explained[2], 1),
    "% variance)"
  ),
  main = "PCA of cfRNA Biomarker Expression",
  pch = 19
)

# Add sample labels
text(
  pca_result$x[, 1],
  pca_result$x[, 2],
  labels = expression_data$Sample_ID,
  pos = 3,
  cex = 0.7
)

# Save PCA plot
png(
  "C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech/PCA.png",
  width = 2400,
  height = 2000,
  res = 300
)

plot(
  pca_result$x[, 1],
  pca_result$x[, 2],
  xlab = paste0(
    "PC1 (",
    round(variance_explained[1], 1),
    "% variance)"
  ),
  ylab = paste0(
    "PC2 (",
    round(variance_explained[2], 1),
    "% variance)"
  ),
  main = "PCA of cfRNA Biomarker Expression",
  pch = 19
)

text(
  pca_result$x[, 1],
  pca_result$x[, 2],
  labels = expression_data$Sample_ID,
  pos = 3,
  cex = 0.7
)

install.packages("ggplot2")
library(ggplot2)

# Create PCA plot using ggplot2


# Load transformed expression data
expression_data <- read.csv(
  "C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech/log2_transformed_expression.csv",
  check.names = FALSE
)

# Load original dataset to obtain Sample IDs
original_data <- read.csv(
  "C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech/synthetic_cfRNA_expression_clinical_data.csv",
  check.names = FALSE
)

# Biomarkers to include in PCA
biomarkers <- c(
  "miR451a",
  "miR486",
  "miR21",
  "miR155",
  "miR205",
  "SCGB2A2",
  "AGR2",
  "TFF1",
  "MALAT1",
  "RPLP0"
)

# Extract expression data
pca_data <- expression_data[, biomarkers]

# Convert to matrix
pca_matrix <- as.matrix(pca_data)


# Run PCA
pca_result <- prcomp(
  pca_matrix,
  scale. = TRUE
)

# PCA summary
summary(pca_result)

# Calculate percentage variance explained
variance_explained <- (
  pca_result$sdev^2 /
    sum(pca_result$sdev^2)
) * 100

# Create PCA data frame
pca_plot_data <- data.frame(
  Sample_ID = original_data$Sample_ID,
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2]
)

# Check number of samples
print(nrow(pca_plot_data))

# Create PCA plot
pca_plot <- ggplot(
  pca_plot_data,
  aes(x = PC1, y = PC2)
) +
  geom_point(size = 3) +
  geom_text(
    aes(label = Sample_ID),
    vjust = -0.7,
    size = 3
  ) +
  labs(
    title = "PCA of cfRNA Biomarker Expression",
    x = paste0(
      "PC1 (",
      round(variance_explained[1], 1),
      "% variance)"
    ),
    y = paste0(
      "PC2 (",
      round(variance_explained[2], 1),
      "% variance)"
    )
  ) +
  theme_minimal()

# Display PCA in RStudio
print(pca_plot)

# Save PCA
ggsave(
  filename = "C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech/PCA2.png",
  plot = pca_plot,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)

cat("PCA plot saved successfully.\n")

file.info("C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech/PCA2.png")$size
