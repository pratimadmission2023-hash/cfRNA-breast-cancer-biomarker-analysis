# Heatmap

# install heatmap
install.packages("pheatmap")

# load heatmap package
library(pheatmap)

# Load the log2-transformed expression data
expression_data <- read.csv(
  "C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech//log2_transformed_expression.csv",
  row.names = NULL
)

# Select biomarker features
biomarkers <- c(
  "miR451a",
  "miR486",
  "miR21",
  "miR155",
  "miR205",
  "SCGB2A2",
  "AGR2",
  "TFF1"
)

heatmap_data <- expression_data[, biomarkers]

#  Convert data into a matrix
heatmap_matrix <- as.matrix(heatmap_data)

#Give samples meaningful row names
original_data <- read.csv(
  "C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech//log2_transformed_expression.csv"
)
rownames(heatmap_matrix) <- original_data$Sample_ID

# generate heat map
pheatmap(
  heatmap_matrix,
  scale = "column",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Synthetic cfRNA Biomarker Expression"
)

# save heat map
png(
  "C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech//heatmap.jpeg",  
  width = 1600,
  height = 1400,
  res = 200
)

pheatmap(
  heatmap_matrix,
  scale = "column",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Synthetic cfRNA Biomarker Expression"
)

# checking existence of file
file.exists("C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech//heatmap.jpeg")

# size of heatmap
file.info("C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech//heatmap.jpeg")$size

# opening directly from R
browseURL(normalizePath("C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech//heatmap.jpeg"))

#Test whether R can create a PNG at all
png("test.png", width = 1000, height = 800, res = 150)

plot(
  1:10,
  1:10,
  main = "PNG Test",
  xlab = "X",
  ylab = "Y"
)

dev.off()

# size of test
file.info("test.png")

# testing heatmap
library(pheatmap)

test_matrix <- matrix(
  rnorm(100),
  nrow = 10,
  ncol = 10
)

png(
  "test_heatmap.png",
  width = 2000,
  height = 1600,
  res = 300
)

pheatmap(
  test_matrix,
  main = "Test Heatmap"
)

dev.off()
file.info("test_heatmap.png")$size

# closing other graphics device
while (!is.null(dev.list())) {
  dev.off()
}

# test
library(pheatmap)

test_matrix <- matrix(
  rnorm(100),
  nrow = 10,
  ncol = 10
)

pheatmap(
  test_matrix,
  main = "Test Heatmap",
  filename = "test_heatmap.png",
  width = 8,
  height = 6
)
file.info("test_heatmap.png")$size

# fix actual heat map
library(pheatmap)

# Load transformed expression data
expression_data <- read.csv(
  "C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech//log2_transformed_expression.csv",
  check.names = FALSE
)

# Select biomarkers
biomarkers <- c(
  "miR451a",
  "miR486",
  "miR21",
  "miR155",
  "miR205",
  "SCGB2A2",
  "AGR2",
  "TFF1"
)

# Select expression data
heatmap_data <- expression_data[, biomarkers]

# Convert to matrix
heatmap_matrix <- as.matrix(heatmap_data)

# Add sample IDs
rownames(heatmap_matrix) <- expression_data$Sample_ID

biomarkers

# column names
colnames(expression_data)

# Create and save heatmap
pheatmap(
  heatmap_matrix,
  scale = "column",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "cfRNA Biomarker Expression Heatmap",
  filename = "C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech//heatmap2.png",
  width = 10,
  height = 8
)

file.info("C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech//heatmap2.png")$size