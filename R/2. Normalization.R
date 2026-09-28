# Normalization

# load data
data <- read.csv(
  "C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech/synthetic_cfRNA_expression_clinical_data.csv"
)

# view column names
colnames(data)

# Select RNA expression features

rna_features <- c(
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

expression_data <- data[, rna_features]

# Check for missing values
colSums(is.na(expression_data))

# Replace missing expression values: For this synthetic demonstration, missing values are replaced with zero before transformation.
expression_data[is.na(expression_data)] <- 0

# Log2 transformation

# The +1 prevents log2(0) from causing an undefined value.
log_expression <- log2(expression_data + 1)

# Inspect transformed data
head(log_expression)
summary(log_expression)

# Compare raw and transformed expression
par(mfrow = c(1, 2))

hist(
  expression_data$miR451a,
  main = "Raw miR451a",
  xlab = "Expression"
)
hist(
  log_expression$miR451a,
  main = "Log2-transformed miR451a",
  xlab = "Log2(Expression + 1)"
)

par(mfrow = c(1, 1))

# Save transformed dataset
write.csv(
  log_expression,
  "C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech//log2_transformed_expression.csv",
  row.names = FALSE
)
