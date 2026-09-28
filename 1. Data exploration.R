# Data exploration

# load data
data <- read.csv(
  "C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech/synthetic_cfRNA_expression_clinical_data.csv"
)
# Check dimensions
dim(data)

# View first rows
head(data)

# Check structure
str(data)

# Check missing values
colSums(is.na(data))

# Summary statistics
summary(data)