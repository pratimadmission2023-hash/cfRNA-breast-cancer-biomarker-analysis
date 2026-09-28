# Clinical association analysis

# Load data
data <- read.csv(
  "C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech/synthetic_cfRNA_expression_clinical_data.csv",
  check.names = FALSE
)

# Check the data
head(data)

# Create biomarker expression groups
# miR451a median
miR451a_median <- median(
  data$miR451a,
  na.rm = TRUE
)

data$miR451a_group <- ifelse(
  data$miR451a >= miR451a_median,
  "High",
  "Low"
)

data$miR451a_group <- as.factor(data$miR451a_group)


# miR486 median
miR486_median <- median(
  data$miR486,
  na.rm = TRUE
)

data$miR486_group <- ifelse(
  data$miR486 >= miR486_median,
  "High",
  "Low"
)

data$miR486_group <- as.factor(data$miR486_group)

# Check clinical variables

table(data$HER2, useNA = "ifany")

table(data$PR, useNA = "ifany")

table(data$LVSI, useNA = "ifany")

table(data$Lymph_Node_Status, useNA = "ifany")

table(data$Tumour_Grade, useNA = "ifany")

# HER2 association

her2_miR451a <- table(
  data$miR451a_group,
  data$HER2
)

print(her2_miR451a)

her2_miR451a_test <- fisher.test(
  her2_miR451a
)

print(her2_miR451a_test)


her2_miR486 <- table(
  data$miR486_group,
  data$HER2
)

print(her2_miR486)

her2_miR486_test <- fisher.test(
  her2_miR486
)

print(her2_miR486_test)

# PR association

pr_miR451a <- table(
  data$miR451a_group,
  data$PR
)

print(pr_miR451a)

pr_miR451a_test <- fisher.test(
  pr_miR451a
)

print(pr_miR451a_test)


pr_miR486 <- table(
  data$miR486_group,
  data$PR
)

print(pr_miR486)

pr_miR486_test <- fisher.test(
  pr_miR486
)

print(pr_miR486_test)

# LVSI association

lvsi_miR451a <- table(
  data$miR451a_group,
  data$LVSI
)

print(lvsi_miR451a)

lvsi_miR451a_test <- fisher.test(
  lvsi_miR451a
)

print(lvsi_miR451a_test)


lvsi_miR486 <- table(
  data$miR486_group,
  data$LVSI
)

print(lvsi_miR486)

lvsi_miR486_test <- fisher.test(
  lvsi_miR486
)

print(lvsi_miR486_test)

# Lymph-node status association

node_miR451a <- table(
  data$miR451a_group,
  data$Lymph_Node_Status
)

print(node_miR451a)

node_miR451a_test <- fisher.test(
  node_miR451a
)

print(node_miR451a_test)


node_miR486 <- table(
  data$miR486_group,
  data$Lymph_Node_Status
)

print(node_miR486)

node_miR486_test <- fisher.test(
  node_miR486
)

print(node_miR486_test)

# Tumour grade association

grade_miR451a <- table(
  data$miR451a_group,
  data$Tumour_Grade
)

print(grade_miR451a)

grade_miR451a_test <- chisq.test(
  grade_miR451a
)

print(grade_miR451a_test)


grade_miR486 <- table(
  data$miR486_group,
  data$Tumour_Grade
)

print(grade_miR486)

grade_miR486_test <- chisq.test(
  grade_miR486
)

print(grade_miR486_test)

# Save statistical results

clinical_results <- data.frame(
  
  Biomarker = c(
    "miR451a",
    "miR486",
    "miR451a",
    "miR486",
    "miR451a",
    "miR486",
    "miR451a",
    "miR486",
    "miR451a",
    "miR486"
  ),
  
  Clinical_Variable = c(
    "HER2",
    "HER2",
    "PR",
    "PR",
    "LVSI",
    "LVSI",
    "Lymph_Node_Status",
    "Lymph_Node_Status",
    "Tumour_Grade",
    "Tumour_Grade"
  ),
  
  Test = c(
    "Fisher's exact test",
    "Fisher's exact test",
    "Fisher's exact test",
    "Fisher's exact test",
    "Fisher's exact test",
    "Fisher's exact test",
    "Fisher's exact test",
    "Fisher's exact test",
    "Chi-square test",
    "Chi-square test"
  ),
  
  P_value = c(
    her2_miR451a_test$p.value,
    her2_miR486_test$p.value,
    pr_miR451a_test$p.value,
    pr_miR486_test$p.value,
    lvsi_miR451a_test$p.value,
    lvsi_miR486_test$p.value,
    node_miR451a_test$p.value,
    node_miR486_test$p.value,
    grade_miR451a_test$p.value,
    grade_miR486_test$p.value
  )
)

print(clinical_results)

# Save results
write.csv(
  clinical_results,
  "C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech/clinical_association_results.csv",
  row.names = FALSE
)

# Completion message
cat(
  "Clinical association analysis completed successfully.\n"
)