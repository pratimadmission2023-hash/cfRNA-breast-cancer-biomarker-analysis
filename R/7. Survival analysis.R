# survival analysis

# Install survival packages if needed
install.packages("survival")
install.packages("survminer")

library(survival)
library(survminer)

# load data
data <- read.csv(
  "C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech/synthetic_cfRNA_expression_clinical_data.csv",
  check.names = FALSE
)

head(data)

# Check survival variables

summary(data$Survival_Months)
table(data$Survival_Status)

# Create survival object

survival_object <- Surv(
  time = data$Survival_Months,
  event = data$Survival_Status
)

# column data
class(data$Survival_Status)
table(data$Survival_Status, useNA = "ifany")

str(data$Survival_Status)

data$Survival_Status <- ifelse(
  data$Survival_Status == "Deceased",
  1,
  0
)

table(data$Survival_Status)

survival_object <- Surv(
  time = data$Survival_Months,
  event = data$Survival_Status
)

# Convert survival status to 0/1
# Censored = 0
# Event = 1
data$Survival_Status <- ifelse(
  data$Survival_Status == "Event",
  1,
  0
)

# Check conversion
table(data$Survival_Status, useNA = "ifany")

data <- read.csv(
  "C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech/synthetic_cfRNA_expression_clinical_data.csv",
  check.names = FALSE
)

# Check the original survival status
table(data$Survival_Status, useNA = "ifany")

# Censored = 0
# Event = 1

data$Survival_Status <- ifelse(
  data$Survival_Status == "Event",
  1,
  0
)

# Check conversion
table(data$Survival_Status, useNA = "ifany")

# survival object
library(survival)

survival_object <- Surv(
  time = data$Survival_Months,
  event = data$Survival_Status
)

survival_object

# Create miR451a expression groups

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

table(data$miR451a_group)

# Kaplan-Meier model for miR451a
km_miR451a <- survfit(
  survival_object ~ miR451a_group,
  data = data
)

print(km_miR451a)

# Log-rank test
logrank_miR451a <- survdiff(
  survival_object ~ miR451a_group,
  data = data
)

print(logrank_miR451a)

# Create miR451a survival plot
plot_miR451a_survival <- ggsurvplot(
  km_miR451a,
  data = data,
  pval = TRUE,
  risk.table = TRUE,
  xlab = "Survival time (months)",
  ylab = "Survival probability",
  title = "Kaplan-Meier Survival Analysis by miR451a Expression"
)

print(plot_miR451a_survival)

# Save miR451a plot
ggsave(
  "C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech/miR451a_KM_survival.png",
  plot = plot_miR451a_survival$plot,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)

# Create miR486 expression groups
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

table(data$miR486_group)

#  Kaplan-Meier model for miR486
km_miR486 <- survfit(
  survival_object ~ miR486_group,
  data = data
)

print(km_miR486)

# Log-rank test
logrank_miR486 <- survdiff(
  survival_object ~ miR486_group,
  data = data
)

print(logrank_miR486)

# Create miR486 survival plot
plot_miR486_survival <- ggsurvplot(
  km_miR486,
  data = data,
  pval = TRUE,
  risk.table = TRUE,
  xlab = "Survival time (months)",
  ylab = "Survival probability",
  title = "Kaplan-Meier Survival Analysis by miR486 Expression"
)

print(plot_miR486_survival)

# Save miR486 plot
ggsave(
  "C:/Users/t08pb25/OneDrive - University of Aberdeen/miR486_KM_survival.png",
  plot = plot_miR486_survival$plot,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)

# Save survival group information
survival_groups <- data[, c(
  "Sample_ID",
  "Survival_Months",
  "Survival_Status",
  "miR451a",
  "miR451a_group",
  "miR486",
  "miR486_group"
)]

write.csv(
  survival_groups,
  "C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech/survival_analysis_groups.csv",
  row.names = FALSE
)

cat("Kaplan-Meier survival analysis completed successfully.\n")
  