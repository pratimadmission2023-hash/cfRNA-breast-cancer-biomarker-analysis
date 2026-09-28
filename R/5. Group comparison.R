# Exploratory group comparison of cfRNA biomarkers

# Load packages
library(ggplot2)

# Load data
data <- read.csv(
  "C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech/synthetic_cfRNA_expression_clinical_data.csv",
  check.names = FALSE
)

# Check the data
head(data)

# Select biomarkers
biomarkers <- c(
  "miR451a",
  "miR486"
)

#  Check tumour grade groups
table(data$Tumour_Grade)

#  Convert tumour grade to factor
data$Tumour_Grade <- as.factor(data$Tumour_Grade)

# Summary 
summary(data[, biomarkers])

# Mean expression by tumour grade
aggregate(
  data[, biomarkers],
  by = list(Tumour_Grade = data$Tumour_Grade),
  FUN = mean,
  na.rm = TRUE
)

# miR451a group comparison plot
plot_miR451a <- ggplot(
  data,
  aes(
    x = Tumour_Grade,
    y = miR451a
  )
) +
  geom_boxplot(
    na.rm = TRUE
  ) +
  geom_jitter(
    width = 0.15,
    size = 2,
    na.rm = TRUE
  ) +
  labs(
    title = "miR451a Expression by Tumour Grade",
    x = "Tumour Grade",
    y = "miR451a expression"
  ) +
  theme_minimal()

print(plot_miR451a)

# Save figure
ggsave(
  "C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech/miR451a_tumour_grade.png",
  plot = plot_miR451a,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)

# miR486 group comparison plot
plot_miR486 <- ggplot(
  data,
  aes(
    x = Tumour_Grade,
    y = miR486
  )
) +
  geom_boxplot(
    na.rm = TRUE
  ) +
  geom_jitter(
    width = 0.15,
    size = 2,
    na.rm = TRUE
  ) +
  labs(
    title = "miR486 Expression by Tumour Grade",
    x = "Tumour Grade",
    y = "miR486 expression"
  ) +
  theme_minimal()

print(plot_miR486)

#Save figure
ggsave(
  "C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech/miR486_tumour_grade.png",
  plot = plot_miR486,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)
# Statistical comparison

# Kruskal-Wallis test for miR451a
kruskal_miR451a <- kruskal.test(
  miR451a ~ Tumour_Grade,
  data = data
)

print(kruskal_miR451a)

# Kruskal-Wallis test for miR486
kruskal_miR486 <- kruskal.test(
  miR486 ~ Tumour_Grade,
  data = data
)

print(kruskal_miR486)


# Save statistical results

results <- data.frame(
  Biomarker = c("miR451a", "miR486"),
  Test = c(
    "Kruskal-Wallis",
    "Kruskal-Wallis"
  ),
  P_value = c(
    kruskal_miR451a$p.value,
    kruskal_miR486$p.value
  )
)

write.csv(
  results,
  "C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech/group_comparison_results.csv",
  row.names = FALSE
)

cat("Group comparison analysis completed successfully.\n")

read.csv("C:/Users/t08pb25/OneDrive - University of Aberdeen/biotech/group_comparison_results.csv")

