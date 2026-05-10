# ==========================================
# Survival Analysis: WPBC Dataset
# ==========================================

#  Load necessary libraries
# If not installed, use: install.packages(c("survival", "survminer", "dplyr"))
library(survival)
library(survminer)
library(dplyr)

# Load the dataset
# Replace the path below with your local path
#file_path <- "D:/githhub_project/Survival Analysis for Breast Cancer/data/wpbc.data"
file_path <- "data/wpbc.data"
data <- read.csv(file_path, header = FALSE, na.strings = "?")

# Assign proper column names
feature_names <- c("radius", "texture", "perimeter", "area", "smoothness", 
                   "compactness", "concavity", "concave_points", "symmetry", "fractal_dim")

cols <- c("ID", "Outcome", "Time", 
          paste0(rep(feature_names, each=3), c("_mean", "_se", "_worst")),
          "TumorSize", "LymphNodes")
colnames(data) <- cols

#  Data Cleaning & Pre-processing
# 1. Remove rows with missing LymphNodes (standard for this dataset)
data <- data %>% filter(!is.na(LymphNodes))

# 2. Convert Outcome to numeric status (R = 1, N = 0)
data$status <- ifelse(data$Outcome == "R", 1, 0)

# 3. Create a categorical variable for grouping (e.g., Tumor Size > Median)
# This replaces the 'treatment' variable from your original plan
data$Size_Group <- ifelse(data$TumorSize > median(data$TumorSize), "Large", "Small")

# Create the Survival Object
# time = Time to recurrence or censoring
# event = status (1 for recurrence)
surv_object <- Surv(time = data$Time, event = data$status)

#  Kaplan-Meier Survival Curves
# Overall Survival
km_fit_total <- survfit(surv_object ~ 1, data = data)
ggsurvplot(km_fit_total, data = data, title = "Overall Recurrence-Free Survival")

# Comparison by Tumor Size Group
km_fit_size <- survfit(surv_object ~ Size_Group, data = data)
ggsurvplot(
  km_fit_size, 
  data = data,
  pval = TRUE,             # Log-Rank test p-value
  risk.table = TRUE,       # Show number at risk
  conf.int = TRUE,         # Show confidence intervals
  palette = c("#E7B800", "#2E9FDF"),
  xlab = "Time (Months)",
  legend.title = "Tumor Size",
  title = "Survival by Tumor Size (Median Split)"
)

#  Log-Rank Test
# Statistical test to see if curves are significantly different
surv_diff <- survdiff(surv_object ~ Size_Group, data = data)
print(surv_diff)

#  Cox Proportional Hazards Model
# We test if TumorSize and LymphNodes are significant predictors
cox_model <- coxph(surv_object ~ TumorSize + LymphNodes + radius_mean, data = data)
summary(cox_model)

#  Cox Model Diagnostics
# Test the Proportional Hazards assumption
test_ph <- cox.zph(cox_model)
print(test_ph)
ggcoxzph(test_ph) # If lines are relatively flat, the model is valid

#  Adjusted Survival Curves
# Shows predicted survival for an average patient in each group
ggadjustedcurves(cox_model, data = data, variable = "Size_Group")