######
# link http://www.sthda.com/english/wiki/survival-analysis-basics
library(tidyverse)
#install.packages(c("survival", "survminer"))
library(survival)
library(survminer)
library(readr)
library(readxl)
library(lubridate)
#install.packages("survminer")
df <- read_tsv("/Users/kumarr9/Downloads/survival.tsv")
# Create a 'Surv' object
surv_obj <- with(df, Surv(as.numeric(difftime(as.Date(`date of death`, "%m/%d/%y"), as.Date(`Date of Diagnosis`, "%m/%d/%y")), units = "days")))
# Fit the survival model
fit <- survfit(surv_obj ~ Cluster, data = df)
# Plot the Kaplan-Meier survival curve
ggsurvplot(fit, data = df, risk.table = TRUE, pval = TRUE, ncensor.plot = TRUE,surv.median.line = "hv",title = "Kaplan-Meier Survival Curve")

###pairwise difference
pairwise_survdiff(
  formula = Surv(as.numeric(difftime(as.Date(`date of death`, "%m/%d/%y"), as.Date(`Date of Diagnosis`, "%m/%d/%y")), units = "days")) ~ Cluster, 
  data = df, 
  p.adjust.method = "fdr"
)

### 2nd try
# Create a 'Surv' object
surv_obj <- with(df, Surv(as.numeric(difftime(as.Date(`date of relapse`, "%m/%d/%y"), as.Date(`Date of Diagnosis`, "%m/%d/%y")), units = "days")))

# Fit the Cox proportional hazards model
cox_model <- coxph(surv_obj ~ Cluster, data = df)

# Extract hazard ratios
hr_summary <- summary(cox_model)

# Print hazard ratios
cat("\nHazard Ratios:\n")
print(hr_summary)

# Plot the Kaplan-Meier survival curve with hazard ratios
ggsurvplot(survfit(surv_obj ~ Cluster, data = df), data = df, risk.table = TRUE, pval = TRUE,
           title = "Kaplan-Meier Survival Curve with Hazard Ratios",
           hr.at = c(1, 2, 3), # Specify which clusters to compare
           conf.int = TRUE)   # Show confidence intervals for hazard ratios

### 3rd try
df$`Date of Diagnosis` <- as.Date(df$`Date of Diagnosis`, format = "%m/%d/%y")
df$`date of relapse` <- as.Date(df$`date of relapse`, format = "%m/%d/%y")

# Create a 'Surv' object
surv_obj <- with(df, Surv(as.numeric(difftime(`date of relapse`, `Date of Diagnosis`, units = "days")),
                          event = ifelse(!is.na(`date of relapse`), 1, 0)))

# Fit the Cox proportional hazards model
cox_model <- coxph(surv_obj ~ Cluster, data = df)

# Extract hazard ratios
hr_summary <- summary(cox_model)

# Print hazard ratios
cat("\nHazard Ratios:\n")
print(hr_summary)

# Plot the Kaplan-Meier survival curve with censored observations and hazard ratios
km_plot <- ggsurvplot(survfit(surv_obj ~ Cluster, data = df), data = df, risk.table = TRUE, pval = TRUE,
                      title = "Kaplan-Meier Survival Curve with Hazard Ratios",
                      hr.at = c(1, 2, 3), # Specify which clusters to compare
                      conf.int = TRUE)   # Show confidence intervals for hazard ratios

# Combine the survival plot and risk table
combined_plot <- ggarrange(km_plot$plot, km_plot$table, ncol = 1, heights = c(0.8, 0.2))

plot(combined_plot)

#### trying multivariate analysis
# Create a 'Surv' object for time-to-event data
df<- read_tsv("/Users/kumarr9/Downloads/survival_multivariate.tsv")
# Create a 'Surv' object for time-to-event data
surv_obj <- with(df, Surv(as.numeric(difftime(as.Date(date_of_death, "%m/%d/%y"), as.Date(Date_of_Diagnosis, "%m/%d/%y")), units = "days")))

# Fit the Cox Proportional-Hazards Model
cox_model <- coxph(surv_obj ~ Cluster + Stage_at_Diagnosis + Days_btwn_Diagnosis_Bx + Platinum_S_R, data = df)

# Visualize the results of the Cox model using ggforest
ggforest(cox_model, data = df) +
  labs(title = "Cox Proportional-Hazards Model Results")

cox_model <- coxph(surv_obj ~ Cluster + Stage_at_Diagnosis + Days_btwn_Diagnosis_Bx + Age_Dx  + Platinum_S_R, data = df)

# Visualize the results of the Cox model using ggforest
ggforest(cox_model, data = df) +
  labs(title = "Cox Proportional-Hazards Model Results by Cluster")




#### testing and trying code ####
library(tidyverse)
df <- read_tsv("/Users/kumarr9/Downloads/survival_new.tsv")
# Convert date columns to Date format
# Convert date columns to Date format
df$Date_of_Diagnosis <- as.Date(df$Date_of_Diagnosis, format = "%m/%d/%y")
df$date_of_death <- as.Date(df$date_of_death, format = "%m/%d/%y")

# Create a survival object
surv_obj <- with(df, Surv(time = as.numeric(difftime(date_of_death, Date_of_Diagnosis, units = "days")), 
                          event = as.numeric(date_of_death != "")))

# Fit a survival model
fit <- survfit(surv_obj ~ Cluster, data = df)
###
df <- read_tsv("/Users/kumarr9/Downloads/latency.tsv")
ggplot(df, aes(x = growth_rate, fill = Cluster)) +
  geom_density(alpha = 0.4) +
  labs(x = "Average Rate of Growth (mm3/day)", y = "Density") +
  # facet_wrap(~Cluster, scales = "free_y") +
  theme_minimal()

ggplot(df, aes(x = latency, fill = Cluster)) +
  geom_density(alpha = 0.4) +
  labs(x = "Average Latency (days)", y = "Density") +
  # facet_wrap(~Cluster, scales = "free_y") +
  theme_minimal()

#### Survival analysis January 30, 2024
#data <- read_excel("/Users/kumarr9/Downloads/survival_new_jan_30.xlsx", sheet = 2)
data <- read_tsv("/Users/kumarr9/Downloads/survival_new_jan_30.tsv")
#colnames(data)
km_trt_fit <- survfit(Surv(Survival_since_diagnosis, Event) ~ Cluster, data=data)
ggsurvplot(km_trt_fit, data = data, risk.table = TRUE, pval = TRUE, ncensor.plot = TRUE,surv.median.line = "hv",title = "Kaplan-Meier Survival Curve")

## first 250 days survival 
data$risk <- ifelse(data$Survival_since_diagnosis < 400, "High Risk", "Low Risk")
surv_obj <- Surv(time = data$Survival_since_diagnosis, event = data$Event)
#km_trt_fit <- survfit(Surv(Survival_since_diagnosis, Event) ~ Cluster, data=data)
km_fit <- survfit(surv_obj ~ Cluster + risk, data = data)
ggsurvplot(km_fit, data = data, risk.table = TRUE, pval = TRUE, ncensor.plot = TRUE,surv.median.line = "hv",title = "Kaplan-Meier Survival Curve")

## for high risk only

# Subset high-risk values
high_risk_data <- data[data$risk == "High Risk", ]

# Create a survival object
surv_obj_high_risk <- Surv(time = high_risk_data$Survival_since_diagnosis, event = high_risk_data$Event)

# Fit survival curves and create a Kaplan-Meier plot for high-risk values
km_fit_high_risk <- survfit(surv_obj_high_risk ~ Cluster, data = high_risk_data)
ggsurvplot(km_fit_high_risk, data = high_risk_data, pval = TRUE)
## mutivariate
for (variable in colnames(high_risk_data)[!colnames(high_risk_data) %in% c("Survival_since_diagnosis", "Sex", "Average_Rate_of_Growth", "Average_Latency", "Event", "risk", "Cluster")]) {
  multivariate_fit <- survfit(Surv(Survival_since_diagnosis, Event) ~ get(variable), data = high_risk_data)
  
  # Create a multivariate survival plot
  ggsurvplot(multivariate_fit, data = high_risk_data, title = variable, pval = TRUE)
}
 ## building random forest
install.packages("ranger")
library(ranger)
r_fit <- ranger(Surv(time = high_risk_data$Survival_since_diagnosis, event = high_risk_data$Event) ~ Cluster + Passage + 
                  Sex + Average_Rate_of_Growth + Average_Latency + Age_Bx,
                data = high_risk_data,
                mtry = 4,
                importance = "permutation",
                splitrule = "extratrees",
                verbose = TRUE)

# Average the survival models
death_times <- r_fit$unique.death.times 
surv_prob <- data.frame(r_fit$survival)
avg_prob <- sapply(surv_prob,mean)
# Plot the survival models for each patient
plot(r_fit$unique.death.times,r_fit$survival[1,], 
     type = "l", 
     ylim = c(0,1),
     col = "red",
     xlab = "Days",
     ylab = "survival",
     main = "Patient Survival Curves")


#
cols <- colors()
for (n in sample(c(2:dim(vet)[1]), 20)){
  lines(r_fit$unique.death.times, r_fit$survival[n,], type = "l", col = cols[n])
}
lines(death_times, avg_prob, lwd = 2)
legend(500, 0.7, legend = c('Average = black'))


#####

