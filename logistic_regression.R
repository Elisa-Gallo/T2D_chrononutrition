#==============================================================================#
#                                                                              #
#                  CHRONO-NUTRITIONAL BEHAVIORS AND INCIDENCE OF T2D           #
#                                 GCAT Analysis                                #
#                                                                              #                                  
#                           Elisa Gallo & Joana Llauradó                       #
#                                    May 2025                                  #
#                                                                              #
#==============================================================================#

# Description: Get regression models tables of the cohort used in this analysis.


#==============================================================================#
#                             DEFINE PARAMETERS                                #
#==============================================================================#

# Directory where outputs/results will be stored. 
results_folder <- "results/3_regressions"
current_date <- "22082025"
cohort_name <- "gcat"


#==============================================================================#
#                        INIT WORKING ENVIRONMENT                              #----
#==============================================================================#



final_df_gcat <- readRDS("data/diabetes_df_new_05052025.rds")

# Load libraries 

library(tidyverse)
library(readxl)


#==============================================================================#
#                             Define functions                                 #----
#==============================================================================#

# Function to run models
run_model <- function(dep_var, adjust_vars, df) {
  formula <- as.formula(paste("has_t2d ~", paste(c(dep_var, adjust_vars), collapse = " + ")))
  glm(formula, data = df, family = binomial)
}

# Function to extract OR, confidence intervals, P-value, and estimates
extract_results <- function(model, model_name) {
  
  # Extract model summary
  coef_results <- summary(model)$coefficients  # Estimates, SE, p-values
  conf_int_est <- confint(model)  # Confidence intervals for log-odds (estimate scale)
  conf_int_OR <- exp(conf_int_est)  # Confidence intervals for OR
  
  # Compute OR (Odds Ratios)
  OR <- exp(coef(model))
  
  # Create a formatted OR with CI column: "OR [Lower_CI, Upper_CI]"
  OR_CI <- paste0(formatC(OR, format = "f", digits = 2), 
                  " [", formatC(conf_int_OR[, 1], format = "f", digits = 2), 
                  ", ", formatC(conf_int_OR[, 2], format = "f", digits = 2), "]")
  
  # Create a results dataframe
  results <- data.frame(
    Variable = rownames(coef_results),   # Variable names (including Intercept)
    Estimate = formatC(coef_results[, "Estimate"], format = "f", digits = 3),  # Log-odds (coefficients)
    Std_Error = formatC(coef_results[, "Std. Error"], format = "f", digits = 3),  # Standard Error
    Lower_CI_Estimate = formatC(conf_int_est[, 1], format = "f", digits = 3),  # Lower CI (log-odds scale)
    Upper_CI_Estimate = formatC(conf_int_est[, 2], format = "f", digits = 3),  # Upper CI (log-odds scale)
    OR_CI = OR_CI,  # Combined OR with CI
    P_value = formatC(coef_results[, "Pr(>|z|)"], format = "f", digits = 3),  # P-value
    Model = model_name
  )
  
  # Reset row names to avoid numbered index
  rownames(results) <- NULL
  
  return(results)
}


#==============================================================================#
#                                                                              #
#     LOGISTIC REGRESSIONS FIRST MEAL, LAST MEAL AND EATING OCCASIONS          #----
#                                                                              #
#==============================================================================#


# Define all variable groups
sociodemographic <- c("edad", "sexo", "education_grouped")
lifestyle <- c("smoking_habit_2018", "mets_semana", "buckland_r_med_sum", 
               "energiat", "p_gr_g4")
health_status <- c("has_cvd")
sleep_patterns <- c( "sleeping_midpoint_avg") 
bmi <- c("bmi") 


# Exclude variables with missing values in final_df_gcat

subset_df <- final_df_gcat %>% 
  dplyr::select(c("edad", "sexo", "education_grouped", 
                  "smoking_habit_2018", "mets_semana", "buckland_r_med_sum", 
                  "energiat", "p_gr_g4", "has_cvd", "bmi", 
                  "sleeping_midpoint_avg",
                  "first_meal_avg", "last_meal_avg", "eating_occasions_avg",
                  "nighttime_fasting_avg",
                  "eating_jetlag_abs", 
                  "has_t2d", "first_meal_category",
                  "id", "first_meal_category",
                  'time_wakeup_firstmeal',
                  'time_lastmeal_gobed', "eating_midpoint_avg")) %>%
  na.omit()


#------------------------------------------------------------------------------#
#                               Run models                                     #
#------------------------------------------------------------------------------#

# Model 0: Crude models (separately for each meal variable)
model_0_first <- run_model("first_meal_avg", c("edad", "sexo"), subset_df)
model_0_last <- run_model("last_meal_avg", c("edad", "sexo"), subset_df)
model_0_occasions <- run_model("eating_occasions_avg", c("edad", "sexo"), subset_df)

# Model 1: Adjusted for sociodemographic & lifestyle
model_1_first <- run_model("first_meal_avg", c(sociodemographic, lifestyle), subset_df)
model_1_last <- run_model("last_meal_avg", c(sociodemographic, lifestyle), subset_df)
model_1_occasions <- run_model("eating_occasions_avg", c(sociodemographic, lifestyle), subset_df)

# Model 2: Model 1 + sleep patterns
model_3_first <- run_model("first_meal_avg", c(sociodemographic, lifestyle, sleep_patterns), subset_df)
model_3_last <- run_model("last_meal_avg", c(sociodemographic, lifestyle, sleep_patterns), subset_df)
model_3_occasions <- run_model("eating_occasions_avg", c(sociodemographic, lifestyle, sleep_patterns), subset_df)

# Model 3: Model 2 + BMI
model_2_first <- run_model("first_meal_avg", c(sociodemographic, lifestyle, sleep_patterns, bmi), subset_df)
model_2_last <- run_model("last_meal_avg", c(sociodemographic, lifestyle, sleep_patterns, bmi), subset_df)
model_2_occasions <- run_model("eating_occasions_avg", c(sociodemographic, lifestyle, sleep_patterns, bmi), subset_df)

# Model 4: Model 3 + Health status
model_4_first <- run_model("first_meal_avg", c(sociodemographic, lifestyle, sleep_patterns, bmi, health_status), subset_df)
model_4_last <- run_model("last_meal_avg", c(sociodemographic, lifestyle, sleep_patterns, bmi, health_status), subset_df)
model_4_occasions <- run_model("eating_occasions_avg", c(sociodemographic, lifestyle, sleep_patterns, bmi, health_status), subset_df)

# Mutually adjusted models (1b - 4b)
model_0b <- run_model(c("first_meal_avg", "last_meal_avg", "eating_occasions_avg"), c("edad", "sexo"), subset_df)
model_1b <- run_model(c("first_meal_avg", "last_meal_avg", "eating_occasions_avg"), c(sociodemographic, lifestyle), subset_df)
model_2b <- run_model(c("first_meal_avg", "last_meal_avg", "eating_occasions_avg"), c(sociodemographic, lifestyle, sleep_patterns), subset_df)
model_3b <- run_model(c("first_meal_avg", "last_meal_avg", "eating_occasions_avg"), c(sociodemographic, lifestyle, sleep_patterns, bmi), subset_df)
model_4b <- run_model(c("first_meal_avg", "last_meal_avg", "eating_occasions_avg"), c(sociodemographic, lifestyle, sleep_patterns, bmi, health_status), subset_df)


# Extract results for all models
results_list <- list(
  extract_results(model_0_first, "Model 0 - First Meal"),
  extract_results(model_0_last, "Model 0 - Last Meal"),
  extract_results(model_0_occasions, "Model 0 - Eating Occasions"),
  
  extract_results(model_1_first, "Model 1 - First Meal"),
  extract_results(model_1_last, "Model 1 - Last Meal"),
  extract_results(model_1_occasions, "Model 1 - Eating Occasions"),
  
  extract_results(model_2_first, "Model 2 - First Meal"),
  extract_results(model_2_last, "Model 2 - Last Meal"),
  extract_results(model_2_occasions, "Model 2 - Eating Occasions"),
  
  extract_results(model_3_first, "Model 3 - First Meal"),
  extract_results(model_3_last, "Model 3 - Last Meal"),
  extract_results(model_3_occasions, "Model 3 - Eating Occasions"),
  
  extract_results(model_4_first, "Model 4 - First Meal"),
  extract_results(model_4_last, "Model 4 - Last Meal"),
  extract_results(model_4_occasions, "Model 4 - Eating Occasions"),
  
  extract_results(model_0b, "Model 0b - Mutually Adjusted"),
  extract_results(model_1b, "Model 1b - Mutually Adjusted"),
  extract_results(model_2b, "Model 2b - Mutually Adjusted"),
  extract_results(model_3b, "Model 3b - Mutually Adjusted"),
  extract_results(model_4b, "Model 4b - Mutually Adjusted")
)

# Combine all results into one data frame
final_results <- bind_rows(results_list)

# Define file name
file_path <- file.path(out_folder, paste0("regressions_models_gcat_",
                                          current_date, ".csv"))

# Save results to CSV
write_csv(final_results, file = file_path)


# Print a message
print("Model results saved to 'regressions_models_gcat.csv'")




#==============================================================================#
#                                                                              #
#                   LOGISTIC REGRESSIONS NIGHTTIME FASTING                     #
#                                                                              #
#==============================================================================#

# Initialize an empty list to store results
results_list_nf <- list()

#Mutually adjusted models (nighttime fasting- eating occasions + eating midpont)
# Model 0NF: Adjusted for age and sex 
model_0NF <- run_model("nighttime_fasting_avg", c("edad", "sexo", 
                                                  "eating_occasions_avg",
                                                  "eating_midpoint_avg"), subset_df)
results_list_nf[["Model_0NF"]] <- extract_results(model_0NF, "Model_0NF")

# Model 1NF: Adjusted for sociodemographic + lifestyle 
model_1NF <- run_model("nighttime_fasting_avg", c(sociodemographic, lifestyle,
                                                  "eating_occasions_avg",
                                                  "eating_midpoint_avg"), subset_df)
results_list_nf[["Model_1NF"]] <- extract_results(model_1NF, "Model_1NF")

# Model 2NF: Model 1NF + sleep patterns
model_2NF <- run_model("nighttime_fasting_avg", 
                       c(sociodemographic, lifestyle, "eating_occasions_avg",
                         "eating_midpoint_avg", sleep_patterns), 
                       subset_df)
results_list_nf[["Model_2NF"]] <- extract_results(model_2NF, "Model_2NF")

# Model 3NF: Model 2NF + BMI 
model_3NF <- run_model("nighttime_fasting_avg", 
                       c(sociodemographic, lifestyle, "eating_occasions_avg",
                         "eating_midpoint_avg", sleep_patterns, bmi), 
                       subset_df)
results_list_nf[["Model_3NF"]] <- extract_results(model_3NF, "Model_3NF")

# Model 4NF: Model 3NF + health status
model_4NF <- run_model("nighttime_fasting_avg", 
                       c(sociodemographic, lifestyle, "eating_occasions_avg",
                         "eating_midpoint_avg", sleep_patterns, bmi, 
                         health_status), subset_df)
results_list_nf[["Model_4NF"]] <- extract_results(model_4NF, "Model_4NF")



#Single exposure

# Model 0: Adjusted for sage and sex
model_1NF_uni <- run_model("nighttime_fasting_avg", c(sociodemographic, lifestyle), subset_df)
results_list_nf[["Model_1NF_uni"]] <- extract_results(model_1NF_uni, "Model_1NF_uni")

# Model 1: Adjusted for sociodemographic + lifestyle
model_1NF_uni <- run_model("nighttime_fasting_avg", c(sociodemographic, lifestyle), subset_df)
results_list_nf[["Model_1NF_uni"]] <- extract_results(model_1NF_uni, "Model_1NF_uni")

# Model 2: Model 1 + sleep patterns
model_2NF_uni <- run_model("nighttime_fasting_avg", 
                           c(sociodemographic, lifestyle, sleep_patterns), 
                           subset_df)
results_list_nf[["Model_2NF_uni"]] <- extract_results(model_2NF_uni, "Model_2NF_uni")

# Model 3: Model 2 + BMI 
model_3NF_uni <- run_model("nighttime_fasting_avg", 
                           c(sociodemographic, lifestyle, sleep_patterns, bmi), 
                           subset_df)
results_list_nf[["Model_3NF_uni"]] <- extract_results(model_3NF_uni, "Model_3NF_uni")

# Model 4: Model 3 + health status

model_4NF_uni <- run_model("nighttime_fasting_avg", 
                           c(sociodemographic, lifestyle, sleep_patterns, bmi, 
                             health_status), subset_df)
results_list_nf[["Model_4NF_uni"]] <- extract_results(model_4NF_uni, "Model_4NF_uni")



# Combine all results into one data frame
final_results_nf <- bind_rows(results_list_nf)

# Define file path
file_path_nf <- file.path(out_folder, paste0("night_fasting_models_gcat_", current_date, ".csv"))

# Save results to CSV with a semicolon separator
write_csv(final_results_nf, file = file_path_nf)




#==============================================================================#
#                                                                              #
#                       LOGISTIC REGRESSIONS EATING JETLAG                     #----
#                                                                              #
#==============================================================================#

# Initialize an empty list to store results
results_list_ejl <- list()

#Mutually adjusted 

# Model 0: Adjusted for age sex
model_0_ejl <- run_model("eating_jetlag_abs", c("edad", "sexo", "eating_occasions_avg"), subset_df)

# Model 1: Adjusted for sociodemographic + lifestyle
model_1_ejl <- run_model("eating_jetlag_abs", c(sociodemographic, lifestyle, "eating_occasions_avg"), subset_df)

# Model 2: Model 1EJL + sleep patterns 
model_2_ejl <- run_model("eating_jetlag_abs", c(sociodemographic, lifestyle, "eating_occasions_avg", sleep_patterns), subset_df)

# Model 3: Model 2EJL + BMI
model_3_ejl <- run_model("eating_jetlag_abs", c(sociodemographic, lifestyle, "eating_occasions_avg", sleep_patterns, bmi), subset_df)

# Model 4: Model 3EJL + Health status
model_4_ejl <- run_model("eating_jetlag_abs", c(sociodemographic, lifestyle, "eating_occasions_avg", sleep_patterns, bmi, health_status), subset_df)


#Single exposure

# Model 0: Adjusted for age and sex
model_0_ejl_uni <- run_model("eating_jetlag_abs", c("edad", "sexo"), subset_df)

# Model 1: Adjusted for sociodemographic + lifestyle
model_1_ejl_uni <- run_model("eating_jetlag_abs", c(sociodemographic, lifestyle), subset_df)

# Model 2: Model 1 + sleep patterns
model_2_ejl_uni <- run_model("eating_jetlag_abs", c(sociodemographic, lifestyle, sleep_patterns), subset_df)

# Model 3: Model 2 + BMI
model_3_ejl_uni <- run_model("eating_jetlag_abs", c(sociodemographic, lifestyle, sleep_patterns, bmi), subset_df)

# Model 4: Model 3 + healt status
model_4_ejl_uni <- run_model("eating_jetlag_abs", c(sociodemographic, lifestyle, sleep_patterns, bmi, health_status), subset_df)


# Extract results for both models
results_list_ejl <- list(
  extract_results(model_0_ejl, "model_0_ejl"),
  extract_results(model_1_ejl, "model_1_ejl"),
  extract_results(model_2_ejl, "model_2_ejl"),
  extract_results(model_3_ejl, "model_3_ejl"),
  extract_results(model_4_ejl, "model_4_ejl"),
  extract_results(model_4_ejl_uni, "model_4_ejl_uni"),
  extract_results(model_3_ejl_uni, "model_3_ejl_uni"),
  extract_results(model_2_ejl_uni, "model_2_ejl_uni"),
  extract_results(model_1_ejl_uni, "model_1_ejl_uni"),
  extract_results(model_0_ejl_uni, "model_0_ejl_uni")
  
)

# Combine all results into one data frame
final_results_ejl <- bind_rows(results_list_ejl)

# Define file path for saving results
file_path_ejl <- file.path(out_folder, paste0("eating_jetlag_models_gcat_", current_date, ".csv"))

# Save results to CSV with a semicolon separator
write_csv(final_results_ejl, file = file_path_ejl)

# Print a message
print(paste("Eating jet lag model results saved to:", file_path_ejl))



#==============================================================================#
#                                                                              #
#     LOGISTIC REGRESSIONS TIME WAKEUP-FIRST MEAL, TIME LAST MEAL-BEDTIME      #----
#                                                                              #
#==============================================================================#

## TIME WAKEUP-FIRST MEAL

# Initialize an empty list to store results
results_list_wufm <- list()

#Mutually adjusted models

# Model 0: Adjusted for age sex
model_0_wufm <- run_model("time_wakeup_firstmeal", c("edad", "sexo", "eating_occasions_avg"), subset_df)

# Model 1: Adjusted for sociodemographic + lifestyle
model_1_wufm <- run_model("time_wakeup_firstmeal", c(sociodemographic, lifestyle, "eating_occasions_avg"), subset_df)

# Model 2: Model 1wufm + sleep patterns 
model_2_wufm <- run_model("time_wakeup_firstmeal", c(sociodemographic, lifestyle, "eating_occasions_avg", sleep_patterns), subset_df)

# Model 3: Model 2wufm + BMI
model_3_wufm <- run_model("time_wakeup_firstmeal", c(sociodemographic, lifestyle, "eating_occasions_avg", sleep_patterns, bmi), subset_df)

# Model 4: Model 3wufm + health status
model_4_wufm <- run_model("time_wakeup_firstmeal", c(sociodemographic, lifestyle, "eating_occasions_avg", sleep_patterns, bmi, health_status), subset_df)


# Single exposure 

# Model 0: Adjusted for age and sex
model_0_wufm_uni <- run_model("time_wakeup_firstmeal", c("edad", "sexo"), subset_df)

# Model 1: Adjusted for sociodemographic + lifestyle
model_1_wufm_uni <- run_model("time_wakeup_firstmeal", c(sociodemographic, lifestyle), subset_df)

# Model 2: Model 1 + sleep patterns
model_2_wufm_uni <- run_model("time_wakeup_firstmeal", c(sociodemographic, lifestyle, sleep_patterns), subset_df)

# Model 3: Model 2 + BMI
model_3_wufm_uni <- run_model("time_wakeup_firstmeal", c(sociodemographic, lifestyle, sleep_patterns, bmi), subset_df)

# Model 4: Model 3 + healt status
model_4_wufm_uni <- run_model("time_wakeup_firstmeal", c(sociodemographic, lifestyle, sleep_patterns, bmi, health_status), subset_df)


# Extract results for both models
results_list_wufm <- list(
  extract_results(model_0_wufm, "model_0_wufm"),
  extract_results(model_1_wufm, "model_1_wufm"),
  extract_results(model_2_wufm, "model_2_wufm"),
  extract_results(model_3_wufm, "model_3_wufm"),
  extract_results(model_4_wufm, "model_4_wufm"),
  extract_results(model_5_wufm, "model_5_wufm"),
  extract_results(model_4_wufm_uni, "model_4_wufm_uni"),
  extract_results(model_3_wufm_uni, "model_3_wufm_uni"),
  extract_results(model_2_wufm_uni, "model_2_wufm_uni"),
  extract_results(model_1_wufm_uni, "model_1_wufm_uni"),
  extract_results(model_0_wufm_uni, "model_0_wufm_uni")
)

# Combine all results into one data frame
final_results_wufm <- bind_rows(results_list_wufm)

# Define file path for saving results
file_path_wufm <- file.path(out_folder, paste0("time_wakeup_firsmeal_models_gcat_", current_date, ".csv"))

# Save results to CSV with a semicolon separator
write_csv(final_results_wufm, file = file_path_wufm)

# Print a message
print(paste("time_wakeup_firsmeal model results saved to:", file_path_wufm))


## TIME LAST MEAL-BEDTIME---------------------------------------------------------

# Initialize an empty list to store results
results_list_lmgb <- list()

# Model : Adjusted for age sex
model_0_lmgb <- run_model("time_lastmeal_gobed", c("edad", "sexo", "eating_occasions_avg"), subset_df)

# Model 1: Adjusted for sociodemographic + lifestyle
model_1_lmgb <- run_model("time_lastmeal_gobed", c(sociodemographic, lifestyle, "eating_occasions_avg"), subset_df)

# Model 2: Model 1lmgb + sleep patterns 
model_2_lmgb <- run_model("time_lastmeal_gobed", c(sociodemographic, lifestyle, "eating_occasions_avg", sleep_patterns), subset_df)

# Model 3: Model 2lmgb + BMI
model_3_lmgb <- run_model("time_lastmeal_gobed", c(sociodemographic, lifestyle, "eating_occasions_avg", sleep_patterns, bmi), subset_df)

# Model 4: Model 3lmgb + health status
model_4_lmgb <- run_model("time_lastmeal_gobed", c(sociodemographic, lifestyle, "eating_occasions_avg", sleep_patterns, bmi, health_status), subset_df)

# Model 5: Model 4lmgb + nighttime fasting
model_5_lmgb <- run_model("time_lastmeal_gobed", c(sociodemographic, lifestyle, "eating_occasions_avg", sleep_patterns, bmi, health_status, "nighttime_fasting_avg"), subset_df)

# Single exposure

# Model 0: Adjusted for age and sex
model_0_lmgb_uni <- run_model("time_lastmeal_gobed", c("edad", "sexo"), subset_df)

# Model 1: Adjusted for sociodemographic + lifestyle
model_1_lmgb_uni <- run_model("time_lastmeal_gobed", c(sociodemographic, lifestyle), subset_df)

# Model 2: Model 1 + sleep patterns
model_2_lmgb_uni <- run_model("time_lastmeal_gobed", c(sociodemographic, lifestyle, sleep_patterns), subset_df)

# Model 3: Model 2 + BMI
model_3_lmgb_uni <- run_model("time_lastmeal_gobed", c(sociodemographic, lifestyle, sleep_patterns, bmi), subset_df)

# Model 4: Model 3 + healt status
model_4_lmgb_uni <- run_model("time_lastmeal_gobed", c(sociodemographic, lifestyle, sleep_patterns, bmi, health_status), subset_df)

# Extract results for both models
results_list_lmgb <- list(
  extract_results(model_0_lmgb, "model_0_lmgb"),
  extract_results(model_1_lmgb, "model_1_lmgb"),
  extract_results(model_2_lmgb, "model_2_lmgb"),
  extract_results(model_3_lmgb, "model_3_lmgb"),
  extract_results(model_4_lmgb, "model_4_lmgb"),
  extract_results(model_5_lmgb, "model_5_lmgb"),
  extract_results(model_4_lmgb_uni, "model_4_lmgb_uni"),
  extract_results(model_3_lmgb_uni, "model_3_lmgb_uni"),
  extract_results(model_2_lmgb_uni, "model_2_lmgb_uni"),
  extract_results(model_1_lmgb_uni, "model_1_lmgb_uni"),
  extract_results(model_0_lmgb_uni, "model_0_lmgb_uni")
  
  
)

# Combine all results into one data frame
final_results_lmgb <- bind_rows(results_list_lmgb)

# Define file path for saving results
file_path_lmgb <- file.path(out_folder, paste0("time_lastmeal_gobed_models_gcat_", current_date, ".csv"))

# Save results to CSV with a semicolon separator
write_csv(final_results_lmgb, file = file_path_lmgb)

# Print a message
print(paste("lastmeal_gobed model results saved to:", file_path_lmgb))




#==============================================================================#
#                                                                              #
#                   LOGISTIC REGRESSIONS Eating midpoint                       #                                                                              #
#==============================================================================#

# Initialize an empty list to store results
results_list_em <- list()



#Single exposure

# Model 0: Adjusted for age and sex
model_0em_uni <- run_model("eating_midpoint_avg", c("edad", "sexo"), subset_df)
results_list_em[["Model_0em_uni"]] <- extract_results(model_0em_uni, "Model_0em_uni")

# Model 1: Adjusted for sociodemographic + lifestyle
model_1em_uni <- run_model("eating_midpoint_avg", c(sociodemographic, lifestyle), subset_df)
results_list_em[["Model_1em_uni"]] <- extract_results(model_1em_uni, "Model_1em_uni")

# Model 2: Model 1 + sleep patterns
model_2em_uni <- run_model("eating_midpoint_avg", c(sociodemographic, lifestyle, sleep_patterns), subset_df)
results_list_em[["Model_2em_uni"]] <- extract_results(model_2em_uni, "Model_2em_uni")

# Model 3: Model 2 + BMI
model_3em_uni <- run_model("eating_midpoint_avg", c(sociodemographic, lifestyle, sleep_patterns, bmi), subset_df)
results_list_em[["Model_3em_uni"]] <- extract_results(model_3em_uni, "Model_3em_uni")

# Model 4: Model 3 + healt status

model_4em_uni <- run_model("eating_midpoint_avg", c(sociodemographic, lifestyle, sleep_patterns, bmi, health_status), subset_df)
results_list_em[["Model_4em_uni"]] <- extract_results(model_4em_uni, "Model_4em_uni")




# Combine all results into one data frame
final_results_em <- bind_rows(results_list_em)

# Define file path
file_path_em <- file.path(out_folder, paste0("eating_midpoint_models_gcat_", current_date, ".csv"))

# Save results to CSV with a semicolon separator
write_csv(final_results_em, file = file_path_em)

