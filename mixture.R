#==============================================================================#
#                                                                              #
#                  CHRONO-NUTRITIONAL BEHAVIORS AND INCIDENCE OF T2D           #
#                                 GCAT Analysis                                #
#                                                                              #                                  
#                                   Elisa Gallo                                #
#                                    June 2025                                #
#                                                                              #
#==============================================================================#

# Description: Get quantile g-computation models.

library(tidyverse)
library(qgcomp)


# Load data 
final_df_gcat <- readRDS(here::here("data","diabetes_df.rds"))


#Quantile g-comp model
run_noboot <- function(dep_var, adjust_vars, exposure) {
  formula <- as.formula(
    paste("has_t2d ~",
          paste(c(dep_var, adjust_vars), collapse = " + "))
  )
  
  qgcomp.glm.noboot(formula, data = subset_df,
                    expnms = exposure, family = binomial(),
                    q = 4)
}


#Using the saMe db as in the logistic regression
subset_df <- final_df_gcat %>% 
  dplyr::select(c("edad", "sexo", "education_grouped", 
                  "smoking_habit_2018", "mets_semana", "buckland_r_med_sum", 
                  "energiat", "p_gr_g4", "has_cvd", "bmi", 
                  "sleeping_midpoint_avg",
                  "first_meal_avg", "last_meal_avg", "eating_occasions_avg",
                  "nighttime_fasting_avg",
                  "eating_jetlag_abs", 
                  "has_t2d", 
                  "id", "eating_midpoint_avg",
                  'time_wakeup_firstmeal',
                  'time_lastmeal_gobed')) %>% 
  na.omit()



#Covariates

sociodemographic <- c("edad", "sexo", "education_grouped")
lifestyle <- c("smoking_habit_2018", "mets_semana", "buckland_r_med_sum", 
               "energiat", "p_gr_g4")
health_status <- c("has_cvd")
sleep_patterns <- c("sleeping_midpoint_avg") 
bmi <- c("bmi") 


#mixture variables
sleep <- c("eating_occasions_avg",
           "nighttime_fasting_avg",
           'time_wakeup_firstmeal',
           "eating_jetlag_abs",
           "time_lastmeal_gobed")


#Run models
model0 <- run_noboot(sleep, c("edad", "sexo"), sleep)
#Adjusted for sociodemographic & lifestyle
model1 <- run_noboot(sleep, c(sociodemographic, lifestyle), sleep)
#Model 1 + sleep patterns
model2 <- run_noboot(sleep, c(sociodemographic, lifestyle, sleep_patterns), sleep)
# Model 3: Model 2 + BMI
model3 <- run_noboot(sleep, c(sociodemographic, lifestyle, sleep_patterns, bmi), sleep)
# Model 4: Model 3 + Health status
model4 <- run_noboot(sleep, c(sociodemographic, lifestyle, sleep_patterns, bmi, health_status), sleep)

fun_weights <- function(model){
  
  rbind(
    data.frame(exposure = names(model$pos.weights), weight = model$pos.weights),
    data.frame(exposure = names(model$neg.weights), weight = -model$neg.weights)
  )
}

model_names <- c("model0", "model1", "model2", "model3", "model4")

#Weights of each mixture variable in each model
weights_list <- map(list(model0, model1, model2, model3, model4), fun_weights) %>%
  setNames(nm = model_names)

#OR and 95%CI of each model
model_td <- map(list(model0, model1, model2, model3, model4), broom::tidy) %>% 
  setNames(nm = model_names)


table_or <- data.frame(term = model_names,
                  OR = c(exp(model_td[["model0"]]$estimate[[2]]), 
                         exp(model_td[["model1"]]$estimate[[2]]),
                         exp(model_td[["model2"]]$estimate[[2]]),
                         exp(model_td[["model3"]]$estimate[[2]]),
                         exp(model_td[["model4"]]$estimate[[2]])),
                  
                 `Lower CI` = c(exp(model_td[["model0"]]$`Lower CI`[[2]]),
                          exp(model_td[["model1"]]$`Lower CI`[[2]]),
                          exp(model_td[["model2"]]$`Lower CI`[[2]]),
                          exp(model_td[["model3"]]$`Lower CI`[[2]]),
                          exp(model_td[["model4"]]$`Lower CI`[[2]])),
                          
                 `Upper CI` = c(exp(model_td[["model0"]]$`Upper CI`[[2]]),
                          exp(model_td[["model1"]]$`Upper CI`[[2]]),
                          exp(model_td[["model2"]]$`Upper CI`[[2]]),
                          exp(model_td[["model3"]]$`Upper CI`[[2]]),
                          exp(model_td[["model4"]]$`Upper CI`[[2]])))



