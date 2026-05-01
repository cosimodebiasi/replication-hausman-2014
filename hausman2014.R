# Replication of the microeconometric study
# "Corporate Incentives and Nuclear Safety"
# by C. Hausman (2014)
  
library(haven)
library(dplyr)
library(readr)
library(tidyr)    
library(lmtest)
library(sandwich)

# Load data

#---Table 1: Annual reactor-level summary statistics---#
    
# Divide "personrem" by the number of units and convert "generation" to million MWh
table1_data <- data %>%
  mutate(
    personrem = ifelse(is.na(number_units) | number_units == 0, NA, personrem / number_units),
    generation = generation / 10^6
    )
  
# Collapse the data by "reactor_name" and "year" to get annual, reactor-specific measures
collapsed_data <- table1_data %>%
  group_by(reactor_name, year) %>%
  summarise(
    initiating = sum(initiating, na.rm = TRUE),
    fire = sum(fire, na.rm = TRUE),
    escalated = sum(escalated, na.rm = TRUE),
    generation = sum(generation, na.rm = TRUE),
    rem = mean(rem, na.rm = TRUE),
    personrem = mean(personrem, na.rm = TRUE),
    capacity_factor = mean(capacity_factor, na.rm = TRUE)
    )
  
collapsed_data <- collapsed_data %>%
  mutate(
    initiating = ifelse(year < 1988, NA, initiating),
    fire = ifelse(year < 1991, NA, fire),
    escalated = ifelse(year < 1996, NA, escalated)
    )
  
general_stats <- collapsed_data %>%
  filter(year >= 1996) %>%
  ungroup() %>%
  summarise(
    mean_initiating = mean(initiating, na.rm = TRUE),
    sd_initiating = sd(initiating, na.rm = TRUE),
    min_initiating = min(initiating, na.rm = TRUE),
    max_initiating = max(initiating, na.rm = TRUE),
        
    mean_fire = mean(fire, na.rm = TRUE),
    sd_fire = sd(fire, na.rm = TRUE),
    min_fire = min(fire, na.rm = TRUE),
    max_fire = max(fire, na.rm = TRUE),
     
    mean_escalated = mean(escalated, na.rm = TRUE),
    sd_escalated = sd(escalated, na.rm = TRUE),
    min_escalated = min(escalated, na.rm = TRUE),
    max_escalated = max(escalated, na.rm = TRUE),
       
    mean_generation = mean(generation, na.rm = TRUE),
    sd_generation = sd(generation, na.rm = TRUE),
    min_generation = min(generation, na.rm = TRUE),
    max_generation = max(generation, na.rm = TRUE),
     
    mean_capacity_factor = mean(capacity_factor, na.rm = TRUE),
    sd_capacity_factor = sd(capacity_factor, na.rm = TRUE),
    min_capacity_factor = min(capacity_factor, na.rm = TRUE),
    max_capacity_factor = max(capacity_factor, na.rm = TRUE),
      
    mean_rem = mean(rem, na.rm = TRUE),
    sd_rem = sd(rem, na.rm = TRUE),
    min_rem = min(rem, na.rm = TRUE),
    max_rem = max(rem, na.rm = TRUE),
       
    mean_personrem = mean(personrem, na.rm = TRUE),
    sd_personrem = sd(personrem, na.rm = TRUE),
    min_personrem = min(personrem, na.rm = TRUE),
    max_personrem = max(personrem, na.rm = TRUE)
    )
 
summary_stats_table <- data.frame(
  measure = c("Initiating events", "Fires", "Collective worker radiation exposure (person-rems)", "Average worker radiation exposure (rems)", "Escalated enforcement", "Generation (million MWh)", "Capacity factor"),
  mean = c(general_stats$mean_initiating, general_stats$mean_fire, general_stats$mean_personrem, general_stats$mean_rem, general_stats$mean_escalated, general_stats$mean_generation, general_stats$mean_capacity_factor),
  sd = c(general_stats$sd_initiating, general_stats$sd_fire, general_stats$sd_personrem, general_stats$sd_rem, general_stats$sd_escalated, general_stats$sd_generation, general_stats$sd_capacity_factor),
  min = c(general_stats$min_initiating, general_stats$min_fire, general_stats$min_personrem, general_stats$min_rem, general_stats$min_escalated, general_stats$min_generation, general_stats$min_capacity_factor),
  max = c(general_stats$max_initiating, general_stats$max_fire, general_stats$max_personrem, general_stats$max_rem, general_stats$max_escalated, general_stats$max_generation, general_stats$max_capacity_factor)
  )
  
summary_stats_table <- summary_stats_table %>%
  mutate(across(.cols = -measure, .fns = ~ round(.x, 2)))
 
print(summary_stats_table)

#---Table 2: Comparing divested and nondivested nuclear reactors---#

table2_data <- data %>%
  mutate(
    personrem = ifelse(is.na(number_units) | number_units == 0, NA, personrem / number_units),
    generation = generation / 10^6
    )

# Collapse the data to the reactor-year-divestiture status level
table2_data <- table2_data %>%
  group_by(reactor_name, facilityname, year, ever_divested, plantid) %>%
  summarise(
    initiating = sum(initiating, na.rm = TRUE),
    fire = sum(fire, na.rm = TRUE),
    escalated = sum(escalated, na.rm = TRUE),
    generation = sum(generation, na.rm = TRUE),
    rem = mean(rem, na.rm = TRUE),
    personrem = mean(personrem, na.rm = TRUE),
    capacity_factor = mean(capacity_factor, na.rm = TRUE),
    .groups = "drop"
    )

# Radiation exposure variables are measured by plant, so I create a new variable, "byplant",
  # which takes the value 1 for the first "facilityname-year" combination and 0 for the rest.
  # This variable is particularly useful for multi-unit plants, as it helps us pick
  # only one reactor per plant.
table2_data <- table2_data %>%
  group_by(facilityname, year) %>%
  mutate(byplant = ifelse(row_number() == 1, 1, 0)) %>%
  ungroup()

results <- data.frame(
  Variable = character(),
  Mean_NeverDivested = numeric(),
  Mean_LaterDivested = numeric(),
  T_Statistic = numeric(),
  P_Value = numeric()
  )

variables <- c("initiating", "fire", "escalated", "generation", "capacity_factor", "personrem", "rem")
 
# Calculate means, t-statistics, and p-values
for (y in variables) {
  filtered_data <- table2_data %>% filter(year >= 1996 & year <= 1998)
  
  # Filter for radiation exposure variables
  if (y %in% c("personrem", "rem")) {
    filtered_data <- filtered_data %>% filter(byplant == 1)
    }

  # Group means
  mean_never_divested <- mean(filtered_data[[y]][filtered_data$ever_divested == 0], na.rm = TRUE)
  mean_later_divested <- mean(filtered_data[[y]][filtered_data$ever_divested == 1], na.rm = TRUE)
  
  # t-statistics and p-values, clustering the s.e. at the plant level
  model <- lm(as.formula(paste(y, "~ ever_divested")), data = filtered_data)
  coeftest_result <- coeftest(model, vcov = vcovCL(model, ~ plantid))
  t_stat <- coeftest_result["ever_divested", "t value"]
  p_value <- coeftest_result["ever_divested", "Pr(>|t|)"]
 
  results <- rbind(
    results,
    data.frame(
      Variable = y,
      Mean_NeverDivested = round(mean_never_divested, 3),
      Mean_LaterDivested = round(mean_later_divested, 3),
      T_Statistic = round(-t_stat, 2),
      P_Value = round(p_value, 2)
      )
  )
  }

print(results)

# My estimates for the t-statistics are slightly different in module; this is most likely
  # because the author used a different formula for clustering. This same comment also
  # applies to subsequent regression estimates, in case of different s.e.

#---Table 3: The effect of divestiture on nuclear power plant safety---#

# All the methodological justifications are contained in my submitted homework paper.

library(MASS) # for negative binomial regression 

# Convert "generation" to million MWh
table3_data <- data %>% mutate(generation = generation / 10^6)

# As usual, collapse the data
table3_data <- table3_data %>%
  group_by(reactor_name, reactorid, year, facilityname, plantid, pwr, year_start, capacity, ever_divested, state_name, censusregion, year_divest) %>%
  summarise(initiating = sum(initiating, na.rm = TRUE),
            fire = sum(fire, na.rm = TRUE),
            escalated = sum(escalated, na.rm = TRUE),
            capacity_factor = mean(capacity_factor, na.rm = TRUE),
            divested = mean(divested, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(initiating = ifelse(year < 1988, NA, initiating),
         fire = ifelse(year < 1991, NA, fire),
         escalated = ifelse(year < 1996, NA, escalated),
         capacity_factor = ifelse(capacity_factor < 0.01, NA, capacity_factor)
         )

# Regressions of the count variables (columns 1-3 of Table 3)
results <- data.frame()
results_exp <- data.frame()

for (var in c("initiating", "fire", "escalated")) {
  # "Baseline" negative binomial regression, accounting for reactor-level fixed and time effects
  formula <- as.formula(paste(var, "~ divested + factor(year) + factor(reactor_name)"))
  model <- glm.nb(formula, data = table3_data)
  coef_test <- coeftest(model, vcov = vcovCL(model, ~ plantid)) 
  divested_base_coef <- coef_test["divested", ]
  pct_change <- exp(divested_base_coef["Estimate"]) - 1
  results <- rbind(
    results,
    data.frame(
      variable = var,
      estimate = round(divested_base_coef["Estimate"], 3),
      std_error = round(divested_base_coef["Std. Error"], 3),
      t_stat = round(divested_base_coef["Estimate"] / divested_base_coef["Std. Error"], 3),
      p_value = round(2 * (1 - pnorm(abs(divested_base_coef["Estimate"] / divested_base_coef["Std. Error"]))), 3),
      percent_change = round(pct_change, 2) * 100,
      n_obs = nobs(model)
      )
    )
  # Neg. bin. regression with dependent variable normalized by capacity factor:
    # we include "capacity_factor" as an exposure variable and offset it (i.e., the coefficient
    # of the logged variable is fixed at 1). By doing so, capacity factor is simply "fixed,"
    # and we make sure that the dependent variable responds to "divested" *only* for a given
    # level of capacity factor (and so, for a given level of generation).
  formula_exp <- as.formula(paste(var, "~ divested + factor(year) + factor(reactor_name) + offset(log(capacity_factor))"))
  model_exp <- glm.nb(formula_exp, data = table3_data)
  coef_test_exp <- coeftest(model_exp, vcov = vcovCL(model_exp, ~ plantid))
  divested_exp_coef <- coef_test_exp["divested", ]
  pct_change_exp <- exp(divested_exp_coef["Estimate"]) - 1
  results_exp <- rbind(
    results_exp,
    data.frame(
      variable = var,
      estimate = round(divested_exp_coef["Estimate"], 3),
      std_error = round(divested_exp_coef["Std. Error"], 3),
      t_stat = round(divested_exp_coef["Estimate"] / divested_exp_coef["Std. Error"], 3),
      p_value = round(2 * (1 - pnorm(abs(divested_exp_coef["Estimate"] / divested_exp_coef["Std. Error"]))), 3),
      percent_change = round(pct_change_exp, 2) * 100,
      n_obs = nobs(model_exp)
      )
  )
}

print(results)
print(results_exp)

# Statistical significance of the "divestiture" effect on unsafe events:
  # initiating: p = 0.14 ---> not significant
  # fire: p = 0.15 ---> not significant
  # p = 0.08 ---> significant at 10%
  # p = 0.006 ---> significant at 1%
  # p = 0.08 ---> significant at 10%
  # p = 0.046 ---> significant at 5%
	
# The same exact results would be given by: (e.g., for "initiating"-- non-normalized)

  # library(fixest)
   
  # a <- fenegbin(
  #     initiating ~ divested | reactor_name + year,  
  #     data = table3_data,
  #     vcov = ~plantid 
  # )
  # summary(a)

# The data setup is ideal for implementing the Callaway and Sant'Anna estimator for
  # staggered DID (reactors are generally divested at different points in time).
  # Unfortunately, the negative binomial specification and the clustering at plant level
  # prevent us from doing so.

# Regressions of the continuous variables, "personrem" and "rem" (columns 4 and 5 in Table 3)

library(plm)

table3bis_data <- data %>%
  group_by(facilityname, year, plantid, pwr, year_divest, state_name, censusregion) %>%
  summarize(
    personrem = mean(personrem, na.rm = TRUE),
    rem = mean(rem, na.rm = TRUE),
    divested = mean(divested, na.rm = TRUE),
    capacity_factor = mean(capacity_factor, na.rm = TRUE),
    ever_divested = mean(ever_divested, na.rm = TRUE),
    capacity = mean(capacity, na.rm = TRUE),
    .groups = "drop"
    ) %>%
  mutate(
    capacity_factor = ifelse(capacity_factor < 0.01, NA, capacity_factor),
    # Normalized variables
    norm_personrem = personrem / capacity_factor,
    norm_rem = rem / capacity_factor
    )
# We're going to perform simple OLS: we can use plm().
  # I also tried to simply run lm() including the fixed effects using
  # factor(plantid) and factor(year) --without the panel structure--
  # and it returns the same exact results.
table3bis_data <- pdata.frame(table3bis_data, index = c("plantid", "year"))

variables_bis <- c("personrem", "rem", "norm_personrem", "norm_rem")
results_bis <- data.frame()

for (var in variables_bis) {
  formula <- as.formula(paste(var, "~ divested"))
  # Two-way fixed effects model (plant-level fixed effects and year time effects)
  model_bis <- plm(formula, data = table3bis_data, model = "within", effect = "twoways")
  clustered_bis <- coeftest(model_bis, vcov = vcovHC(model_bis, type = "HC1", cluster = "group"))
     
  divested_coef <- clustered_bis["divested", ]
        
  # Change in expected value (percent)
  pred <- predict(model_bis)
  avg_pred <- mean(pred[table3bis_data$divested == 0], na.rm = TRUE)
  pct_change <- (divested_coef["Estimate"] / avg_pred) * 100

  results_bis <- rbind(
    results_bis,
    data.frame(
      variable = var,
      estimate = round(divested_coef["Estimate"], 3),
      std_error = round(divested_coef["Std. Error"], 3),
      p_value = round(divested_coef["Pr(>|t|)"], 3),
      percent_change = floor(pct_change),
      n_obs = nobs(model_bis)
      )
  )
}

print(results_bis)

# No coefficient is statistically significant.

# Percentage change in expected value for continuous variables:
  # beta / E[y_it | d_it = 0; alpha_i, v_t]

# The same exact results would be given by: (e.g., for "rem"-- non-normalized)

  # b <- feols(
  #     rem ~ divested | plantid + year,  
  #     data = table3bis_data,
  #     cluster = ~plantid 
  # )
  # summary(b)

#---Table 4: Robustness Checks---#

library(fixest)
library(broom) 
library(glmmTMB) 

robustness_data <- data %>%
  group_by(reactor_name, reactorid, year, facilityname, plantid, pwr, year_start,
           capacity, ever_divested, state_name, censusregion, year_divest) %>%
  summarise(
    initiating = sum(initiating, na.rm = TRUE),
    fire = sum(fire, na.rm = TRUE),
    escalated = sum(escalated, na.rm = TRUE),
    capacity_factor = mean(capacity_factor, na.rm = TRUE),
    divested = mean(divested, na.rm = TRUE),
    same_owner = mean(same_owner, na.rm = TRUE),
    coowneddum = mean(coowneddum, na.rm = TRUE),
    .groups = "drop"
    )

robustness_data <- robustness_data %>%
  mutate(
    initiating = ifelse(year < 1988, NA, initiating),
    fire = ifelse(year < 1991, NA, fire),
    escalated = ifelse(year < 1996, NA, escalated),
    capacity_factor = ifelse(capacity_factor < 0.01, NA, capacity_factor)
    )

# We will need this later to deal with zero-events reactors
for (v in c("initiating", "fire", "escalated")) {
  robustness_data <- robustness_data %>%
    group_by(reactorid) %>%
    mutate(!!paste0("zero", v) := sum(.data[[v]], na.rm = TRUE)) %>%
    ungroup()
  }

robustness_data <- robustness_data %>%
  mutate(capacity_factor = ifelse(is.na(capacity_factor) | capacity_factor <= 0, 1, capacity_factor))

# Generate fixed effects variables
for (v in c("initiating", "fire", "escalated")) {
  robustness_data <- robustness_data %>%
    mutate(!!paste0("R", v) := as.factor(reactor_name),
           !!paste0("Y", v) := as.factor(year))
}

# I create this function to store regression results (much faster !)-- for count specifications
store <- function(model, var) {
  coef_summary <- tidy(model, conf.int = TRUE) %>% filter(term == "divested")
  pct_change <- exp(coef_summary$estimate) - 1
  data.frame(
    variable = var,
    estimate = round(coef_summary$estimate, 2),
    std_error = round(coef_summary$std.error, 2),
    t_stat = round(coef_summary$statistic, 2),
    p_value = round(coef_summary$p.value, 3),
    percent_change = floor(pct_change * 100)
  )
}

results <- list()

# Poisson regression
for (v in c("initiating", "fire", "escalated")) {
  model_p <- feglm(as.formula(paste(v, "~ divested +", paste0("Y", v), "+", paste0("R", v), "| year")),
                   data = robustness_data %>% filter(get(paste0("zero", v)) != 0),
                   family = poisson, cluster = ~plantid)
  results[[paste0(v, "_poisson")]] <- store(model_p, v)
  }

# Poisson regression with normalized variables
for (v in c("initiating", "fire", "escalated")) {
  model_p_norm <- feglm(as.formula(paste(v, "~ divested +", paste0("Y", v), "+", paste0("R", v), "+ offset(log(capacity_factor))", "| year")),
                        data = robustness_data %>% filter(get(paste0("zero", v)) != 0),
                        family = poisson, cluster = ~plantid)
  results[[paste0(v, "_poisson_norm")]] <- store(model_p_norm, v)
  }

# Negative binomial regression
for (v in c("initiating", "fire", "escalated")) {
  nb <- glm.nb(as.formula(paste(v, "~ divested +", paste0("Y", v), "+", paste0("R", v))),
               data = robustness_data %>% filter(get(paste0("zero", v)) != 0))
  model_nb <- coeftest(nb, vcov = vcovCL(nb, cluster = ~plantid))
  results[[paste0(v, "_nb")]] <- store(model_nb, v)
  }

# Negative binomial regression with normalized variables
for (v in c("initiating", "fire", "escalated")) {
  nb_norm <- glm.nb(as.formula(paste(v, "~ divested +", paste0("Y", v), "+", paste0("R", v), "+ offset(log(capacity_factor))")), 
                    data = robustness_data %>% filter(get(paste0("zero", v)) != 0))
  model_nb_norm <- coeftest(nb_norm, vcov = vcovCL(nb_norm, cluster = ~plantid))
  results[[paste0(v, "_nb_norm")]] <- store(model_nb_norm, v)
}

# Function to compute expected percentage change for OLS
store_ols <- function(model, var, data) {
  coef_summary <- tidy(model, conf.int = TRUE) %>% filter(term == "divested")
  pred <- predict(model)
  avg_pred <- mean(pred[data$divested == 0], na.rm = TRUE)
  pct_change <- (coef_summary$estimate / avg_pred) * 100
  data.frame(
    variable = var,
    estimate = round(coef_summary$estimate, 2),
    std_error = round(coef_summary$std.error, 2),
    t_stat = round(coef_summary$statistic, 2),
    p_value = round(coef_summary$p.value, 3),
    percent_change = floor(pct_change)
  )
}

# OLS regression
for (v in c("initiating", "fire", "escalated")) {
  model_ols <- feols(as.formula(paste(v, "~ divested +", paste0("Y", v), "+", paste0("R", v))), 
                     data = robustness_data, cluster = ~plantid)
  results[[paste0(v, "_ols")]] <- store_ols(model_ols, v, robustness_data)
}

# OLS regression with normalized variables
for(v in c("initiating", "fire", "escalated")) {
  robustness_data <- robustness_data %>%
    mutate(!!paste0("norm_", v) := .data[[v]] / capacity_factor)
  
  model_ols_norm <- feols(as.formula(paste0("norm_", v, " ~ divested + Y", v, " + R", v)),
                          data = robustness_data, cluster = ~plantid)
  results[[paste0(v, "_ols_norm")]] <- store_ols(model_ols_norm, v, robustness_data)
  }

final_results <- bind_rows(results, .id = "model")
print(final_results)

#---Table 5: Heterogeneity by reactor charateristics---#

library(car)

# Count variables (columns 1-3 in Table 5)

table5_data <- data %>%
  group_by(reactor_name, reactorid, year, facilityname, plantid, pwr, year_start, capacity, ever_divested, state_name, censusregion, year_divest) %>%
  summarise(initiating = sum(initiating, na.rm = TRUE),
            fire = sum(fire, na.rm = TRUE),
            escalated = sum(escalated, na.rm = TRUE),
            capacity_factor = mean(capacity_factor, na.rm = TRUE),
            divested = mean(divested, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(initiating = ifelse(year < 1988, NA, initiating),
         fire = ifelse(year < 1991, NA, fire),
         escalated = ifelse(year < 1996, NA, escalated),
         capacity_factor = ifelse(capacity_factor < 0.01, NA, capacity_factor))

# Some new variables (bwr, small, new...) and their interaction with "divested"
table5_data <- table5_data %>%
  mutate(
    div_pwr = divested * pwr,
    div_bwr = divested * (1 - pwr),
    old = ifelse(year_start <= 1979, 1, 0),
    new = ifelse(year_start > 1979, 1, 0),
    div_new = divested * new,
    div_old = divested * old,
    small = ifelse(capacity < 1000, 1, 0),
    large = ifelse(capacity >= 1000, 1, 0),
    div_small = divested * small,
    div_large = divested * large
  )

# There are three regressions for each dependent variable:
# one with "divestiture BWR" and "d. PWR," one with "divestiture older" and "d. newer,"
# and one with "divestiture small" and "d. large."
# For a matter of convenience, I will call these three models "A", "B" and "C", respectively.

library(stargazer)

models_A <- list()
models_B <- list()
models_C <- list()

chisq_results <- list()

for (var in c("initiating", "fire", "escalated")) {
  
  # Model A: y = div_brw + div_pwr + v_t + alpha_i (neg. bin.)
  formula_A <- as.formula(paste(var, "~ div_bwr + div_pwr + factor(year) + factor(reactor_name)"))
  model_A <- glm.nb(formula_A, data = table5_data)
  coeftest_A <- coeftest(model_A, vcov = vcovCL(model_A, ~ plantid)) 
  
  models_A[[var]] <- round(coeftest_A, 2)
  
  # Chi-square statistics for each regression: we compare the two estimated average effects
  chisq_A <- linearHypothesis(model_A, "div_bwr = div_pwr", vcov = vcovCL(model_A, ~ plantid), test = "Chisq")
  chisq_results[[paste(var, "bwr_vs_pwr", sep = "_")]] <- list(
    chisq_stat = chisq_A$Chisq[2],
    p_value = chisq_A$`Pr(>Chisq)`[2]
  )
  
  # Model B: y = div_old + div_new + v_t + alpha_i (neg. bin.)
  formula_B <- as.formula(paste(var, "~ div_old + div_new + factor(year) + factor(reactor_name)"))
  model_B <- glm.nb(formula_B, data = table5_data)
  coeftest_B <- coeftest(model_B, vcov = vcovCL(model_B, ~ plantid)) 
  
  models_B[[var]] <- round(coeftest_B, 2)
  
  chisq_B <- linearHypothesis(model_B, "div_old = div_new", vcov = vcovCL(model_B, ~ plantid), test = "Chisq")
  chisq_results[[paste(var, "old_vs_new", sep = "_")]] <- list(
    chisq_stat = chisq_B$Chisq[2],
    p_value = chisq_B$`Pr(>Chisq)`[2]
  )
  
  formula_C <- as.formula(paste(var, "~ div_small + div_large + factor(year) + factor(reactor_name)"))
  model_C <- glm.nb(formula_C, data = table5_data)
  coeftest_C <- coeftest(model_C, vcov = vcovCL(model_C, ~ plantid)) 
  
  models_C[[var]] <- round(coeftest_C, 2)
  
  chisq_C <- linearHypothesis(model_C, "div_small = div_large", vcov = vcovCL(model_C, ~ plantid), test = "Chisq")
  chisq_results[[paste(var, "small_vs_large", sep = "_")]] <- list(
    chisq_stat = chisq_C$Chisq[2],
    p_value = chisq_C$`Pr(>Chisq)`[2]
  )
}

stargazer(models_A$initiating, models_A$fire, models_A$escalated, 
          type = "text", 
          title = "Model A: Divestiture BWR and divestiture PWR",
          covariate.labels = c("Divestiture, BWR", "Divestiture, PWR"),
          column.labels = c("Initiating events", "Fires", "Escalated enforcement"),
          omit = c("^factor\\(year\\)", "^factor\\(reactor_name\\)", "Constant"),
          omit.stat = c("f", "ser"))

stargazer(models_B$initiating, models_B$fire, models_B$escalated, 
          type = "text", 
          title = "Model B: Divestiture, old and Divestiture, new", 
          covariate.labels = c("Divestiture, old", "Divestiture, new"), 
          column.labels = c("Initiating events", "Fires", "Escalated enforcement"),
          omit = c("^factor\\(year\\)", "^factor\\(reactor_name\\)", "Constant"),
          omit.stat = c("f", "ser"))

stargazer(models_C$initiating, models_C$fire, models_C$escalated, 
          type = "text", 
          title = "Model C: Divestiture, small and Divestiture, large", 
          covariate.labels = c("Divestiture, small", "Divestiture, large"), 
          column.labels = c("Initiating events", "Fires", "Escalated enforcement"),
          omit = c("^factor\\(year\\)", "^factor\\(reactor_name\\)", "Constant"),
          omit.stat = c("f", "ser"))


print(chisq_results) # I show a cleaner output below

# The numbers of observations for each regression are exactly the same as those in Table 3,
  # since we collapsed the dataset in the same way. The same applies to the subsequent
  # regressions (for "personrem" and "rem").

# Continuous variables (columns 4 and 5 in Table 5)

table5bis_data <- data %>% filter(!is.na(rem))

table5bis_data <- table5bis_data %>%
  group_by(facilityname, year, plantid, pwr, year_divest, state_name, censusregion) %>%
  summarise(
    personrem = mean(personrem, na.rm = TRUE),
    rem = mean(rem, na.rm = TRUE),
    divested = mean(divested, na.rm = TRUE),
    capacity_factor = mean(capacity_factor, na.rm = TRUE),
    ever_divested = mean(ever_divested, na.rm = TRUE),
    capacity = min(capacity, na.rm = TRUE),
    year_start = min(year_start, na.rm = TRUE)
  )

table5bis_data <- table5bis_data %>%
  mutate(
    div_pwr = divested * pwr,
    div_bwr = divested * (1 - pwr),
    old = ifelse(year_start <= 1979, 1, 0),
    new = ifelse(year_start > 1979, 1, 0),
    div_new = divested * new,
    div_old = divested * old,
    small = ifelse(capacity < 1000, 1, 0),
    large = ifelse(capacity >= 1000, 1, 0),
    div_small = divested * small,
    div_large = divested * large
  )

# Again, for a matter of convenience I will call the three distinct models "D", "E" and "F".

models_D <- list()
models_E <- list()
models_F <- list()

chisq_results_bis <- list()

# Instead of using the panel data structure, we can as well perform directly lm()
  # and include the fixed effects by year and plant ID. (The dependent variables
  # are alrady clustered by group().) In any case, that's literally the same:
  # I checked the results with both methods and they are identical.
# As usual, for "personrem" and "rem" we cluster the s.e. by plant.

for (var in c("personrem", "rem")) {
  # Model D: y = div_brw + div_pwr + v_t + alpha_i (OLS)
  formula_D <- as.formula(paste(var, "~ div_bwr + div_pwr + factor(year) + factor(plantid)"))
  model_D <- lm(formula_D, data = table5bis_data)
  coeftest_D <- coeftest(model_D, vcov = vcovCL(model_D, ~ plantid))
  
  models_D[[var]] <- round(coeftest_D, 3)
  
  chisq_D <- linearHypothesis(model_D, "div_bwr = div_pwr", vcov = vcovCL(model_D, ~ plantid), test = "Chisq")
  chisq_results_bis[[paste(var, "bwr_vs_pwr", sep = "_")]] <- list(
    chisq_stat = chisq_D$Chisq[2],
    p_value = chisq_D$`Pr(>Chisq)`[2]
  )
  
  # Model E: y = div_old + div_large + v_t + alpha_i (OLS)
  formula_E <- as.formula(paste(var, "~ div_old + div_new + factor(year) + factor(plantid)"))
  model_E <- lm(formula_E, data = table5bis_data)
  coeftest_E <- coeftest(model_E, vcov = vcovCL(model_E, ~ plantid))
  
  models_E[[var]] <- round(coeftest_E, 3)
  
  chisq_E <- linearHypothesis(model_E, "div_old = div_new", vcov = vcovCL(model_E, ~ plantid), test = "Chisq")
  chisq_results_bis[[paste(var, "old_vs_new", sep = "_")]] <- list(
    chisq_stat = chisq_E$Chisq[2],
    p_value = chisq_E$`Pr(>Chisq)`[2]
  )
  
  # Model F: y = div_small + div_large + v_t + alpha_i (OLS)
  formula_F <- as.formula(paste(var, "~ div_small + div_large + factor(year) + factor(plantid)"))
  model_F <- lm(formula_F, data = table5bis_data)
  coeftest_F <- coeftest(model_F, vcov = vcovCL(model_F, ~ plantid))
  
  models_F[[var]] <- round(coeftest_F, 3)
  
  chisq_F <- linearHypothesis(model_F, "div_small = div_large", vcov = vcovCL(model_F, ~ plantid), test = "Chisq")
  chisq_results_bis[[paste(var, "small_vs_large", sep = "_")]] <- list(
    chisq_stat = chisq_F$Chisq[2],
    p_value = chisq_F$`Pr(>Chisq)`[2]
  )
}

stargazer(models_D$personrem, models_D$rem, 
          type = "text", 
          title = "Model D: Divestiture BWR and divestiture PWR",
          covariate.labels = c("Divestiture, BWR", "Divestiture, PWR"),
          column.labels = c("personrem", "rem"),
          omit = c("^factor\\(year\\)", "^factor\\(plantid\\)", "Constant"),
          omit.stat = c("f", "ser"))

stargazer(models_E$personrem, models_E$rem,
          type = "text", 
          title = "Model E: Divestiture, old and Divestiture, new", 
          covariate.labels = c("Divestiture, old", "Divestiture, new"), 
          column.labels = c("personrem", "rem"),
          omit = c("^factor\\(year\\)", "^factor\\(plantid\\)", "Constant"),
          omit.stat = c("f", "ser"))

stargazer(models_F$personrem, models_F$rem,
          type = "text", 
          title = "Model F: Divestiture, small and Divestiture, large", 
          covariate.labels = c("Divestiture, small", "Divestiture, large"), 
          column.labels = c("personrem", "rem"),
          omit = c("^factor\\(year\\)", "^factor\\(plantid\\)", "Constant"),
          omit.stat = c("f", "ser"))

# Chi-square statistics
chisq_table <- data.frame(
  variable = c("Initiating events", "Fires", "Escalated enforcement", "Radiation (person-rems)", "Average radiation (rems)"),
  chisq_bwr_vs_pwr = c(chisq_results$initiating_bwr_vs_pwr$chisq_stat, chisq_results$fire_bwr_vs_pwr$chisq_stat, chisq_results$escalated_bwr_vs_pwr$chisq_stat, chisq_results_bis$personrem_bwr_vs_pwr$chisq_stat, chisq_results_bis$rem_bwr_vs_pwr$chisq_stat),
  p_value1 = c(chisq_results$initiating_bwr_vs_pwr$p_value, chisq_results$fire_bwr_vs_pwr$p_value, chisq_results$escalated_bwr_vs_pwr$p_value, chisq_results_bis$personrem_bwr_vs_pwr$p_value, chisq_results_bis$rem_bwr_vs_pwr$p_value),
  chisq_old_vs_new = c(chisq_results$initiating_old_vs_new$chisq_stat, chisq_results$fire_old_vs_new$chisq_stat, chisq_results$escalated_old_vs_new$chisq_stat, chisq_results_bis$personrem_old_vs_new$chisq_stat, chisq_results_bis$rem_old_vs_new$chisq_stat),
  p_value2 = c(chisq_results$initiating_old_vs_new$p_value, chisq_results$fire_old_vs_new$p_value, chisq_results$escalated_old_vs_new$p_value, chisq_results_bis$personrem_old_vs_new$p_value, chisq_results_bis$rem_old_vs_new$p_value),
  chisq_small_vs_large = c(chisq_results$initiating_small_vs_large$chisq_stat, chisq_results$fire_small_vs_large$chisq_stat, chisq_results$escalated_small_vs_large$chisq_stat, chisq_results_bis$personrem_small_vs_large$chisq_stat, chisq_results_bis$rem_small_vs_large$chisq_stat),
  p_value3 = c(chisq_results$initiating_small_vs_large$p_value, chisq_results$fire_small_vs_large$p_value, chisq_results$escalated_small_vs_large$p_value, chisq_results_bis$personrem_small_vs_large$p_value, chisq_results_bis$rem_small_vs_large$p_value)
  )

chisq_table <- chisq_table %>% mutate(across(.cols = -variable, .fns = ~ round(.x, 3)))

print(chisq_table)

# To be fair, Hausman got one rejection region wrong in her table:
  # The p-value of the chi-square statistic for "div_bwr = div_pwr", on "rem",
  # is 0.009, which is lower than 0.01, so we reject the null at the 1% level,
  # not only at the 5% level.	
# ... the author probably just computed the p-value over the rounded value of 6.808,
  # that is 6.81.

