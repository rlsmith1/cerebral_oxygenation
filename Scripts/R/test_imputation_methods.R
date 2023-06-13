


# libraries ---------------------------------------------------------------

library(tidyverse)
library(RColorBrewer)
library(gridExtra)
library(patchwork)
library(Hmisc)
library(corrr)
library(tidymodels)
library(themis)
library(vip)
library(janitor)
library(kableExtra)



# data --------------------------------------------------------------------


df_hbox <- read_csv("Data/Hans_datatable_exports/malawi key data v04Jan2022.csv") %>% 
  as_tibble() %>% 
  janitor::clean_names() %>% 
  mutate(status = replace(status, status == "HV", "HC")) %>% 
  dplyr::rename("subject_id" = "subject_id_session") %>% 
  
  dplyr::filter(status %in% c("CM", "HC", "UM") & session_no == 1 | subject_id == "TM0003CM01") %>% 
  dplyr::filter(!grepl("blood", subject_id)) %>% 
  mutate(status = factor(status, levels = c("HC", "UM", "CM"))) %>% 
  mutate_at(c("bp_dias", "bp_sys", "glucose", "lactate"), as.numeric) %>% # converts remaining character variables to numerics
  select(subject_id, status, everything()) %>% 
  group_by(status)



# data format -------------------------------------------------------------


df_model <- df_hbox %>% 
  mutate(sex = ifelse(sex == "Female", 1, 0)) %>% 
  filter(status != "HC") %>% 
  mutate(status = factor(status, levels = c("UM", "CM"))) %>% 
  
  # transform non-normal variables (by Shapiro-Wilk) using log-transform
  mutate(log10_cerebral_hb_tot = log10(cerebral_hb_tot),
         log10_resp_rate = log10(resp_rate),
         log10_lactate = log10(lactate)) %>% 
  
  # select variables of interest
  dplyr::select(subject_id, status, log10_cerebral_hb_tot, cerebral_hb_tot_alpha2, hct, log10_lactate, sex)



# impute with median ------------------------------------------------------

set.seed(20221125)

# resample
hbtot_boot <- bootstraps(df_model, strata = status, times = 5000)

### Hb_tot

# recipe
hbtot_rec <- recipe(status ~ log10_cerebral_hb_tot + hct + log10_lactate + sex + subject_id, data = df_model) %>% 
  step_impute_median(all_predictors()) %>% 
  step_normalize(all_predictors()) %>%
  step_zv(all_predictors()) %>%
  step_smote(status) %>% 
  update_role(subject_id, new_role = "id")

# create workflow
hbtot_wf <- workflow() %>% add_recipe(hbtot_rec)

# logistic regression model specifications
hbtot_glm_spec <- logistic_reg() %>% set_engine("glm")

# create model on 5000 different resamples
hbtot_glm_res <- hbtot_wf %>% 
  add_model(hbtot_glm_spec) %>%
  fit_resamples(
    resamples = hbtot_boot,
    metrics = metric_set(roc_auc, accuracy, sensitivity, specificity),
    control = tune::control_resamples(save_pred = TRUE, verbose = TRUE)
  )

# select best model to be used for evaluation
hbtot_glm_best <- hbtot_glm_res %>% select_best("roc_auc")

# finalize model
hbtot_final <- hbtot_wf %>% 
  add_model(hbtot_glm_spec) %>% 
  finalize_workflow(hbtot_glm_best) %>%
  fit(df_model) %>%
  extract_fit_parsnip()

# Odds ratios
df_CI <- confint(hbtot_final$fit) %>% exp() %>% as_tibble(rownames = "term")

hbtot_final %>%
  tidy(exponentiate = TRUE) %>%
  #filter(p.value < 0.05) %>% 
  left_join(df_CI, by = "term") %>%
  mutate(estimate = ifelse(estimate < 10^-2, scientific(estimate, digits = 3), round(estimate, 2)),
         conf.low = ifelse(`2.5 %` < 10^-2, scientific(`2.5 %`, digits = 3), round(`2.5 %`, 2)),
         conf.high = ifelse(`97.5 %` < 10^-2, scientific(`97.5 %`, digits = 3), round(`97.5 %`, 2))) %>%
  mutate(p.value = scientific(p.value, digits = 3),
         statistic = round(statistic, 2)) %>%
  unite("95% CI", c(conf.low, conf.high), sep = ", ") %>%
  dplyr::select(c(term, estimate, `95% CI`, statistic, p.value)) %>% 
  dplyr::rename("Term" = "term", "Odds ratio" = "estimate",
                "t-statistic" = "statistic", "p-value" = "p.value") %>%
  filter(Term != "(Intercept)") %>% 
  
  write.csv("Outputs/test_imputation_methods/cerebral_Hb_tot_median.csv", row.names = FALSE)
  


### Hb_tot alpha

# recipe
hbtot_alpha_rec <- recipe(status ~ cerebral_hb_tot_alpha2 + hct + log10_lactate + sex + subject_id, data = df_model) %>% 
  step_impute_median(all_predictors()) %>% 
  step_normalize(all_predictors()) %>%
  step_zv(all_predictors()) %>%
  step_smote(status) %>% 
  update_role(subject_id, new_role = "id")

# create workflow
hbtot_alpha_wf <- workflow() %>% add_recipe(hbtot_alpha_rec)

# logistic regression model specifications
hbtot_alpha_glm_spec <- logistic_reg() %>% set_engine("glm")

# create model on 5000 different resamples
hbtot_alpha_glm_res <- hbtot_alpha_wf %>% 
  add_model(hbtot_alpha_glm_spec) %>%
  fit_resamples(
    resamples = hbtot_boot,
    metrics = metric_set(roc_auc, accuracy, sensitivity, specificity),
    control = tune::control_resamples(save_pred = TRUE, verbose = TRUE)
  )

# select best model to be used for evaluation
hbtot_alpha_glm_best <- hbtot_alpha_glm_res %>% select_best("roc_auc")

# finalize model
hbtot_alpha_final <- hbtot_alpha_wf %>% 
  add_model(hbtot_alpha_glm_spec) %>% 
  finalize_workflow(hbtot_alpha_glm_best) %>%
  fit(df_model) %>%
  extract_fit_parsnip()

# Odds ratios
df_CI <- confint(hbtot_alpha_final$fit) %>% exp() %>% as_tibble(rownames = "term")

hbtot_alpha_final %>%
  tidy(exponentiate = TRUE) %>%
  #filter(p.value < 0.05) %>% 
  left_join(df_CI, by = "term") %>%
  mutate(estimate = ifelse(estimate < 10^-2, scientific(estimate, digits = 3), round(estimate, 2)),
         conf.low = ifelse(`2.5 %` < 10^-2, scientific(`2.5 %`, digits = 3), round(`2.5 %`, 2)),
         conf.high = ifelse(`97.5 %` < 10^-2, scientific(`97.5 %`, digits = 3), round(`97.5 %`, 2))) %>%
  mutate(p.value = scientific(p.value, digits = 3),
         statistic = round(statistic, 2)) %>%
  unite("95% CI", c(conf.low, conf.high), sep = ", ") %>%
  dplyr::select(c(term, estimate, `95% CI`, statistic, p.value)) %>% 
  dplyr::rename("Term" = "term", "Odds ratio" = "estimate",
                "t-statistic" = "statistic", "p-value" = "p.value") %>%
  filter(Term != "(Intercept)") %>% 
  
  write.csv("Outputs/test_imputation_methods/cerebral_Hb_tot_alpha_median.csv", row.names = FALSE)



# impute with KNN ------------------------------------------------------


### Hb_tot

# recipe
hbtot_rec <- recipe(status ~ log10_cerebral_hb_tot + hct + log10_lactate + sex + subject_id, data = df_model) %>% 
  step_impute_knn(all_predictors()) %>% 
  step_normalize(all_predictors()) %>%
  step_zv(all_predictors()) %>%
  step_smote(status) %>% 
  update_role(subject_id, new_role = "id")

# create workflow
hbtot_wf <- workflow() %>% add_recipe(hbtot_rec)

# logistic regression model specifications
hbtot_glm_spec <- logistic_reg() %>% set_engine("glm")

# create model on 5000 different resamples
hbtot_glm_res <- hbtot_wf %>% 
  add_model(hbtot_glm_spec) %>%
  fit_resamples(
    resamples = hbtot_boot,
    metrics = metric_set(roc_auc, accuracy, sensitivity, specificity),
    control = tune::control_resamples(save_pred = TRUE, verbose = TRUE)
  )

# select best model to be used for evaluation
hbtot_glm_best <- hbtot_glm_res %>% select_best("roc_auc")

# finalize model
hbtot_final <- hbtot_wf %>% 
  add_model(hbtot_glm_spec) %>% 
  finalize_workflow(hbtot_glm_best) %>%
  fit(df_model) %>%
  extract_fit_parsnip()

# Odds ratios
df_CI <- confint(hbtot_final$fit) %>% exp() %>% as_tibble(rownames = "term")

hbtot_final %>%
  tidy(exponentiate = TRUE) %>%
  #filter(p.value < 0.05) %>% 
  left_join(df_CI, by = "term") %>%
  mutate(estimate = ifelse(estimate < 10^-2, scientific(estimate, digits = 3), round(estimate, 2)),
         conf.low = ifelse(`2.5 %` < 10^-2, scientific(`2.5 %`, digits = 3), round(`2.5 %`, 2)),
         conf.high = ifelse(`97.5 %` < 10^-2, scientific(`97.5 %`, digits = 3), round(`97.5 %`, 2))) %>%
  mutate(p.value = scientific(p.value, digits = 3),
         statistic = round(statistic, 2)) %>%
  unite("95% CI", c(conf.low, conf.high), sep = ", ") %>%
  dplyr::select(c(term, estimate, `95% CI`, statistic, p.value)) %>% 
  dplyr::rename("Term" = "term", "Odds ratio" = "estimate",
                "t-statistic" = "statistic", "p-value" = "p.value") %>%
  filter(Term != "(Intercept)") %>% 
  
  write.csv("Outputs/test_imputation_methods/cerebral_Hb_tot_knn.csv", row.names = FALSE)



### Hb_tot alpha


# recipe
hbtot_alpha_rec <- recipe(status ~ cerebral_hb_tot_alpha2 + hct + log10_lactate + sex + subject_id, data = df_model) %>% 
  step_impute_knn(all_predictors()) %>% 
  step_normalize(all_predictors()) %>%
  step_zv(all_predictors()) %>%
  step_smote(status) %>% 
  update_role(subject_id, new_role = "id")

# create workflow
hbtot_alpha_wf <- workflow() %>% add_recipe(hbtot_alpha_rec)

# logistic regression model specifications
hbtot_alpha_glm_spec <- logistic_reg() %>% set_engine("glm")

# create model on 5000 different resamples
hbtot_alpha_glm_res <- hbtot_alpha_wf %>% 
  add_model(hbtot_alpha_glm_spec) %>%
  fit_resamples(
    resamples = hbtot_boot,
    metrics = metric_set(roc_auc, accuracy, sensitivity, specificity),
    control = tune::control_resamples(save_pred = TRUE, verbose = TRUE)
  )

# select best model to be used for evaluation
hbtot_alpha_glm_best <- hbtot_alpha_glm_res %>% select_best("roc_auc")

# finalize model
hbtot_alpha_final <- hbtot_alpha_wf %>% 
  add_model(hbtot_alpha_glm_spec) %>% 
  finalize_workflow(hbtot_alpha_glm_best) %>%
  fit(df_model) %>%
  extract_fit_parsnip()

# Odds ratios
df_CI <- confint(hbtot_alpha_final$fit) %>% exp() %>% as_tibble(rownames = "term")

hbtot_alpha_final %>%
  tidy(exponentiate = TRUE) %>%
  #filter(p.value < 0.05) %>% 
  left_join(df_CI, by = "term") %>%
  mutate(estimate = ifelse(estimate < 10^-2, scientific(estimate, digits = 3), round(estimate, 2)),
         conf.low = ifelse(`2.5 %` < 10^-2, scientific(`2.5 %`, digits = 3), round(`2.5 %`, 2)),
         conf.high = ifelse(`97.5 %` < 10^-2, scientific(`97.5 %`, digits = 3), round(`97.5 %`, 2))) %>%
  mutate(p.value = scientific(p.value, digits = 3),
         statistic = round(statistic, 2)) %>%
  unite("95% CI", c(conf.low, conf.high), sep = ", ") %>%
  dplyr::select(c(term, estimate, `95% CI`, statistic, p.value)) %>% 
  dplyr::rename("Term" = "term", "Odds ratio" = "estimate",
                "t-statistic" = "statistic", "p-value" = "p.value") %>%
  filter(Term != "(Intercept)") %>% 
  
  write.csv("Outputs/test_imputation_methods/cerebral_Hb_tot_alpha_knn.csv", row.names = FALSE)



# remove NA values --------------------------------------------------------



### Hb_tot

# recipe
hbtot_rec <- recipe(status ~ log10_cerebral_hb_tot + hct + log10_lactate + sex + subject_id, data = df_model) %>% 
  step_naomit(all_predictors()) %>% 
  step_normalize(all_predictors()) %>%
  step_zv(all_predictors()) %>%
  step_smote(status) %>% 
  update_role(subject_id, new_role = "id")

# create workflow
hbtot_wf <- workflow() %>% add_recipe(hbtot_rec)

# logistic regression model specifications
hbtot_glm_spec <- logistic_reg() %>% set_engine("glm")

# create model on 5000 different resamples
hbtot_glm_res <- hbtot_wf %>% 
  add_model(hbtot_glm_spec) %>%
  fit_resamples(
    resamples = hbtot_boot,
    metrics = metric_set(roc_auc, accuracy, sensitivity, specificity),
    control = tune::control_resamples(save_pred = TRUE, verbose = TRUE)
  )

# select best model to be used for evaluation
hbtot_glm_best <- hbtot_glm_res %>% select_best("roc_auc")

# finalize model
hbtot_final <- hbtot_wf %>% 
  add_model(hbtot_glm_spec) %>% 
  finalize_workflow(hbtot_glm_best) %>%
  fit(df_model) %>%
  extract_fit_parsnip()

hbtot_final %>% 
  tidy %>% 
  filter(term != "(Intercept)") %>% 
  mutate_if(is.numeric, ~ round(.x, 3))  # AIC 81.92, every term is significant

# Odds ratios
df_CI <- confint(hbtot_final$fit) %>% exp() %>% as_tibble(rownames = "term")

hbtot_final %>%
  tidy(exponentiate = TRUE) %>%
  #filter(p.value < 0.05) %>% 
  left_join(df_CI, by = "term") %>%
  mutate(estimate = ifelse(estimate < 10^-2, scientific(estimate, digits = 3), round(estimate, 2)),
         conf.low = ifelse(`2.5 %` < 10^-2, scientific(`2.5 %`, digits = 3), round(`2.5 %`, 2)),
         conf.high = ifelse(`97.5 %` < 10^-2, scientific(`97.5 %`, digits = 3), round(`97.5 %`, 2))) %>%
  mutate(p.value = scientific(p.value, digits = 3),
         statistic = round(statistic, 2)) %>%
  unite("95% CI", c(conf.low, conf.high), sep = ", ") %>%
  dplyr::select(c(term, estimate, `95% CI`, statistic, p.value)) %>% 
  dplyr::rename("Term" = "term", "Odds ratio" = "estimate",
                "t-statistic" = "statistic", "p-value" = "p.value") %>%
  filter(Term != "(Intercept)") %>% 
  
  write.csv("Outputs/test_imputation_methods/cerebral_Hb_tot_naomit.csv", row.names = FALSE)



### Hb_tot alpha


# recipe
hbtot_alpha_rec <- recipe(status ~ cerebral_hb_tot_alpha2 + hct + log10_lactate + sex + subject_id, data = df_model) %>% 
  step_naomit(all_predictors()) %>% 
  step_normalize(all_predictors()) %>%
  step_zv(all_predictors()) %>%
  step_smote(status) %>% 
  update_role(subject_id, new_role = "id")

# create workflow
hbtot_alpha_wf <- workflow() %>% add_recipe(hbtot_alpha_rec)

# logistic regression model specifications
hbtot_alpha_glm_spec <- logistic_reg() %>% set_engine("glm")

# create model on 5000 different resamples
hbtot_alpha_glm_res <- hbtot_alpha_wf %>% 
  add_model(hbtot_alpha_glm_spec) %>%
  fit_resamples(
    resamples = hbtot_boot,
    metrics = metric_set(roc_auc, accuracy, sensitivity, specificity),
    control = tune::control_resamples(save_pred = TRUE, verbose = TRUE)
  )

# select best model to be used for evaluation
hbtot_alpha_glm_best <- hbtot_alpha_glm_res %>% select_best("roc_auc")

# finalize model
hbtot_alpha_final <- hbtot_alpha_wf %>% 
  add_model(hbtot_alpha_glm_spec) %>% 
  finalize_workflow(hbtot_alpha_glm_best) %>%
  fit(df_model) %>%
  extract_fit_parsnip()

# Odds ratios
df_CI <- confint(hbtot_alpha_final$fit) %>% exp() %>% as_tibble(rownames = "term")

hbtot_alpha_final %>%
  tidy(exponentiate = TRUE) %>%
  #filter(p.value < 0.05) %>% 
  left_join(df_CI, by = "term") %>%
  mutate(estimate = ifelse(estimate < 10^-2, scientific(estimate, digits = 3), round(estimate, 2)),
         conf.low = ifelse(`2.5 %` < 10^-2, scientific(`2.5 %`, digits = 3), round(`2.5 %`, 2)),
         conf.high = ifelse(`97.5 %` < 10^-2, scientific(`97.5 %`, digits = 3), round(`97.5 %`, 2))) %>%
  mutate(p.value = scientific(p.value, digits = 3),
         statistic = round(statistic, 2)) %>%
  unite("95% CI", c(conf.low, conf.high), sep = ", ") %>%
  dplyr::select(c(term, estimate, `95% CI`, statistic, p.value)) %>% 
  dplyr::rename("Term" = "term", "Odds ratio" = "estimate",
                "t-statistic" = "statistic", "p-value" = "p.value") %>%
  filter(Term != "(Intercept)") %>% 
  
  write.csv("Outputs/test_imputation_methods/cerebral_Hb_tot_alpha_naomit.csv", row.names = FALSE)


