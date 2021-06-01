

# model 3: THCa only -----------------------------------------------------------------


  # select only THC_alpha
  df_hbox_model3 <- df_hbox_voi_01 %>% 
    dplyr::select(c(Status, THCα)) %>% 
    mutate(`1/10 * THCα` = THCα/10) %>% 
    dplyr::select(-THCα)
  
  # resample
  set.seed(8910)
  hbox_boot3 <- bootstraps(df_hbox_model3, times = 1000)
  
  # recipe
  hbox_rec3 <- recipe(Status ~ ., data = df_hbox_model3) %>% 
    step_BoxCox(all_predictors()) %>% 
    step_nzv(all_predictors())
  
  # workflow
  hbox_wf3 <- workflow() %>% 
    add_recipe(hbox_rec3)
  
  # logistic regression model spec
  glm_spec3 <- logistic_reg() %>% 
    set_engine("glm")
  
  glm_res3 <- hbox_wf3 %>% 
    add_model(glm_spec3) %>%
    fit_resamples(
      resamples = hbox_boot3,
      metrics = metric_set(roc_auc, accuracy, sensitivity, specificity),
      control = tune::control_resamples(save_pred = TRUE, verbose = TRUE)
    )
  
  # select best model to be used for evaluation
  glm_best3 <- glm_res3 %>% select_best("roc_auc")
  
  # finalize model
  hbox_final3 <- hbox_wf3 %>% 
    add_model(glm_spec3) %>% 
    finalize_workflow(glm_best3) %>%
    fit(df_hbox_model3) %>%
    pull_workflow_fit()
  
  ## Evaluate
  
  collect_metrics(glm_res3) %>% 
    dplyr::select(.metric, mean, std_err) %>% 
    mutate(mean = mean*100, std_err = std_err*100) %>% 
    mutate(Metric = c("Accuracy", "ROC AUC", "Sensitivity", "Specificity")) %>% 
    dplyr::select(c("Metric", "mean", "std_err")) %>% 
    
    dplyr::rename("Mean across resamples" = "mean", "Standard error" = "std_err") %>% 
    mutate_if(is.numeric, round, 2) 
  
  hbox_final3 %>% 
    tidy() %>% 
    mutate(p.value = as.character(p.value)) %>% 
    mutate_if(is.numeric, round, 2) %>% 
    mutate(p.value = as.numeric(p.value)) %>% 
    mutate(p.value = scientific(p.value, digits = 3)) %>% 
    dplyr::rename("Term" = "term", "Coefficient" = "estimate", "Standard error" = "std.error", 
                  "Statistic" = "statistic", "p-value" = "p.value")
    
  
  
  

# model 4: THC + THCa --------------------------------------------------------------


  # select  THC_alpha & THC
  df_hbox_model4 <- df_hbox_voi_01 %>% 
    dplyr::select(c(Status, THC, THC_alpha)) %>% 
    dplyr::rename("THCα" = "THC_alpha") %>% 
    mutate(`1/10 * THCα` = THCα/10) %>% 
    dplyr::select(-THCα)
  
  # resample
  set.seed(1011)
  hbox_boot4 <- bootstraps(df_hbox_model4, times = 100)
  
  # recipe
  hbox_rec4 <- recipe(Status ~ ., data = df_hbox_model4) %>% 
    step_BoxCox(all_predictors()) %>% 
    step_nzv(all_predictors())
  
  # workflow
  hbox_wf4 <- workflow() %>% 
    add_recipe(hbox_rec4)
  
  # logistic regression model spec
  glm_spec4 <- logistic_reg() %>% 
    set_engine("glm")
  
  glm_res4 <- hbox_wf4 %>% 
    add_model(glm_spec4) %>%
    fit_resamples(
      resamples = hbox_boot4,
      metrics = metric_set(roc_auc, accuracy, sensitivity, specificity),
      control = tune::control_resamples(save_pred = TRUE, verbose = TRUE)
    )
  
  # select best model to be used for evaluation
  glm_best4 <- glm_res4 %>% select_best("roc_auc")
  
  # finalize model
  hbox_final4 <- hbox_wf4 %>% 
    add_model(glm_spec4) %>% 
    finalize_workflow(glm_best4) %>%
    fit(df_hbox_model4) %>%
    pull_workflow_fit()
  
  ## Evaluate
  
  collect_metrics(glm_res4) %>% 
    dplyr::select(.metric, mean, std_err) %>% 
    mutate(mean = mean*100, std_err = std_err*100) %>% 
    mutate(Metric = c("Accuracy", "ROC AUC", "Sensitivity", "Specificity")) %>% 
    dplyr::select(c("Metric", "mean", "std_err")) %>% 
    
    dplyr::rename("Mean across resamples" = "mean", "Standard error" = "std_err") %>% 
    mutate_if(is.numeric, round, 2) # model performs worse
  
  hbox_final4 %>% 
    tidy() %>% 
    mutate(p.value = as.character(p.value)) %>% 
    mutate_if(is.numeric, round, 2) %>% 
    mutate(p.value = as.numeric(p.value)) %>% 
    mutate(p.value = scientific(p.value, digits = 3)) %>% 
    dplyr::rename("Term" = "term", "Coefficient" = "estimate", "Standard error" = "std.error", 
                  "Statistic" = "statistic", "p-value" = "p.value") # THC ns
  
  
  
# model 5: Cerebral O2 sat + cerebral O2 sat a + THC + THCa --------------------------------------------------------------
  
  
  # select  THC_alpha & THC
  `df_hbox_model5` <- df_hbox_voi_01 %>% 
    mutate(`1/10 * THCα` = THCα/10) %>% 
    dplyr::select(-THCα)
  
  # resample
  set.seed(1012)
  hbox_boot5 <- bootstraps(df_hbox_model5, times = 100)
  
  # recipe
  hbox_rec5 <- recipe(Status ~ ., data = df_hbox_model5) %>% 
    step_BoxCox(all_predictors()) %>% 
    step_nzv(all_predictors())
  
  # workflow
  hbox_wf5 <- workflow() %>% 
    add_recipe(hbox_rec5)
  
  # logistic regression model spec
  glm_spec5 <- logistic_reg() %>% 
    set_engine("glm")
  
  glm_res5 <- hbox_wf5 %>% 
    add_model(glm_spec5) %>%
    fit_resamples(
      resamples = hbox_boot5,
      metrics = metric_set(roc_auc, accuracy, sensitivity, specificity),
      control = tune::control_resamples(save_pred = TRUE, verbose = TRUE)
    )
  
  # select best model to be used for evaluation
  glm_best5 <- glm_res5 %>% select_best("roc_auc")
  
  # finalize model
  hbox_final5 <- hbox_wf5 %>% 
    add_model(glm_spec5) %>% 
    finalize_workflow(glm_best5) %>%
    fit(df_hbox_model5) %>%
    pull_workflow_fit()
  
  ## Evaluate
  
  collect_metrics(glm_res5) %>% 
    dplyr::select(.metric, mean, std_err) %>% 
    mutate(mean = mean*100, std_err = std_err*100) %>% 
    mutate(Metric = c("Accuracy", "ROC AUC", "Sensitivity", "Specificity")) %>% 
    dplyr::select(c("Metric", "mean", "std_err")) %>% 
    
    dplyr::rename("Mean across resamples" = "mean", "Standard error" = "std_err") %>% 
    mutate_if(is.numeric, round, 2) # model performs worse
  
  hbox_final5 %>% 
    tidy() %>% 
    mutate(p.value = as.character(p.value)) %>% 
    mutate_if(is.numeric, round, 2) %>% 
    mutate(p.value = as.numeric(p.value)) %>% 
    mutate(p.value = scientific(p.value, digits = 3)) %>% 
    dplyr::rename("Term" = "term", "Coefficient" = "estimate", "Standard error" = "std.error", 
                  "Statistic" = "statistic", "p-value" = "p.value") # THC ns
  
  
  
  
  
# model 6: THC + THCa + hematocrit + lactate --------------------------------------------------------------
  
  df_hbox_model6 <- df_hbox_impute %>% 
    dplyr::select(c(2,4:5,19:22)) %>% 
    mutate(`1/10 * THCα` = THCα/10) %>% 
    dplyr::select(-THCα) %>% 
    mutate(Status = as.factor(ifelse(Status == "CM", 1, 0))) %>% 
    
    dplyr::select(-c(`Cerebral O2 sat`, `Cerebral O2 sat α`))
  
  
  glm(Status ~ ., data = df_hbox_model6, family = "binomial") %>% summary()
  

  # resample
  set.seed(1013)
  hbox_boot6 <- bootstraps(df_hbox_model6, times = 100)
  
  # recipe
  hbox_rec6 <- recipe(Status ~ ., data = df_hbox_model6) %>% 
    step_BoxCox(all_predictors()) %>% 
    step_nzv(all_predictors())
  
  # workflow
  hbox_wf6 <- workflow() %>% 
    add_recipe(hbox_rec6)
  
  # logistic regression model spec
  glm_spec6 <- logistic_reg() %>% 
    set_engine("glm")
  
  glm_res6 <- hbox_wf6 %>% 
    add_model(glm_spec6) %>%
    fit_resamples(
      resamples = hbox_boot6,
      metrics = metric_set(roc_auc, accuracy, sensitivity, specificity),
      control = tune::control_resamples(save_pred = TRUE, verbose = TRUE)
    )
  
  # select best model to be used for evaluation
  glm_best6 <- glm_res6 %>% select_best("roc_auc")
  
  # finalize model
  hbox_final6 <- hbox_wf6 %>% 
    add_model(glm_spec6) %>% 
    finalize_workflow(glm_best6) %>%
    fit(df_hbox_model6) %>%
    pull_workflow_fit()
  
  ## Evaluate
  
  collect_metrics(glm_res6) %>% 
    dplyr::select(.metric, mean, std_err) %>% 
    mutate(mean = mean*100, std_err = std_err*100) %>% 
    mutate(Metric = c("Accuracy", "ROC AUC", "Sensitivity", "Specificity")) %>% 
    dplyr::select(c("Metric", "mean", "std_err")) %>% 
    
    dplyr::rename("Mean across resamples" = "mean", "Standard error" = "std_err") %>% 
    mutate_if(is.numeric, round, 2) # model performs worse
  
  hbox_final6 %>% 
    tidy() %>% 
    mutate(p.value = as.character(p.value)) %>% 
    mutate_if(is.numeric, round, 2) %>% 
    mutate(p.value = as.numeric(p.value)) %>% 
    mutate(p.value = scientific(p.value, digits = 3)) %>% 
    dplyr::rename("Term" = "term", "Coefficient" = "estimate", "Standard error" = "std.error", 
                  "Statistic" = "statistic", "p-value" = "p.value") # THC ns
  
  

  
# model 7: THC + THCa + hematocrit --------------------------------------------------------------
  
  df_hbox_model7 <- df_hbox_impute %>% 
    dplyr::select(c(2,4:5,19:22)) %>% 
    mutate(`1/10 * THCα` = THCα/10) %>% 
    dplyr::select(-THCα) %>% 
    mutate(Status = as.factor(ifelse(Status == "CM", 1, 0))) %>% 
    
    dplyr::select(-c(`Cerebral O2 sat`, `Cerebral O2 sat α`, `Lactate (mmol/L)`))
  
  
  
  # resample
  set.seed(1014)
  hbox_boot7 <- bootstraps(df_hbox_model7, times = 100)
  
  # recipe
  hbox_rec7 <- recipe(Status ~ ., data = df_hbox_model7) %>% 
    step_BoxCox(all_predictors()) %>% 
    step_nzv(all_predictors())
  
  # workflow
  hbox_wf7 <- workflow() %>% 
    add_recipe(hbox_rec7)
  
  # logistic regression model spec
  glm_spec7 <- logistic_reg() %>% 
    set_engine("glm")
  
  glm_res7 <- hbox_wf7 %>% 
    add_model(glm_spec7) %>%
    fit_resamples(
      resamples = hbox_boot7,
      metrics = metric_set(roc_auc, accuracy, sensitivity, specificity),
      control = tune::control_resamples(save_pred = TRUE, verbose = TRUE)
    )
  
  # select best model to be used for evaluation
  glm_best7 <- glm_res7 %>% select_best("roc_auc")
  
  # finalize model
  hbox_final7 <- hbox_wf7 %>% 
    add_model(glm_spec7) %>% 
    finalize_workflow(glm_best7) %>%
    fit(df_hbox_model7) %>%
    pull_workflow_fit()
  
  ## Evaluate
  
  collect_metrics(glm_res7) %>% 
    dplyr::select(.metric, mean, std_err) %>% 
    mutate(mean = mean*100, std_err = std_err*100) %>% 
    mutate(Metric = c("Accuracy", "ROC AUC", "Sensitivity", "Specificity")) %>% 
    dplyr::select(c("Metric", "mean", "std_err")) %>% 
    
    dplyr::rename("Mean across resamples" = "mean", "Standard error" = "std_err") %>% 
    mutate_if(is.numeric, round, 2) # model performs worse
  
  hbox_final7 %>% 
    tidy() %>% 
    mutate(p.value = as.character(p.value)) %>% 
    mutate_if(is.numeric, round, 2) %>% 
    mutate(p.value = as.numeric(p.value)) %>% 
    mutate(p.value = scientific(p.value, digits = 3)) %>% 
    dplyr::rename("Term" = "term", "Coefficient" = "estimate", "Standard error" = "std.error", 
                  "Statistic" = "statistic", "p-value" = "p.value") 
  
  
  
  
  
  

# compare models ----------------------------------------------------------

  collect_metrics(glm_res2)  # THCa + THC sd
  collect_metrics(glm_res3)  # THCa
  collect_metrics(glm_res4)  # THC + THCa
  collect_metrics(glm_res5)  # THC + THCa + Ceb O2 sat + Ceb O2 sat a
  collect_metrics(glm_res6)  # THC + THCa + hematocrit + lactate
  collect_metrics(glm_res7)  # THC + THCa + hematocrit
  
  
  tidy(hbox_final2)
  tidy(hbox_final3) #
  tidy(hbox_final4)
  tidy(hbox_final5)
  tidy(hbox_final6) #
  tidy(hbox_final7) #
  

  glm(Status ~ ., data = df_hbox_voi_01, family = "binomial") %>% summary()
  glm(Status ~ THC + THCα + `Cerebral O2 sat`, data = df_hbox_voi_01, family = "binomial") %>% summary()
  glm(Status ~ THC + THCα, data = df_hbox_voi_01, family = "binomial") %>% summary()
  glm(Status ~ THCα, data = df_hbox_voi_01, family = "binomial") %>% summary()
  