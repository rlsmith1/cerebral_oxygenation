

# model 3: THCa only -----------------------------------------------------------------


  # select only THC_alpha
  df_hbox_model3 <- df_hbox_voi_01 %>% 
    dplyr::select(c(Status, THC_alpha)) %>% 
    dplyr::rename("THCα" = "THC_alpha") %>% 
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
  
  
  

# compare models ----------------------------------------------------------

  collect_metrics(glm_res2)  # THCa + THC sd
  collect_metrics(glm_res3)  # THCa
  collect_metrics(glm_res4)  # THC + THCa
  
  
  tidy(hbox_final2)
  tidy(hbox_final3)
  tidy(hbox_final4)

  glm(Status ~ THC + THC_alpha, data = df_hbox_voi_01, family = "binomial") %>% summary()
  

                     