

### not used in final manuscript because we want to be able to evaluate coefficients


# libraries ---------------------------------------------------------------


    require(readxl)
    library(tidyverse)
    library(tidymodels)
    library(themis)
    library(vip)  
    library(glmnet)



# data --------------------------------------------------------------------


  df_model <- read.csv("Data/master_datatable.csv") %>% 
    as_tibble() %>% 
    select(c(subject_id, Status, Hematocrit, Lactate,
             cerebral_hb_tot, cerebral_hbtot_long_a, cerebral_hb_oxy, cerebral_hboxy_long_a))


# Build model -------------------------------------------------------------

   
  # Remove HC and look at missing data
  df_model <- df_model %>% filter(Status != "HC") %>% mutate(Status = factor(Status, levels = c("CM", "UM")))
  df_model %>% is.na %>% colSums
  df_model %>% filter(is.na(cerebral_hboxy_long_a))
  
  require(naniar)
  df_model %>% select(-c(subject_id, Status)) %>% gg_miss_upset() # one patient is missing Hct and Lactate, 3 patients missing Hb_oxy alpha
  df_model %>% count(Status) 
  
  
    # create recipe
    my_rec <- recipe(Status ~ ., data = df_model) %>% 
      step_impute_knn(all_predictors()) %>% 
      step_normalize(all_predictors()) %>%
      step_zv(all_predictors()) %>%
      step_smote(Status) %>% 
      update_role(subject_id, new_role = "id")
    
    # prep recipe
    my_prep <- my_rec %>% prep(strings_as_factors = FALSE)
    
    # model specification
    lasso_spec <- logistic_reg(penalty = tune(), mixture = 1) %>%
      set_engine("glmnet")
    
    # create workflow
    lasso_wf <- workflow() %>% 
      add_recipe(my_rec) %>% 
      add_model(lasso_spec)
    

        
# Tune LASSO parameters ---------------------------------------------------


    # set grid
    lambda_grid <- grid_regular(penalty(), levels = 50)
    
    # resample
    set.seed(234)
    my_boot <- bootstraps(df_model, strata = Status)
    
    # run across resamples
    doParallel::registerDoParallel()
    
    set.seed(2021)  
    lasso_grid <- tune_grid(
      lasso_wf,
      resamples = my_boot,
      metrics = metric_set(roc_auc, accuracy, sensitivity, specificity),
      control = control_grid(save_pred = TRUE),
      grid = lambda_grid
    )  
    
    # look at metrics
    lasso_grid %>% 
      collect_metrics() %>% 
      
      ggplot(aes(penalty, mean, color = .metric)) +
      geom_errorbar(aes(ymin = mean - std_err,
                        ymax = mean + std_err)) +
      geom_line(size = 1.5, show.legend = FALSE) +
      facet_wrap(~.metric, nrow = 2) +
      scale_x_log10() +
      theme(legend.position = "none")
  

    
# choose the final model --------------------------------------------------


      
    # select best penalty value
    highest_roc_auc <- lasso_grid %>% 
      select_best("roc_auc") %>% 
      select(penalty)
    
    # show metrics
    lasso_grid %>% 
      collect_metrics() %>% 
      filter(penalty == highest_roc_auc$penalty)

    # finalize wf
    final_lasso <- finalize_workflow(lasso_wf, highest_roc_auc)
    
    # fit model
    final_lasso_fit <- final_lasso %>% 
      fit(df_model) %>% 
      extract_fit_parsnip()
    
    # extract results
    final_lasso_fit %>% tidy()
    
    
    
    
# Explore results with plots ----------------------------------------------


    # ROC curve
    lasso_grid %>% 
      collect_predictions(parameters = highest_roc_auc) %>% 
      roc_curve(truth = Status, .pred_class) %>% 
      
      ggplot(aes(x = 1 - specificity, y = sensitivity)) +
      geom_line() +
      geom_abline(intercept = 0, slope = 1, lty = 2)
    
    # variable importance bar chart
    final_lasso_fit %>% 
      vi(lambda = highest_roc_auc$penalty) %>% 
      mutate(Variable = fct_reorder(Variable, Importance), 
             Importance = ifelse(Sign == "NEG", Importance*-1, Importance),
             Sign = factor(Sign, levels = c("POS", "NEG"))) %>% 
      
      ggplot(aes(x = Importance, y = Variable, fill = Sign)) +
      geom_col() +
      theme_bw()
    
    # confusion matrix
    lasso_grid %>% 
      collect_predictions(parameters = highest_roc_auc) %>% 
      count(Status, .pred_class) %>% 
      group_by(Status) %>% 
      mutate(total_true = sum(n), freq = n/total_true, `Percent occurence` = round(freq*100, 2)) %>% 
      
      ggplot(aes(x = Status, y = .pred_class, fill = `Percent occurence`, label = `Percent occurence`)) +
      geom_tile() +
      geom_text(size = 6) +
      scale_fill_gradient2(low = "#4575b4", high = "#d73027") +
      theme_bw() 
    

    
    
    


