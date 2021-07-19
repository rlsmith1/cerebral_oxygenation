


# libraries ---------------------------------------------------------------


    library(tidyverse)
    library(tidymodels)
    library(vip)  
    library(glmnet)



# Train model -------------------------------------------------------------


    # create recipe
    my_rec <- recipe(Status ~ ., data = df_model) %>% 
      step_BoxCox(all_predictors()) %>%
      step_nzv(all_predictors())
    
    # model specification
    lasso_spec <- logistic_reg(penalty = 0.1, mixture = 1) %>% # picking specific penalty (tune in next section)
      set_engine("glmnet")
    
    # create workflow
    lasso_wf <- workflow() %>% 
      add_recipe(my_rec)
    
    # fit the workflow to model + data
    lasso_fit <- lasso_wf %>% 
      add_model(lasso_spec) %>% 
      fit(data = df_model)
    
    # extract results
    lasso_fit %>% 
      pull_workflow_fit() %>% 
      tidy()



    
    
# Tune LASSO parameters ---------------------------------------------------


    
    # resample
    set.seed(234)
    my_boot <- bootstraps(df_model, strata = group, times = 1000)
    
    # new specification
    tune_spec <- logistic_reg(penalty = tune(), mixture = 1) %>% # tune model to determine optimal penalty value
      set_engine("glmnet")
    
    # set grid
    lambda_grid <- grid_regular(penalty(),
                                levels = 50)
    
    # run across resamples
    doParallel::registerDoParallel()
    
    set.seed(2021)  
    lasso_grid <- tune_grid(
      lasso_wf %>% add_model(tune_spec),
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
      facet_wrap(~.metric, scales = "free", nrow = 2) +
      scale_x_log10() +
      theme(legend.position = "none")
    
    # select best penalty value
    highest_roc_auc <- lasso_grid %>% 
      select_best("roc_auc") %>% 
      select(penalty)

    # finalize wf
    final_lasso <- finalize_workflow(lasso_wf %>% add_model(tune_spec),
                                     highest_roc_auc)
    
    # fit model
    final_lasso_fit <- final_lasso %>% 
      fit(df_model) %>% 
      pull_workflow_fit()
    
    # extract results
    final_lasso_fit %>% tidy()
    
    
    
    
    
# Explore results with plots ----------------------------------------------


    
    # ROC curve
    lasso_grid %>% 
      collect_predictions(parameters = highest_roc_auc) %>% 
      roc_curve(truth = group, .pred_class1) %>% 
      
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
      count(group, .pred_class) %>% 
      group_by(group) %>% 
      mutate(total_true = sum(n), freq = n/total_true, `Percent occurence` = round(freq*100, 2)) %>% 
      
      ggplot(aes(x = group, y = .pred_class, fill = `Percent occurence`, label = `Percent occurence`)) +
      geom_tile() +
      geom_text(size = 6) +
      scale_fill_gradient2(low = "#4575b4", high = "#d73027") +
      theme_bw() 
    

    
    
    
    
    
    

