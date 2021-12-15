


# Libraries ---------------------------------------------------------------

        require(readxl)
        library(tidyverse)
        library(nortest)
        library(FSA)
        library(tidymodels)
        library(themis)
        library(factoextra)
        library(plotly)
        library(usemodels)
        library(ranger)
        library(vip)



# Data --------------------------------------------------------------------

        # read in data
        df_master <- read_xlsx("Data/master_datatable.csv") %>% as_tibble() %>% select(-X)
        

# Format data -------------------------------------------------------------
        
        
        # Trim for variables of interest (VOI) - variables that are significantly different between UM and CM
        df_master_voi <- df_master %>% 
                dplyr::select(c(subject_id, Status,
                                cerebral_hb_oxy, cerebral_hboxy_long_a, cerebral_hb_tot, cerebral_hbtot_long_a,
                                Hematocrit, Lactate)) 
        
        # impute missing data with median by group
        f_impute <- function(x, na.rm = TRUE) (replace(x, is.na(x), median(x, na.rm = na.rm)))
        
        df_master_impute <- df_master_voi %>% 
          group_by(Status) %>% 
          mutate_at(3:ncol(.), f_impute) %>% 
          ungroup() %>% 
          filter(Status != "HC")
        
        # pivot long
        df_master_voi_long <- df_master_voi %>%
          filter(Status != "HC") %>% 
          pivot_longer(cols = 3:ncol(.), names_to = "variable", values_to = "value")
        


        
# Explore data with plots -------------------------------------------------------

        
        # histogram to explore normality of data
        df_master_voi_long %>% 
                
                ggplot(aes(x = value)) +
                geom_histogram(bins = 20) +
                facet_grid(Status ~ variable, scales = "free")
        
        # plot all
        df_master_voi_long %>% 
                
                ggplot(aes(x = Status, y = value)) +
                geom_violin(aes(color = Status)) +
                geom_jitter(position = position_jitter(0.2), shape = 1) +
                stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
                facet_wrap(~variable, scales = "free_y") +
                
                theme_bw() +
                #ggtitle("cerebral.hbdeox") +
                theme(legend.position = "none")
   
        
        
        

# Use PCA to explore if these variables cluster Status --------------------

        # convert tibble to df
        df_master_impute_df <- df_master_impute %>% 
          select(-c(subject_id, Status)) %>% 
          as.data.frame(df_master_impute) 
        rownames(df_master_impute_df) <- df_master_impute$subject_id
        
        # compute PCA
        my_pca <- prcomp(df_master_impute_df, scale = TRUE)
        
        # Scree plot (visualize eigenvalues)
        fviz_eig(my_pca)
        
        # Graph of individual patients
        fviz_pca_ind(my_pca,
                     col.ind = "cos2", # Color by the quality of representation
                     gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                     repel = TRUE     # Avoid text overlapping
        )
        
        # Graph of variables
        fviz_pca_var(my_pca,
                     col.var = "contrib", # Color by contributions to the PC
                     gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                     repel = TRUE     # Avoid text overlapping
        )
        
        # Biplot of individuals and variables
        groups <- as.factor(df_master_impute$Status)
        p_biplot <- fviz_pca_biplot(my_pca, repel = TRUE,
                                    col.var = "black", # Variables color
                                    col.ind = groups  # Individuals color
        )
        
        p_biplot$layers[[3]] <- NULL
        p_biplot$layers[[2]] <- NULL
        p_biplot$layers[[1]]$geom_params$size <- 3 # didn't work
        
        # Eigenvalues
        pca_eig_val <- get_eigenvalue(my_pca)
        
        # Results for Variables
        pca_res_var <- get_pca_var(my_pca)
        pca_res_var$coord          # Coordinates
        pca_res_var$contrib        # Contributions to the PCs
        pca_res_var$cos2           # Quality of representation 
        
        # Results for individuals
        pca_res_ind <- get_pca_ind(my_pca)
        pca_res_ind$coord          # Coordinates
        pca_res_ind$contrib        # Contributions to the PCs
        pca_res_ind$cos2           # Quality of representation 
        
        # 2D PCA
        p_pca <- fviz_pca_ind(my_pca,
                              col.ind = groups, # color by groups
                              palette = c("#00AFBB",  "#FC4E07", "#7CAE00"),
                              addEllipses = TRUE, # Concentration ellipses
                              ellipse.type = "confidence",
                              legend.title = "Groups",
                              repel = TRUE
        )
        
        p_pca$layers[[4]] <- NULL # get rid of text labels
        
        # 3D PCA
        df_pca <- pca_res_ind$coord %>% 
                as_tibble(rownames = "subject") %>% 
                mutate(Status = df_hbox_impute$Status)
        
        plot_ly(x = df_pca$Dim.1, y = df_pca$Dim.2, z = df_pca$Dim.3, 
                type = "scatter3d", mode = "markers", color = df_pca$Status)
        
        # biggest contributors
        df_pca_contrib <- pca_res_var$contrib %>% 
                as_tibble(rownames = "variable") %>% 
                arrange(-Dim.1, -Dim.2)
        
        df_contrib <- tibble(Dim.1 = paste0(df_pca_contrib$variable, 
                                            sep = ": ", 
                                            round(df_pca_contrib$Dim.1, 2)),
                             Dim.2 = paste0(arrange(df_pca_contrib, -Dim.2)$variable, 
                                            sep = ": ", 
                                            round(arrange(df_pca_contrib, -Dim.2)$Dim.2, 2)),
                             Dim.3 = paste0(arrange(df_pca_contrib, -Dim.3)$variable, 
                                            sep = ": ", 
                                            round(arrange(df_pca_contrib, -Dim.3)$Dim.3, 2)))
        






# variable correlations ---------------------------------------------------

        f_calc_cor <- function(df, status) {
          
          df %>%
            dplyr::filter(Status == status) %>% 
            dplyr::select(-c(subject_id, Status)) %>% 
            cor(method = "spearman")
          
        }
        
        # create run correlations of potential regressors
        m_cor_hv <- df_master_impute %>% f_calc_cor(status = "HC")
        
        m_cor_um <- df_master_impute %>% f_calc_cor(status = "UM")
        
        m_cor_cm <- df_master_impute %>% f_calc_cor(status = "CM")
        
        # correlation matrix
        
        # write plot function
        f_plot_cor <- function(mat) {
          
          reshape2::melt(mat) %>% 
            
            mutate_if(is.factor, factor, levels = c("Lactate", "Hematocrit", 
                                                    "cerebral_hb_oxy", "cerebral_hb_tot", 
                                                    "cerebral_hboxy_long_a", "cerebral_hbtot_long_a")) %>% 
            
            ggplot(aes(x = Var1, y = Var2, fill = value, label = round(value, 2))) + 
            geom_tile() +
            scale_fill_gradient2(low = "#4575b4", high = "#d73027") +
            geom_text() +
            
            xlab("") +
            ylab("") +
            theme_bw()
            
          
        }
        
        # extract legend
        p_leg <- m_cor_hv %>% f_plot_cor() %>% get_legend()
        
        # create plots
        m_cor_um %>% f_plot_cor
        m_cor_cm %>% f_plot_cor

        
# Logistic regression model --------------------------------------------------------------


        # Data
        df_model <- df_master_voi %>% select(-subject_id)
        
        ## write function to create model on given df
        
        f_log_reg_model <- function(df_model, times) {
          
          # resample
          hbox_boot <- bootstraps(df_model, strata = Status, times = times)
          
          # recipe
          hbox_rec <- recipe(Status ~ ., data = df_model) %>% 
            step_BoxCox(all_predictors()) %>% # Box-Cox transformation of data
            step_nzv(all_predictors()) %>% # remove variables with non-zero variance
            step_impute_knn(all_predictors()) %>% # impute using KNN
            step_smote(Status) # upsample to deal with unbalanced sample sizes
          
          # create workflow
          hbox_wf <- workflow() %>% add_recipe(hbox_rec)
          
          # logistic regression model specifications
          glm_spec <- logistic_reg() %>% set_engine("glm")
          
          # create model on resamples
          glm_res <- hbox_wf %>% 
            add_model(glm_spec) %>%
            fit_resamples(
              resamples = hbox_boot,
              metrics = metric_set(roc_auc, accuracy, sensitivity, specificity),
              control = tune::control_resamples(save_pred = TRUE, verbose = TRUE)
            )
          
          # select best model to be used for evaluation
          glm_best <- glm_res %>% select_best("roc_auc")
          
          # finalize model
          hbox_final <- hbox_wf %>% 
            add_model(glm_spec) %>% 
            finalize_workflow(glm_best) %>%
            fit(df_model) %>%
            pull_workflow_fit()
          
          # view coefficients & p-values
          model_coefs <- hbox_final %>% tidy()
          
          # return list including model and a table of best model variable statistics
          list(glm_res, model_coefs)
          
        }
        
        


        ## Perform functions on all iterations of model
        
        set.seed(123)
        
        # select starting variables
        df_hbox_model1 <- df_hbox_log_reg %>% 
          dplyr::select(c(Status, Hematocrit, Lactate, avg_Hb_o2sat, Hb_o2sat_second_alpha, avg_Hb_conc, Hb_conc_second_alpha)) %>% 
          mutate(Status = factor(Status, levels = c("CM", "UM")))
        
        
        
        # Finalize model
                
                # finalize workflow
                final_rf <- ranger_workflow %>% 
                        finalize_workflow(select_best(ranger_tune, metric = "accuracy")) 
                
                hbox_fit <- last_fit(final_rf, hbox_split) # fitting to training data, evaluating on testing
                
                collect_metrics(hbox_fit) # computed on test set
                
                # collect predictions on test set
                collect_predictions(hbox_fit) %>% 
                        
                        filter(Status != .pred_class) # only got one wrong
                
        # Determine variable importance

                # new specification for ranger model
                imp_spec <- ranger_spec %>% 
                        finalize_model(select_best(ranger_tune, metric = "accuracy")) %>% 
                        set_engine("ranger", importance = "permutation")
                
                # plot importance 
                workflow() %>% 
                        add_recipe(ranger_recipe) %>% 
                        add_model(imp_spec) %>% 
                        fit(hbox_train) %>% 
                        pull_workflow_fit() %>% 
                        vip(aesthetics = list(alpha = 0.8, fill = "midnightblue"))
                


        
# Same procedure but only with NIRS data --------------------------------------------------------------
                
                
        # Data
        hbox_df2 <- df_hbox_impute %>% select(c(Status, DeoxyHb, DeoxyHb_alpha, OxyHb, OxyHb_alpha, O2sat, O2sat_alpha, THC, THC_alpha, THC_sd))
                
        # Build a model
                
                # split data set
                set.seed(123)
                hbox_split2 <- initial_split(hbox_df2, strata = Status)
                hbox_train2 <- training(hbox_split2)
                hbox_test2 <- testing(hbox_split2)  
                
                # resamples (not a lot of data, need to use cross validation)
                set.seed(234)
                vfold_cv(hbox_train, strata = Status) # test too small
                hbox_folds2 <- bootstraps(hbox_train2, strata = Status)
                
                # Scaffolding for setting up common types of models
                use_ranger(Status ~ ., data = hbox_df2)
                
                # Copy usemodels code...
                
                # Create recipe
                ranger_recipe2 <- 
                        recipe(formula = Status ~ ., data = hbox_df2) %>% 
                        
                        # BoxCox transformation to normality
                        step_BoxCox(all_predictors(), -all_outcomes()) %>% 
                        
                        # remove variables with non-zero variance
                        step_nzv(all_predictors(), -all_outcomes()) %>% 
                        
                        # impute missing data using K-nearest neighbors
                        step_knnimpute(all_predictors(), -all_outcomes()) 
                
                # Model specifications (set for tuning)
                ranger_spec2 <- 
                        rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>% 
                        set_mode("classification") %>% 
                        set_engine("ranger") 
                
                # Create workflow
                ranger_workflow2 <- 
                        workflow() %>% 
                        add_recipe(ranger_recipe2) %>% 
                        add_model(ranger_spec2) 
                
        # Tune model
                
                # Initial tune on bootstraps
                set.seed(84374)
                ranger_tune2 <-
                        tune_grid(ranger_workflow2, 
                                  resamples = hbox_folds2, 
                                  grid = 11)
                
                # Explore tuning results
                show_best(ranger_tune2, metric = "roc_auc")
                show_best(ranger_tune2, metric = "accuracy")
                
                # Visualize results
                autoplot(ranger_tune2)
                
        # Finalize model
                
                # finalize workflow
                final_rf2 <- ranger_workflow2 %>% 
                        finalize_workflow(select_best(ranger_tune2, metric = "accuracy")) 
                
                hbox_fit2 <- last_fit(final_rf2, hbox_split2) # fitting to training data, evaluating on testing
                
                collect_metrics(hbox_fit2) # computed on test set
                
                # collect predictions on test set
                collect_predictions(hbox_fit2) %>% 
                        
                        filter(Status != .pred_class) # 6 wrong
                
        # Determine variable importance
                
                # new specification for ranger model
                imp_spec2 <- ranger_spec2 %>% 
                        finalize_model(select_best(ranger_tune2, metric = "accuracy")) %>% 
                        set_engine("ranger", importance = "permutation")
                
                # plot importance 
                workflow() %>% 
                        add_recipe(ranger_recipe2) %>% 
                        add_model(imp_spec2) %>% 
                        fit(hbox_train2) %>% 
                        pull_workflow_fit() %>% 
                        vip(aesthetics = list(alpha = 0.8, fill = "midnightblue"))
                


                
# Same procedure but only with other data --------------------------------------------------------------
                
                
        # Data
        hbox_df3 <- df_hbox_impute %>% select(-c(subject_id, DeoxyHb, DeoxyHb_alpha, OxyHb, OxyHb_alpha, O2sat, O2sat_alpha, THC, THC_alpha, THC_sd))
                
        # Build a model
                
                # split data set
                set.seed(1234)
                hbox_split3 <- initial_split(hbox_df3, strata = Status)
                hbox_train3 <- training(hbox_split3)
                hbox_test3 <- testing(hbox_split3)  
                
                # resamples (not a lot of data, need to use cross validation)
                set.seed(2345)
                vfold_cv(hbox_train3, strata = Status) # test too small
                hbox_folds3 <- bootstraps(hbox_train3, strata = Status)
                
                # Scaffolding for setting up common types of models
                use_ranger(Status ~ ., data = hbox_df3)
                
                # Copy usemodels code...
                
                # Create recipe
                ranger_recipe3 <- 
                        recipe(formula = Status ~ ., data = hbox_df3) %>% 
                        
                        # BoxCox transformation to normality
                        step_BoxCox(all_predictors(), -all_outcomes()) %>% 
                        
                        # remove variables with non-zero variance
                        step_nzv(all_predictors(), -all_outcomes()) %>% 
                        
                        # impute missing data using K-nearest neighbors
                        step_knnimpute(all_predictors(), -all_outcomes()) 
                
                # Model specifications (set for tuning)
                ranger_spec3 <- 
                        rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>% 
                        set_mode("classification") %>% 
                        set_engine("ranger") 
                
                # Create workflow
                ranger_workflow3 <- 
                        workflow() %>% 
                        add_recipe(ranger_recipe3) %>% 
                        add_model(ranger_spec3) 
                
        # Tune model
                
                # Initial tune on bootstraps
                set.seed(8434)
                ranger_tune3 <-
                        tune_grid(ranger_workflow3, 
                                  resamples = hbox_folds3, 
                                  grid = 11)
                
                # Explore tuning results
                show_best(ranger_tune3, metric = "roc_auc")
                show_best(ranger_tune3, metric = "accuracy")
                
                # Visualize results
                autoplot(ranger_tune3)
                
        # Finalize model
                
                # finalize workflow
                final_rf3 <- ranger_workflow3 %>% 
                        finalize_workflow(select_best(ranger_tune3, metric = "accuracy")) 
                
                hbox_fit3 <- last_fit(final_rf3, hbox_split3) # fitting to training data, evaluating on testing
                
                collect_metrics(hbox_fit3) # computed on test set
                
                # collect predictions on test set
                collect_predictions(hbox_fit3) %>% 
                        
                        filter(Status != .pred_class) # 6 wrong
                
        # Determine variable importance
                
                # new specification for ranger model
                imp_spec3 <- ranger_spec3 %>% 
                        finalize_model(select_best(ranger_tune3, metric = "accuracy")) %>% 
                        set_engine("ranger", importance = "permutation")
                
                # plot importance 
                workflow() %>% 
                        add_recipe(ranger_recipe3) %>% 
                        add_model(imp_spec3) %>% 
                        fit(hbox_train3) %>% 
                        pull_workflow_fit() %>% 
                        vip(aesthetics = list(alpha = 0.8, fill = "midnightblue"))
                
                
                
                
                
                
                