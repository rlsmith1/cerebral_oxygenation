


# Libraries ---------------------------------------------------------------

        library(tidyverse)
        library(nortest)
        library(FSA)
        library(tidymodels)
        library(factoextra)
        library(plotly)
        library(usemodels)
        library(ranger)
        library(vip)



# Data --------------------------------------------------------------------

        # read in data
        df_hbox <- read.csv("Data/98 first visit records with cerebral oxygenation stats_trimmed.csv") %>% as_tibble()
        

# Format data -------------------------------------------------------------

        
        # get rid of "Other"
        df_hbox <- df_hbox %>% 
                filter(Status != "Other") %>% 
                mutate(Status = factor(Status, levels = c("HV", "UM", "CM")))
        
        # identify duplicated rows (TM0003, TM2001)
        df_hbox <- df_hbox %>% filter(duplicated(Subject.ID..NIAID.) == FALSE)
        
        # Convert age to numeric
        df_hbox <- df_hbox %>% 
                mutate(Age = gsub(" Years, ", "_", Age)) %>% 
                mutate(Age = gsub(" Months ", "", Age)) %>% 
                separate(Age, into = c("Years", "Months"), sep = "_") %>% 
                mutate(Years = as.numeric(Years)) %>% 
                mutate(Months = as.numeric(Months)/12) %>% 
                mutate(Age = Years + Months) %>% 
                select(-c(Years, Months))
        
        # Convert sex to factor
        df_hbox <- df_hbox %>% mutate(Sex = as.factor(Sex))
        
        # Trim for variables of interest (VOI)
        df_hbox_trim <- df_hbox %>% 
                dplyr::select(c(Subject.ID..NIAID., Status, Glucose, Hematocrit, Lactate, Temperature,
                                Arginine.umol.L, Haptoglobin..mg.dl., Hemoglobin..uM., Whole.Blood.Nitrite, 
                                cerebral.hbdeox, cerebral.hbdeox.alpha, cerebral.hbox, cerebral.hbox.alpha,
                                cerebral.sat, cerebral.sat.alpha, cerebral.thc, cerebral.thc.alpha, cerebral.thc.norm.std)) %>% 
                
                dplyr::rename("subject_id" = "Subject.ID..NIAID.", "Arginine" = "Arginine.umol.L", "Haptoglobin" = "Haptoglobin..mg.dl.", 
                              "Hemoglobin" = "Hemoglobin..uM.", "WB_Nitrite" = "Whole.Blood.Nitrite", "DeoxyHb" = "cerebral.hbdeox", 
                              "DeoxyHb_alpha" = "cerebral.hbdeox.alpha", "OxyHb" = "cerebral.hbox", "OxyHb_alpha" = "cerebral.hbox.alpha", 
                              "O2sat" = "cerebral.sat", "O2sat_alpha" = "cerebral.sat.alpha", "THC" = "cerebral.thc", 
                              "THC_alpha" = "cerebral.thc.alpha", "THC_sd" = "cerebral.thc.norm.std") 
        
        # impute missing data with median by group
        f_impute <- function(x, na.rm = TRUE) (replace(x, is.na(x), median(x, na.rm = na.rm)))
        
        df_hbox_impute <- df_hbox_trim %>% 
                group_by(Status) %>% 
                mutate_at(3:ncol(.), f_impute) %>% 
                ungroup()

        # pivot long
        df_hbox_long <- df_hbox_trim %>% pivot_longer(cols = 3:ncol(.), names_to = "measurement", values_to = "value")
        

        
        
        
# Explore data with plots -------------------------------------------------------

        
        # histogram to explore normality of data
        df_hbox_long %>% 
                
                ggplot(aes(x = value)) +
                geom_histogram(bins = 20) +
                facet_grid(Status ~ measurement, scales = "free")
        
        # plot all
        df_hbox_long %>% 
                
                ggplot(aes(x = Status, y = value)) +
                geom_violin(aes(color = Status)) +
                geom_jitter(position = position_jitter(0.2), shape = 1) +
                stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
                facet_wrap(~measurement, scales = "free_y") +
                
                theme_bw() +
                #ggtitle("cerebral.hbdeox") +
                theme(legend.position = "none")
   
        
        
        

# Use PCA to explore if these variables cluster Status --------------------

     

        # convert tibble to df
        df_hbox_impute_df <- df_hbox_impute %>% select(-c(subject_id, Status)) %>% as.data.frame(df_hbox_impute) 
        rownames(df_hbox_impute_df) <- df_hbox_impute$subject_id
        
        # compute PCA
        my_pca <- prcomp(df_hbox_impute_df, scale = TRUE)
        
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
        groups <- as.factor(df_hbox_impute$Status)
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
        




# Tidymodels: Random forest predictive modeling --------------------------------------------------------------


        # Data
        hbox_df <- df_hbox_trim %>% select(-subject_id)
        
        # Build a model
        
                # split data set
                set.seed(123)
                hbox_split <- initial_split(hbox_df, strata = Status)
                hbox_train <- training(hbox_split)
                hbox_test <- testing(hbox_split)  
                
                # resamples (not a lot of data, need to use cross validation)
                set.seed(234)
                vfold_cv(hbox_train, strata = Status) # test too small
                hbox_folds <- bootstraps(hbox_train, strata = Status)
                
                # Scaffolding for setting up common types of models
                use_ranger(Status ~ ., data = hbox_df)
                
                # Copy usemodels code...
                
                # Create recipe
                ranger_recipe <- 
                        recipe(formula = Status ~ ., data = hbox_df) %>% 

                        # BoxCox transformation to normality
                        step_BoxCox(all_predictors(), -all_outcomes()) %>% 
                        
                        # remove variables with non-zero variance
                        step_nzv(all_predictors(), -all_outcomes()) %>% 
                        
                        # impute missing data using K-nearest neighbors
                        step_knnimpute(all_predictors(), -all_outcomes()) 
                
                # Model specifications (set for tuning)
                ranger_spec <- 
                        rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>% 
                        set_mode("classification") %>% 
                        set_engine("ranger") 
                
                # Create workflow
                ranger_workflow <- 
                        workflow() %>% 
                        add_recipe(ranger_recipe) %>% 
                        add_model(ranger_spec) 
                
        # Tune model
                
                # Initial tune on bootstraps
                set.seed(84374)
                ranger_tune <-
                        tune_grid(ranger_workflow, 
                                  resamples = hbox_folds, 
                                  grid = 11)
                
                # Explore tuning results
                show_best(ranger_tune, metric = "roc_auc")
                show_best(ranger_tune, metric = "accuracy")
                
                # Visualize results
                autoplot(ranger_tune)
                
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
        hbox_df2 <- df_hbox_impute %>% select(-c(subject_id, Glucose, Hematocrit, Lactate, Temperature, Arginine, Haptoglobin, Hemoglobin, WB_Nitrite))
                
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
                
                
                
                
                
                
                
                

                
                
                
                
                
                
                
                
                