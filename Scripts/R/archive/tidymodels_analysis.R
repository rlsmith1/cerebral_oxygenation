

### Run the LR model on variables of interest to identify variables that drive discrimination between UM and CM


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
        library(janitor)



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

        
      
        df_model <- read_csv("Data/Hans_datatable_exports/malawi key data v04Jan2022.csv") %>% 
          as_tibble() %>% 
          janitor::clean_names() %>% 
          dplyr::rename("subject_id" = "subject_id_session") %>% 
          
          dplyr::filter(status %in% c("CM", "UM") & session_no == 1 | subject_id == "TM0003CM01") %>% 
          dplyr::filter(!grepl("blood", subject_id)) %>% 
          mutate(status = factor(status, levels = c("CM", "UM"))) %>% 
          select(subject_id, status, hct, lactate, cerebral_hb_tot, cerebral_hb_tot_alpha2, cerebral_hb_oxy, cerebral_hb_oxy_alpha2) %>% 
          mutate_at(3:8, as.numeric)
        
        # look at missing data
        df_model %>% is.na %>% colSums
        df_model %>% filter(is.na(cerebral_hb_oxy_alpha2))
        
        require(naniar)
        df_model %>% select(-c(subject_id, status)) %>% gg_miss_upset() # 12 pts missing all 4 NIRS voi
        df_model %>% count(status) 
        
        # remove patients with all 4 NIRS variables missing
        df_model <- df_model %>% drop_na(contains("hb_tot"))
        df_model %>% select(-c(subject_id, status)) %>% gg_miss_upset() # 3 patients missing hb_oxy, 3 missing hb_oxy_alpha, 3 missing both, 2 missing lactate, 1 missing hct & lactate

        # add 1/10 alpha as variable
        df_model <- df_model %>% 
          mutate(one_tenth_cerebral_hb_tot_alpha = 1/10 * cerebral_hb_tot_alpha2,
                 one_tenth_cerebral_hb_oxy_alpha = 1/10 * cerebral_hb_oxy_alpha2) %>% 
          select(-c(cerebral_hb_tot_alpha2, cerebral_hb_oxy_alpha2))
        
        # resample
        set.seed(234)
        hbox_boot <- bootstraps(df_model, strata = status, times = 5000)
        
        # recipe
        hbox_rec <- recipe(status ~ ., data = df_model) %>% 
          step_impute_knn(all_predictors()) %>% 
          step_normalize(all_predictors()) %>%
          step_zv(all_predictors()) %>%
          step_smote(status) %>% 
          update_role(subject_id, new_role = "id")
        
        # create workflow
        hbox_wf <- workflow() %>% add_recipe(hbox_rec)
        
        # logistic regression model specifications
        glm_spec <- logistic_reg() %>% set_engine("glm")
        
        # create model on 5000 different resamples
        require(doParallel)
        doParallel::registerDoParallel()
        
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
          extract_fit_parsnip()
        
        # metrics
        collect_metrics(glm_res) %>% 
          dplyr::select(.metric, mean, std_err) %>% 
          mutate(mean = mean*100, std_err = std_err*100)
        
        # coefficients
        df_CI <- confint(hbox_final$fit) %>% exp() %>% as_tibble(rownames = "term")
        hbox_final %>%
          tidy(exponentiate = TRUE) %>% 
          left_join(df_CI, by = "term") %>% 
          select(term, estimate, `2.5 %`, `97.5 %`, statistic, p.value) %>% 
          filter(term != "(Intercept)") %>% 
          arrange(p.value)
        
        
        save(glm_res, hbox_final, file = "objects/glm_model.Rdata")     
                
                
                
                