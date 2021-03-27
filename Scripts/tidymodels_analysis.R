


# libraries ---------------------------------------------------------------


        library(tidyverse)
        library(nortest)
        library(FSA)
        library(tidymodels)
        library(factoextra)
        library(plotly)



# data --------------------------------------------------------------------


        # set wd
        setwd("~")
        setwd("Documents/cerebral oxygenation/")
        
        # read in data
        df_hbox <- read.csv("data exports for analysis/98 first visit records with cerebral oxygenation stats_trimmed.csv") %>% as_tibble()
        

# format data -------------------------------------------------------------

        
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
        
        # Trim for interesting variables
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
        

        
        
# Exploratory plots -------------------------------------------------------

        
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
        
        # PCA (can we group by these variables?)
        
                # convert tibble to df
                df_hbox_impute_df <- df_hbox_impute %>% select(-c(subject_id, Status)) %>% as.data.frame(df_hbox_impute) 
                rownames(df_hbox_impute_df) <- df_hbox_impute$subject_id
                
                # compute PCA
                my_pca <- prcomp(df_hbox_impute_df, scale = TRUE)
                
                # visualize eigenvalues
                fviz_eig(my_pca)
                
                # Graph of individuals
                fviz_pca_ind(my_pca,
                             col.ind = "cos2", # Color by the quality of representation
                             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                             repel = TRUE     # Avoid text overlapping
                )
                
                # Graph of variables
                p_var <- fviz_pca_var(my_pca,
                             col.var = "contrib", # Color by contributions to the PC
                             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                             repel = TRUE     # Avoid text overlapping
                )
                
                # Biplot of individuals and variables
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
                
                # plot
                groups <- as.factor(df_hbox_impute$Status)
                p_pca <- fviz_pca_ind(my_pca,
                             col.ind = groups, # color by groups
                             palette = c("#00AFBB",  "#FC4E07", "#7CAE00"),
                             addEllipses = TRUE, # Concentration ellipses
                             ellipse.type = "confidence",
                             legend.title = "Groups",
                             repel = TRUE
                )
                
                p_pca$layers[[4]] <- NULL # get rid of text labels
                
                df_pca_contrib <- pca_res_var$contrib %>% 
                        as_tibble(rownames = "variable") %>% 
                        arrange(-Dim.1, -Dim.2)
                
                # 3D PCA
                df_pca <- pca_res_ind$coord %>% 
                        as_tibble(rownames = "subject") %>% 
                        mutate(Status = df_hbox_impute$Status)
                
                plot_ly(x = df_pca$Dim.1, y = df_pca$Dim.2, z = df_pca$Dim.3, 
                        type = "scatter3d", mode = "markers", color = df_pca$Status)
                
                # biggest contributors
                df_contrib <- tibble(Dim.1 = paste0(df_pca_contrib$variable, 
                                      sep = ": ", 
                                      round(df_pca_contrib$Dim.1, 2)),
                       Dim.2 = paste0(arrange(df_pca_contrib, -Dim.2)$variable, 
                                      sep = ": ", 
                                      round(arrange(df_pca_contrib, -Dim.2)$Dim.2, 2)),
                       Dim.3 = paste0(arrange(df_pca_contrib, -Dim.3)$variable, 
                                      sep = ": ", 
                                      round(arrange(df_pca_contrib, -Dim.3)$Dim.3, 2)))

        # SVM
        
                # overview
                glm(Status ~ ., data = df_hbox_trim, family = quasibinomial) %>% summary()
                


# Boxcox transformation - mass package? which metrics lead me to one status vs another, ordinal regression, 
# svm to separate different variables, k-nearest neighbors, PCA
# do parametric tests on transformed data
# random forest - how well can you predict these things
# won't find interaction if you do one at a time
df_hbox_long %>% filter(measurement == "cerebral.hbdeox") %>% 
        
        ggplot(aes(x = value^2)) +
        geom_histogram()


count(df_hbox, Status)
