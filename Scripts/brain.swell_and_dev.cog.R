


# libraries ---------------------------------------------------------------


  library(tidyverse)



# data --------------------------------------------------------------------


  df_all_data <- read.csv("Data/malawi key data v11Nov2021.csv") %>% as_tibble() 

  
  

# explore data ------------------------------------------------------------

  
   sum(!is.na(df_all_data$brain.swell)) # 40 swelling scores
  
  
  df_data_filt <- df_all_data %>% 
    dplyr::select(c(Subject.ID..NIAID., Status, brain.swell, devcog01, devcog06, devcog12,
                    Hematocrit, cerebral_hb_tot, cerebral_hb_tot_alpha, cerebral_hb_tot_alpha2)) %>% 
    filter(!is.na(brain.swell) & cerebral_hb_tot != 0)
  
  df_data_filt %>% pivot_longer(3:5, names_to = "variable") %>% 
    
    # filter(Status == "CM") %>% 
    
    ggplot(aes(x = Status, y = value)) +
    geom_violin() +
    geom_point(aes(color = Subject.ID..NIAID.), shape = 1) +
    facet_wrap(~ variable, scales = "free_y") +
    theme(legend.position = "none")
 
  

# correlations ------------------------------------------------------------


  require(corrr)
  require(ggpmisc)
  
  
  # high cerebral hb_tot outliers??
  df_data_filt_outliers <- df_data_filt %>% 
    mutate(Q1_hbtot = quantile(cerebral_hb_tot, 0.25),
           Q3_hbtot = quantile(cerebral_hb_tot, 0.75),
           IQR_hbtot = Q3_hbtot - Q1_hbtot,
           outlier_hbtot = case_when(
             
             cerebral_hb_tot < Q1_hbtot - IQR_hbtot ~ 1,
             cerebral_hb_tot > Q3_hbtot + IQR_hbtot ~ 1,
             TRUE ~ 0
             
           )
           
    ) %>% 
    filter(outlier_hbtot != 1)
  
  
  # brain swelling and hb_tot (significant)
  df_data_filt_outliers %>% 
    ggplot(aes(x = cerebral_hb_tot, y = brain.swell)) +
    geom_point(shape = 1, size = 2) +
    geom_smooth(method = "lm") +
    stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., ..p.value.label.., sep = "~~~")),
                 parse = TRUE) +
    theme_bw()
  
  # hb_tot and alpha
  df_data_filt_outliers %>% 
    ggplot(aes(x = cerebral_hb_tot, y = cerebral_hb_tot_alpha2)) +
    geom_point(shape = 1, size = 2) +
    geom_smooth(method = "lm") +
    stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., ..p.value.label.., sep = "~~~")),
                 parse = TRUE) +
    theme_bw()

  # hb_tot alpha and brain.swell
  df_data_filt_outliers %>% 
    ggplot(aes(x = cerebral_hb_tot_alpha2, y = brain.swell)) +
    geom_point(shape = 1, size = 2) +
    geom_smooth(method = "lm") +
    stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., ..p.value.label.., sep = "~~~")),
                 parse = TRUE) +
    theme_bw()
  
  
  # devcog01 and brain.swell
  df_data_filt_outliers %>% 
    ggplot(aes(x = brain.swell, y = devcog01)) +
    geom_point(shape = 1, size = 2) +
    geom_smooth(method = "lm") +
    stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., ..p.value.label.., sep = "~~~")),
                 parse = TRUE) +
    theme_bw()
  
  # devcog01 and hb_tot
  df_data_filt_outliers %>% 
    ggplot(aes(x = cerebral_hb_tot, y = devcog01)) +
    geom_point(shape = 1, size = 2) +
    geom_smooth(method = "lm") +
    stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., ..p.value.label.., sep = "~~~")),
                 parse = TRUE) +
    theme_bw()
  
  # devcog01 and hb_tot alpha
  df_data_filt_outliers %>% 
    ggplot(aes(x = cerebral_hb_tot_alpha2, y = devcog01)) +
    geom_point(shape = 1, size = 2) +
    geom_smooth(method = "lm") +
    stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., ..p.value.label.., sep = "~~~")),
                 parse = TRUE) +
    geom_vline(xintercept = 0.5, lty = 2, alpha = 0.5, color = "black") +
    theme_bw()
  
  # combined models
  lm(brain.swell ~ cerebral_hb_tot + Hematocrit, data = df_data_filt_outliers) %>% summary()
  lm(devcog01 ~ cerebral_hb_tot + cerebral_hb_tot_alpha2 + Hematocrit, data = df_data_filt_outliers) %>% summary()

  # brain.swell, hb_tot, & hematocrit
  df_data_filt_outliers %>% 
    ggplot(aes(x = cerebral_hb_tot, y = brain.swell)) +
    geom_point(aes(color = Hematocrit, size = Hematocrit), shape = 1) +
    geom_smooth(method = "lm") +
    stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., ..p.value.label.., sep = "~~~")),
                 parse = TRUE) +
    theme_bw()
  
  
  # categorize alpha?
  df_data_filt_outliers <- df_data_filt_outliers %>% 
    
    mutate(hb_tot_alpha_cat = case_when(
      
      cerebral_hb_tot_alpha2 < 0.4 ~ "anti-correlated",
      cerebral_hb_tot_alpha2 > 0.4 & cerebral_hb_tot_alpha2 < 0.6 ~ "white noise",
      cerebral_hb_tot_alpha2 > 0.6 & cerebral_hb_tot_alpha2 < 0.9 ~ "correlated",
      cerebral_hb_tot_alpha2 > 0.9 & cerebral_hb_tot_alpha2 < 1.1 ~ "pink noise",
      cerebral_hb_tot_alpha2 > 1.1 ~ "non-stationary",
      
    ))
    
  df_data_filt_outliers %>% 
    ggplot(aes(x = cerebral_hb_tot, y = brain.swell, color = hb_tot_alpha_cat)) +
    geom_point(shape = 1, size = 2) +
    # geom_smooth(method = "lm") +
    # stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., ..p.value.label.., sep = "~~~")),
    #              parse = TRUE) +
    theme_bw()
  
  
  
  
  
  
  
  
  

