


# libraries ---------------------------------------------------------------


  library(tidyverse)
  require(readxl)



# data --------------------------------------------------------------------


  df_all_data <- read_xlsx("Data/malawi key data v13Dec2021.xlsx") %>% as_tibble() 
  df_alpha_data <- read.csv("Data/master_datatable.csv") %>% as_tibble() %>% select(-X)
  
  

# explore data ------------------------------------------------------------

  
   sum(!is.na(df_all_data$brain.swell)) # 40 swelling scores
  df_all_data %>% colnames()
  
  # extract brain swell, devcog scores, and parasite data
  df_data_voi <- df_all_data %>% 
    dplyr::select(c(Subject.ID..NIAID., Status, session.no,
                    Hematocrit, parasitemia, hrp2, Lactate,
                    brain.swell, devcog01, devcog06, devcog12, outcome)) %>% 
    filter(Status == "CM" & session.no == 1) %>% 
    rename("subject_id" = "Subject.ID..NIAID.")
  
  # extract alpha data
  df_alpha_voi <- df_alpha_data %>% 
    dplyr::select(c(subject_id, Status,
                    cerebral_hb_oxy, cerebral_hboxy_long_a, cerebral_hb_tot, cerebral_hbtot_long_a)) %>% 
    filter(Status == "CM") %>% 
    mutate(subject_id = str_replace(subject_id, "CM01", ""))
  
  # combine alpha and other data
  df_voi <- df_data_voi %>% 
    left_join(df_alpha_voi, by = c("subject_id", "Status")) %>% 
    dplyr::select(-session.no) %>% 
    mutate(Lactate = as.numeric(Lactate))
  
  
  df_voi %>% pivot_longer(3:ncol(.), names_to = "variable") %>% 
    
    # filter(Status == "CM") %>% 
    
    ggplot(aes(x = Status, y = value)) +
    geom_violin() +
    geom_point(aes(color = subject_id), shape = 1) +
    facet_wrap(~ variable, scales = "free_y") +
    theme(legend.position = "none")
 
  

# correlations ------------------------------------------------------------


  require(corrr)
  require(ggpmisc)
  
  # hrp2 vs parasitemia
  df_voi %>% 
    ggplot(aes(x = hrp2, y = parasitemia)) +
    geom_point()
  

  # high cerebral hb_tot outliers??
  df_voi %>% 
    ggplot(aes(x = cerebral_hb_tot)) +
    geom_histogram()
  
  df_voi_outliers <- df_voi %>% 
    mutate(Q1_hbtot = quantile(cerebral_hb_tot, 0.25, na.rm = TRUE),
           Q3_hbtot = quantile(cerebral_hb_tot, 0.75, na.rm = TRUE),
           IQR_hbtot = Q3_hbtot - Q1_hbtot,
           outlier_hbtot = case_when(
             
             cerebral_hb_tot < Q1_hbtot - IQR_hbtot ~ 1,
             cerebral_hb_tot > Q3_hbtot + IQR_hbtot ~ 1,
             TRUE ~ 0
             
           )
           
    ) %>% 
    filter(outlier_hbtot != 1)
  
  
  # brain swelling and hb_tot (almost significant)
  df_voi_outliers %>% 
    ggplot(aes(x = cerebral_hb_tot, y = brain.swell)) +
    geom_point(shape = 1, size = 2) +
    geom_smooth(method = "lm") +
    stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., ..p.value.label.., sep = "~~~")),
                 parse = TRUE) +
    theme_bw()
  
  # hb_tot and alpha (ns)
  df_voi_outliers %>% 
    ggplot(aes(x = cerebral_hb_tot, y = cerebral_hbtot_long_a)) +
    geom_point(shape = 1, size = 2) +
    geom_smooth(method = "lm") +
    stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., ..p.value.label.., sep = "~~~")),
                 parse = TRUE) +
    theme_bw()

  # hb_tot alpha and brain.swell (ns)
  df_voi_outliers %>% 
    ggplot(aes(x = cerebral_hbtot_long_a, y = brain.swell)) +
    geom_point(shape = 1, size = 2) +
    geom_smooth(method = "lm") +
    stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., ..p.value.label.., sep = "~~~")),
                 parse = TRUE) +
    theme_bw()
  
  
  # devcog01 and brain.swell
  df_voi_outliers %>% 
    ggplot(aes(x = brain.swell, y = devcog01)) +
    geom_point(shape = 1, size = 2) +
    geom_smooth(method = "lm") +
    stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., ..p.value.label.., sep = "~~~")),
                 parse = TRUE) +
    theme_bw()
  
  # devcog01 and hb_tot
  df_voi_outliers %>% 
    ggplot(aes(x = cerebral_hb_tot, y = devcog01)) +
    geom_point(shape = 1, size = 2) +
    geom_smooth(method = "lm") +
    stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., ..p.value.label.., sep = "~~~")),
                 parse = TRUE) +
    theme_bw()
  
  # devcog01 and hb_tot alpha
  df_voi_outliers %>% 
    ggplot(aes(x = cerebral_hbtot_long_a, y = devcog01)) +
    geom_point(shape = 1, size = 2) +
    geom_smooth(method = "lm") +
    stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., ..p.value.label.., sep = "~~~")),
                 parse = TRUE) +
    geom_vline(xintercept = 0.5, lty = 2, alpha = 0.5, color = "black") +
    theme_bw()
  
  # brain.swell and outcome
  df_voi_outliers %>% 
    ggplot(aes(x = brain.swell, y = outcome)) +
    geom_point(shape = 1, size = 2) +
    geom_smooth(method = "lm") +
    stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., ..p.value.label.., sep = "~~~")),
                 parse = TRUE) +
    geom_vline(xintercept = 0.5, lty = 2, alpha = 0.5, color = "black") +
    theme_bw()

  

  # combined models
  lm(brain.swell ~ cerebral_hb_tot + Hematocrit, data = df_voi_outliers) %>% summary()
  lm(brain.swell ~ cerebral_hb_tot + Hematocrit + hrp2, data = df_voi_outliers) %>% summary()
  lm(brain.swell ~ cerebral_hb_tot + Hematocrit + parasitemia, data = df_voi_outliers) %>% summary()
  
  lm(devcog01 ~ cerebral_hb_tot + cerebral_hbtot_long_a + Hematocrit, data = df_voi_outliers) %>% summary()

  lm(outcome ~ cerebral_hbtot_long_a, data = df_voi_outliers) %>% summary()
  

  # brain.swell, hb_tot, & hematocrit
  df_voi_outliers %>% 
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
  
  
  
  
  
  
  
  
  

