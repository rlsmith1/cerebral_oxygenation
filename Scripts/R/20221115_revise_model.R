

#### Plot correlations of clinical variables across each patient group to determine what variables to include in the LR model


# libraries -----------------------------------------------------------------

  
  library(tidyverse)
  library(RColorBrewer)
  library(gridExtra)
  library(patchwork)
  library(Hmisc)
  library(corrr)


# set theme for plots -----------------------------------------------------

  theme_set(theme_bw() +
              theme(plot.title = element_text(size = 18),
                    axis.title = element_text(size = 15),
                    axis.text = element_text(size = 12),
                    strip.text = element_text(size = 15),
                    legend.title = element_text(size = 15),
                    legend.text = element_text(size = 12)))
  



# data --------------------------------------------------------------------


  df_hbox <- read_csv("Data/Hans_datatable_exports/malawi key data v04Jan2022.csv") %>% 
    as_tibble() %>% 
    janitor::clean_names() %>% 
    mutate(status = replace(status, status == "HV", "HC")) %>% 
    dplyr::rename("subject_id" = "subject_id_session") %>% 
    
    dplyr::filter(status %in% c("CM", "HC", "UM") & session_no == 1 | subject_id == "TM0003CM01") %>% 
    dplyr::filter(!grepl("blood", subject_id)) %>% 
    mutate(status = factor(status, levels = c("HC", "UM", "CM"))) %>% 
    mutate_at(c("bp_dias", "bp_sys", "glucose", "lactate"), as.numeric) %>% # converts remaining character variables to numerics
    select(subject_id, status, everything()) %>% 
    group_by(status)


# format data -------------------------------------------------------------


# define clinical variables of interest
  clinical_voi <- c("age_calc", "temperature", "hr", "bp_sys", "bp_dias", "resp_rate", "pulse_ox", "hct", "lactate", "glucose", "aa_arg")
  hbox_voi <- c("cerebral_hb_tot", "cerebral_hb_tot_alpha2")

# select clinical VOI
  df_model <- df_hbox %>% 
    dplyr::select(c(status, sex, all_of(hbox_voi), all_of(clinical_voi))) %>% 
    mutate(sex_coded = ifelse(sex == "Female", 1, 0)) %>% 
    dplyr::select(-sex)
  
# impute missing data with median by group
  f_impute <- function(x, na.rm = TRUE) (replace(x, is.na(x), median(x, na.rm = na.rm)))
  
  df_model_impute <- df_model %>% 
    group_by(status) %>% 
    mutate_at(2:ncol(.), f_impute) %>% 
    ungroup() 
  
 

# plot alpha data ---------------------------------------------------------


  
# create color vector
  my_col <- c("CM" = "#00BFC4", "HC" = "#7CAE00", "UM" = "#F8766D")
  
# extract hb_tot and hb_oxy for the brain and muscle
  df_fig2 <- df_hbox %>% 
    select(contains("alpha2")) %>%
    ungroup() %>% 
    pivot_longer(2:ncol(.), names_to = "variable", values_to = "value") %>% 
    separate(variable, into = c("tissue", "variable"), sep = "_", extra = "merge") %>% 
    mutate(variable = str_remove(variable, "2"))
  
# replot - adjust y-axis
  p1 <- df_fig2 %>% 
    
    ggplot(aes(x = status, y = value)) +
    geom_violin(aes(color = status)) +
    geom_jitter(position = position_jitter(0.15), shape = 1, size = 2) +
    stat_summary(fun = "median", geom = "crossbar", aes(color = status), size = 0.2, width = 0.5) +
    facet_grid(variable ~ tissue, scales = "free_y") +
    
    # geom_hline(yintercept = 0.4, lty = 2, alpha = 0.7) +
    # geom_hline(yintercept = 0.6, lty = 2, alpha = 0.7) +
    # geom_hline(yintercept = 1.0, lty = 2, alpha = 0.7) +
    # geom_hline(yintercept = 1.5, lty = 2, alpha = 0.7) +
    
    scale_color_manual(values = my_col) +
    ylim(c(0, 1.2)) +
    labs(x = "") +
    ggtitle("Fig 2 adjust y-axis") +
    theme_bw() +
    theme(plot.title = element_text(size = 18),
          legend.position = "none",
          strip.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.title.y = element_blank(),
          axis.text = element_text(size = 15))
   
# transparent fill?  
  p2 <- df_fig2 %>% 
    
    ggplot(aes(x = status, y = value)) +
    geom_violin(aes(color = status), fill = "transparent") +
    geom_jitter(position = position_jitter(0.15), shape = 1, size = 2) +
    stat_summary(fun = "median", geom = "crossbar", aes(color = status), size = 0.2, width = 0.5) +
    facet_grid(variable ~ tissue, scales = "free_y") +
    
    # geom_hline(yintercept = 0.4, lty = 2, alpha = 0.7) +
    # geom_hline(yintercept = 0.6, lty = 2, alpha = 0.7) +
    # geom_hline(yintercept = 1.0, lty = 2, alpha = 0.7) +
    # geom_hline(yintercept = 1.5, lty = 2, alpha = 0.7) +
    
    scale_color_manual(values = my_col) +
    ylim(c(0, 1.2)) +
    labs(x = "") +
    ggtitle("Fig 2 adjust y-axis & transparent violins") +
    theme_bw() +
    theme(plot.title = element_text(size = 18),
          legend.position = "none",
          strip.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.title.y = element_blank(),
          axis.text = element_text(size = 15))
  
  
  

# plot distributions ------------------------------------------------------


  p3 <- df_model %>% 
    pivot_longer(2:ncol(.), names_to = "variable", values_to = "value") %>% 
    filter(!(variable == "age_calc" & value > 20000)) %>% 
    
    ggplot(aes(x = value)) +
    geom_density(aes(fill = status), alpha = 0.5) +
    facet_wrap(~ variable, scales = "free") +
    ggtitle("Distributions of variables of interest")
    
  df_hbox %>% filter(age_calc > 20000)  # ??????

  
  
# correlations ------------------------------------------------------------

  
  df_cor_r <- df_model_impute %>% 
    group_by(status) %>% 
    nest() %>% 
    mutate(corr = map(.x = data, 
                      .f = ~ .x %>% 
                        correlate %>% 
                        dplyr::select(term, contains("cerebral_hb_tot")) %>% 
                        filter(term %in% clinical_voi) 
    )
    ) %>% 
    unnest(cols = c(corr)) %>% 
    pivot_longer(4:5, names_to = "voi", values_to = "r") %>% 
    mutate(voi = ifelse(voi == "cerebral_hb_tot", "Hb_tot", "Hb_tot_a"))
  
  df_cor_p <- df_model_impute %>% 
    group_by(status) %>% 
    nest() %>% 
    mutate(corr = map(.x = data, 
                      .f = ~ .x %>% 
                        as.data.frame %>% 
                        as.matrix %>% 
                        rcorr %>% 
                        .$P %>% 
                        corrr::as_cordf() %>% 
                        dplyr::select(term, contains("cerebral_hb_tot")) %>% 
                        filter(term %in% clinical_voi) 
    )
    ) %>% 
    unnest(cols = c(corr)) %>% 
    pivot_longer(4:5, names_to = "voi", values_to = "p") %>% 
    mutate(voi = ifelse(voi == "cerebral_hb_tot", "Hb_tot", "Hb_tot_a"))
  

  p4 <- ggplot() +
    geom_tile(data = df_cor_p %>% filter(p < 0.05),
              mapping = aes(x = voi, y = term), 
              color = "black") +
    geom_tile(data = df_cor_r,
              mapping = aes(x = voi, 
                            y = term,
                            fill = r),
              height = 0.95, width = 0.95) +
    geom_text(data = df_cor_r,
              mapping = aes(x = voi, 
                            y = term,
                            label = round(r, 2))) +
    facet_wrap( ~ status) +
    labs(x = "", y = "") +
    scale_fill_gradient2(high = "#B2182B", mid = "white", low = "#2166AC",
                         limits = c(-1, 1)) +
    ggtitle("Correlations of clinical variables with cerebral Hb_tot") +
    theme_classic() +
    theme(plot.title = element_text(size = 18),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12))
  
  # plot significant correlations
  df_cor_p %>% 
    filter(p < 0.05) %>% 
    left_join(df_cor_r)
  
  p_um1 <- df_hbox %>% 
    filter(status == "UM") %>% 
    dplyr::select(cerebral_hb_tot, hr) %>% 
    
    ggplot(aes(x = cerebral_hb_tot, y = hr, color = status)) +
    geom_point(size = 2, show.legend = FALSE) +
    geom_smooth(method = "lm", color = "black", lty = 2, se = FALSE) +
    annotate(geom = "text", 
             label = "r = -0.44; p < 0.01",
             x = 300, y = 160,
             size = 7) + 
    scale_color_manual(values = my_col) +
    ggtitle("UM HR & Hb_tot")
    
  p_um2 <- df_hbox %>% 
    filter(status == "UM") %>% 
    dplyr::select(cerebral_hb_tot, hct) %>% 
    
    ggplot(aes(x = cerebral_hb_tot, y = hct, color = status)) +
    geom_point(size = 2, show.legend = FALSE) +
    geom_smooth(method = "lm", color = "black", lty = 2, se = FALSE) +
    annotate(geom = "text",
             label = "r = 0.39; p = 0.03",
             x = 300, y = 27,
             size = 7) +
    scale_color_manual(values = my_col) +
    ggtitle("UM Hematocrit & Hb_tot")
  
  p_hc1 <- df_hbox %>% 
    filter(status == "HC") %>% 
    dplyr::select(cerebral_hb_tot, aa_arg) %>% 
    
    ggplot(aes(x = cerebral_hb_tot, y = aa_arg, color = status)) +
    geom_point(size = 2, show.legend = FALSE) +
    geom_smooth(method = "lm", color = "black", lty = 2, se = FALSE) +
    annotate(geom = "text",
             label = "r = 0.51; p < 0.01",
             x = 95, y = 45,
             size = 7) +
    scale_color_manual(values = my_col) +
    ggtitle("HC arg & Hb_tot")
  
  p_hc2 <- df_hbox %>% 
    filter(status == "HC") %>% 
    dplyr::select(cerebral_hb_tot_alpha, aa_arg) %>% 
    
    ggplot(aes(x = cerebral_hb_tot_alpha, y = aa_arg, color = status)) +
    geom_point(size = 2, show.legend = FALSE) +
    geom_smooth(method = "lm", color = "black", lty = 2, se = FALSE) +
    annotate(geom = "text",
             label = "r = -0.43; p = 0.02",
             x = 1, y = 45,
             size = 7) +
    scale_color_manual(values = my_col) +
    ggtitle("HC arg & Hb_tot_a")
  
  p5 <- (p_um1 + p_um2) / (p_hc1 + p_hc2)
  

# save all 
  pdf("Outputs/20221115_revised_plots.pdf", width = 12, height = 10)
  list(p1, p2, p3, p4, p5)
  dev.off()

    
    
    
  
