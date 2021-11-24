

# libraries ---------------------------------------------------------------


  require(tidyverse)
  require(purrr)
  require(data.table)



# read in data ------------------------------------------------------------

  # Hb_tot
  hbtot_path <- "Data/filtered_signals/0.2s_mean_0.001_1_filt/Hb_tot"
  hbtot_files <- list.files(hbtot_path)
  l_hbtot_filt <- 1:length(hbtot_files) %>% purrr::map(~read.table(paste0(hbtot_path, sep = "/", hbtot_files[.x])))
  names(l_hbtot_filt) <- sub("*.txt", "", hbtot_files)
  
  # Hb_oxy
  hboxy_path <- "Data/filtered_signals/0.2s_mean_0.001_1_filt/Hb_oxy"
  hboxy_files <- list.files(hboxy_path)
  l_hboxy_filt <- 1:length(hboxy_files) %>% purrr::map(~read.table(paste0(hboxy_path, sep = "/", hboxy_files[.x])))
  names(l_hboxy_filt) <- sub("*.txt", "", hboxy_files)
  


# format data -------------------------------------------------------------

  
  # Hb_tot
    l_hbtot_filt <- 1:length(l_hbtot_filt) %>% 
      purrr::map(~as_tibble(l_hbtot_filt[[.x]]) %>%
                   mutate(Time = seq(0, nrow(.)/50, by = 0.02)[-1]) %>% 
                   dplyr::rename("hbtot_filt" = "V1") %>% 
                   mutate(subject_id = names(l_hbtot_filt[.x])) %>% 
                   mutate(mean = mean(hbtot_filt),
                          subtract_mean = hbtot_filt - mean,
                          cumsum = cumsum(subtract_mean))) %>% 
      rbindlist() %>% 
      as_tibble() %>% 
      group_by(subject_id) %>% 
      filter(grepl("1$", subject_id), n() > 200*50) %>% 
      split(.$subject_id)
  
    require(doParallel)
    registerDoParallel()
    
    sampling_freq <- 50
    l_dfa_hbtot <- 1:length(l_hbtot_filt) %>% purrr::map(~f_dfa(l_hbtot_filt[[.x]], 10))
  
  # Hb_oxy
    l_hboxy_filt <- 1:length(l_hboxy_filt) %>% 
      purrr::map(~as_tibble(l_hboxy_filt[[.x]]) %>%
                   mutate(Time = seq(0, nrow(.)/50, by = 0.02)[-1]) %>% 
                   dplyr::rename("hboxy_filt" = "V1") %>% 
                   mutate(subject_id = names(l_hboxy_filt[.x])) %>% 
                   mutate(mean = mean(hboxy_filt),
                          subtract_mean = hboxy_filt - mean,
                          cumsum = cumsum(subtract_mean))) %>% 
      rbindlist() %>% 
      as_tibble() %>% 
      group_by(subject_id) %>% 
      filter(grepl("1$", subject_id), n() > 200*50) %>% 
      split(.$subject_id)
    

  l_dfa_hboxy <- 1:length(l_hboxy_filt) %>% purrr::map(~f_dfa(l_hboxy_filt[[.x]], 10))
  
  # save DFA results
    save(l_dfa_hboxy, l_dfa_hbtot, file = "0.01_filt_signal_dfa.Rdata")
    save(l_dfa_hboxy_0.001, l_dfa_hbtot_0.001, file = "0.001_filt_signal_dfa.Rdata")
    save(l_dfa_hboxy_0.001_1, l_dfa_hbtot_0.001_1, file = "0.001_1_filt_signal_dfa.Rdata")
  
    # USE THESE FOR ANALYSIS/FIGURES!!!
    save(l_dfa_hboxy, l_dfa_hbtot, file = "bandpass_filtered_dfa_results.Rdata")
    
    
# plots n stuff -----------------------------------------------------------


  1:length(l_dfa_hboxy) %>% 
    purrr::map(~l_dfa_hboxy[[.x]][[1]] %>% 
                 mutate(subject_id = sub("*.txt", "", hboxy_files[.x]))) %>% 
    rbindlist() %>% 
    as_tibble() %>% 
    dplyr::mutate(Status = substr(subject_id, start = 7, stop = 8)) %>% 
    #filter(window_size < 1/0.001) %>% 
    
    ggplot(aes(x = log10(window_size), y = log10(avg_fluctuation))) +
    geom_point() +
    facet_wrap(~subject_id)
  
  1:length(l_dfa_hbtot) %>% 
    purrr::map(~l_dfa_hbtot[[.x]][[1]] %>% 
                 mutate(subject_id = sub("*.txt", "", hbtot_files[.x]))) %>% 
    rbindlist() %>% 
    as_tibble() %>% 
    dplyr::mutate(Status = substr(subject_id, start = 7, stop = 8)) %>% 
    # filter(window_size < 1/0.001) %>% 

    ggplot(aes(x = log10(window_size), y = log10(avg_fluctuation))) +
    geom_point() +
    facet_wrap(~subject_id)

  
  # alphas
  
      require(segmented)
      require(gghalves)
      
      f_segment_slope <- function(l_dfa, filter_thresh = 0.001) {
        
        df <- l_dfa[[1]] %>% 
          mutate(log10_avg_fluct = log10(avg_fluctuation), log10_window_size = log10(window_size)) #%>% 
          #dplyr::filter(window_size < 1/filter_thresh)
        fit <- lm(log10_avg_fluct ~ log10_window_size, data = df) 
        segfit <- segmented(fit)
        
        if(length(segfit$coefficients) == 2){
          
          short_a <- NA
          long_a <- NA
          breakpoint <- NA
          
        } else if (length(segfit$coefficients) == 4) {
          
          short_a <- slope(segfit)$log10_window_size[1,1]
          long_a <- slope(segfit)$log10_window_size[2,1]
          breakpoint <- summary(segfit)$psi[2]
          
        }
        
        return(tibble(overall_a = l_dfa[[2]], short_a = short_a, long_a = long_a, breakpoint = breakpoint))
        
      }
      
      df_hbtot_all_alphas <- 1:length(l_dfa_hbtot) %>% 
        purrr::map(~f_segment_slope(l_dfa_hbtot[[.x]])) %>% 
        rbindlist() %>% 
        as_tibble()
      
      df_hboxy_all_alphas <- 1:length(l_dfa_hboxy) %>% 
        purrr::map(~f_segment_slope(l_dfa_hboxy[[.x]])) %>% 
        rbindlist() %>% 
        as_tibble()
  
  # overall
  
      # Hb_tot
      df_hbtot_all_alphas <- df_hbtot_all_alphas %>% 
        mutate(subject_id = sub("*.txt", "", names(l_hbtot_filt)),
               Status = substr(subject_id, start = 7, stop = 8),
               Status = ifelse(Status == "HV", "HC", Status),
               Status = factor(Status, levels = c("HC", "UM", "CM"))) 
      
      df_hbtot_all_alphas %>% 
        ggplot(aes(x = Status, y = overall_a, color = Status)) +  
        geom_half_point(shape = 1) +
        geom_half_boxplot() +
        ylim(0, 1.5) +
        theme_bw()
    
      require(DescTools)
      kruskal.test(overall_a ~ Status, df_hbtot_all_alphas)
      DunnTest(overall_a ~ Status, df_hbtot_all_alphas)
      
      # Hb_oxy
      df_hboxy_all_alphas <- df_hboxy_all_alphas %>% 
        mutate(subject_id = sub("*.txt", "", names(l_hboxy_filt)),
               Status = substr(subject_id, start = 7, stop = 8),
               Status = ifelse(Status == "HV", "HC", Status),
               Status = factor(Status, levels = c("HC", "UM", "CM"))) 
      
      df_hboxy_all_alphas %>% 
        ggplot(aes(x = Status, y = overall_a, color = Status)) +  
        geom_half_point(shape = 1) +
        geom_half_boxplot() +
        ylim(0, 1.5) +
        theme_bw()
      
      require(DescTools)
      kruskal.test(overall_a ~ Status, df_hboxy_all_alphas_0.001_1)
      # DunnTest(overall_a ~ Status, df_hboxy_all_alphas_0.001_1)
      
  # short
      
      # Hb_tot
      df_hbtot_all_alphas %>% 
        ggplot(aes(x = Status, y = short_a, color = Status)) +  
        geom_half_point(shape = 1) +
        geom_half_boxplot() +
        ylim(0, 2) +
        theme_bw()
      
      kruskal.test(short_a ~ Status, df_hbtot_all_alphas)
      DunnTest(short_a ~ Status, df_hbtot_all_alphas)
      
      # Hb_oxy
      df_hboxy_all_alphas %>% 
        ggplot(aes(x = Status, y = short_a, color = Status)) +  
        geom_half_point(shape = 1) +
        geom_half_boxplot() +
        ylim(0, 2) +
        theme_bw()
      
      kruskal.test(short_a ~ Status, df_hboxy_all_alphas)
      

  # long
      
      # Hb_tot
      df_hbtot_all_alphas %>% 
        ggplot(aes(x = Status, y = long_a, color = Status)) +  
        geom_half_point(shape = 1) +
        geom_half_boxplot() +
        ylim(0, 1.5) +
        theme_bw()
      
      kruskal.test(long_a ~ Status, df_hbtot_all_alphas)
      DunnTest(long_a ~ Status, df_hbtot_all_alphas)
      
      # Hb_oxy
      df_hboxy_all_alphas %>% 
        ggplot(aes(x = Status, y = long_a, color = Status)) +  
        geom_half_point(shape = 1) +
        geom_half_boxplot() +
        ylim(0, 1.5) +
        theme_bw()
      
      kruskal.test(long_a ~ Status, df_hboxy_all_alphas)
      
  # medians
      df_hboxy_all_alphas %>% 
        group_by(Status) %>% 
        summarise(median(overall_a, na.rm = TRUE),
                  median(short_a, na.rm = TRUE),
                  median(long_a, na.rm = TRUE))
      
      df_hbtot_all_alphas %>% 
        group_by(Status) %>% 
        summarise(median(overall_a, na.rm = TRUE),
                  median(short_a, na.rm = TRUE),
                  median(long_a, na.rm = TRUE))
      
  

# test against sims -------------------------------------------------------


      load("sims.Rdata")
      l_wn_dfa %>% f_segment_slope()
  
      wn_alpha <- 1:length(l_wn_dfa) %>% 
        purrr::map(~l_wn_dfa[[.x]][[2]]) %>% 
        unlist() 
      
      cm_hbtot_long_alpha <- df_hbtot_all_alphas %>% filter(Status == "CM") %>% .$long_a
      cm_hboxy_long_alpha <- df_hboxy_all_alphas%>% filter(Status == "CM") %>% .$long_a
      
      wilcox.test(wn_alpha, cm_hboxy_long_alpha)
      wilcox.test(wn_alpha, cm_hbtot_long_alpha)
  
  
  