

### Run the full analysis pipeline on UM patient signals that had been missing from the original datatable


# libraries ---------------------------------------------------------------


    library(tidyverse)
    library(purrr)
    library(data.table)
    library(zoo)



# read in data -------------------------------------------------------------------------


    # get names
    my_files <- list.files(path = "Data/Raw data/missing_UM", pattern = "*.txt", full.names = TRUE, recursive = FALSE)

    # read in
    l_raw_data <- 1:length(my_files) %>% 
      purrr::map(~read.delim2(file = my_files[.x], header = FALSE, sep = "\t") %>% as_tibble()) 
    
    # name with subject_id
    l_names <- 1:length(my_files) %>% 
      purrr::map(~substr(my_files[.x], start = 26, stop = 35)) 
    
    names(l_raw_data) <- l_names



# format data -------------------------------------------------------------

  
    # format to just THC and o2 sat in workable form
    l_raw_data <- 1:length(l_raw_data) %>% 
      
      purrr::map(~mutate(l_raw_data[[.x]], subject_id = l_names[[.x]])) %>% 
      purrr::map(~dplyr::filter(.x, !grepl("Patient|Data|[R|r]aw", V1))) %>% 
      purrr::map(~dplyr::select(.x, c(ncol(.x), 2:4))) %>% 
      purrr::map(~dplyr::rename(.x, "Time" = "V2", "hb_oxy" = "V3", "hb_tot" = "V4")) %>% 
      purrr::map(~mutate(.x, Time = as.numeric(Time))) %>% 
      purrr::map(~mutate(.x, hb_tot = as.numeric(hb_tot))) %>% 
      purrr::map(~mutate(.x, hb_oxy = as.numeric(hb_oxy)))
    
    # renumber patients
    l_raw_data <- 1:length(l_raw_data) %>% 
      purrr::map(~mutate(l_raw_data[[.x]], number = c(1:length(l_raw_data))[.x])) 
    
    # condense list to df
    df_raw_data <- rbindlist(l_raw_data) %>% as_tibble() 
    
    # add status
    df_raw_data <- df_raw_data %>% 
      mutate(Status = substr(subject_id, start = 7, stop = 8)) %>% 
      dplyr::select(c(number, subject_id, Status, everything()))
    
    # make sure status is UM
    df_raw_data %>% count(Status)




# remove "blips" ----------------------------------------------------------



    # remove signals that aren't at least 200s  
    sampling_rate <- 1/.02
    df_raw_data %>% group_by(number) %>% count() %>% dplyr::filter(n < 200*sampling_rate) # 8 not long enough
    df_raw_data <- df_raw_data %>% filter(number != 8)
    df_raw_data %>% count(number)
    
    
  # plotting functions
    
    # hb_tot
    f_plot_hbtot <- function(df, num) {
      
      df %>% 
        dplyr::filter(number == num) %>% 
        ggplot(aes(x = Time, y = hb_tot)) +
        geom_line() 
      
    }
    
    # hb_oxy
    f_plot_hboxy <- function(df, num) {
      
      df %>% 
        dplyr::filter(number == num) %>% 
        ggplot(aes(x = Time, y = hb_oxy)) +
        geom_line()
      
    }
    
    # segment signal function
    f_segment_signal <- function(df, num, start, end) {
      
      df %>% 
        dplyr::filter(number == num) %>% 
        dplyr::filter(Time > start & Time < end) # time is in SECONDS
    }
    
    

  # TM1015UM01
    # hb_tot
    df_raw_data %>% f_plot_hbtot(1) + theme_bw()
    df_hbtot1 <- df_raw_data %>% dplyr::select(-hb_oxy) %>% filter(hb_tot > 25 & hb_tot < 150)
    df_hbtot1 %>% f_plot_hbtot(1)
    df_hbtot1$hb_tot %>% var(na.rm = TRUE) # check var less than 1000
    
    # hb_oxy
    df_raw_data %>% f_plot_hboxy(1)
    df_hboxy1 <- df_raw_data %>% dplyr::select(-hb_tot) %>% filter(hb_oxy > -100 & hb_oxy < 100)
    df_hboxy1 %>% f_plot_hboxy(1)
    df_hboxy1$hb_oxy %>% var(na.rm = TRUE)
    
  # TM1023UM01
    # signal too noisy :(
    
  # TM1024UM01
    # hb_tot
    df_raw_data %>% f_plot_hbtot(3) + theme_bw()
    df_hbtot3 <- df_raw_data %>% dplyr::select(-hb_oxy) %>% f_segment_signal(3, start = 0, end = 750)
    df_hbtot3 %>% f_plot_hbtot(3)
    df_hbtot3$hb_tot %>% mean(na.rm = TRUE) # check var less than 3000
    
    # hb_oxy
    df_raw_data %>% f_plot_hboxy(3)
    df_hboxy3 <- df_raw_data %>% dplyr::select(-hb_tot) %>% f_segment_signal(3, start = 0, end = 750) %>% filter(hb_oxy < 400 & hb_oxy > -190)
    df_hboxy3 %>% f_plot_hboxy(3)
    df_hboxy3$hb_oxy %>% var(na.rm = TRUE)
    
  # TM1025UM01
    # hb_tot
    df_raw_data %>% f_plot_hbtot(4) + theme_bw()
    df_hbtot4 <- df_raw_data %>% dplyr::select(-hb_oxy) %>% f_segment_signal(4, start = 500, end = 1000)
    df_hbtot4 %>% f_plot_hbtot(4)
    df_hbtot4$hb_tot %>% var(na.rm = TRUE) # check var less than 4000
    
    # hb_oxy
    df_raw_data %>% f_plot_hboxy(4)
    df_hboxy4 <- df_raw_data %>% dplyr::select(-hb_tot) %>% f_segment_signal(4, start = 700, end = 950) %>% filter(hb_oxy < 200 & hb_oxy > -100)
    df_hboxy4 %>% f_plot_hboxy(4)
    df_hboxy4$hb_oxy %>% var(na.rm = TRUE)
    
  # TM1026UM01
    # signal too noisy :(
    
  # TM1027UM01
    # signal too noisy :(

  # TM1028UM01
    # hb_tot
    df_raw_data %>% f_plot_hbtot(7) + theme_bw()
    df_hbtot7 <- df_raw_data %>% dplyr::select(-hb_oxy) %>% f_segment_signal(7, start = 210, end = 800) %>% filter(hb_tot > -100)
    df_hbtot7 %>% f_plot_hbtot(7)
    df_hbtot7$hb_tot %>% mean(na.rm = TRUE)
    
    # hb_oxy
    df_raw_data %>% f_plot_hboxy(7)
    df_hboxy7 <- df_raw_data %>% dplyr::select(-hb_tot) %>% f_segment_signal(7, start = 200, end = 490) %>% filter(hb_oxy < 250 & hb_oxy > -250)
    df_hboxy7 %>% f_plot_hboxy(7)
    df_hboxy7$hb_oxy %>% mean(na.rm = TRUE)
    
  # TM1030UM01
    # hb_tot
    df_raw_data %>% f_plot_hbtot(9) + theme_bw()
    df_hbtot9 <- df_raw_data %>% dplyr::select(-hb_oxy) %>% f_segment_signal(9, start = 0, end = 1250) %>% filter(hb_tot > 50)
    df_hbtot9 %>% f_plot_hbtot(9)
    df_hbtot9$hb_tot %>% mean(na.rm = TRUE)
    
    # hb_oxy
    df_raw_data %>% f_plot_hboxy(9)
    df_hboxy9 <- df_raw_data %>% dplyr::select(-hb_tot) %>% f_segment_signal(9, start = 0, end = 1250) %>% filter(hb_oxy < 200 & hb_oxy > -190)
    df_hboxy9 %>% f_plot_hboxy(9)
    df_hboxy9$hb_oxy %>% mean(na.rm = TRUE)
    
  # TM1031UM01
    # hb_tot - signal too noisy
    
    # hb_oxy
    df_raw_data %>% f_plot_hboxy(10)
    df_hboxy10 <- df_raw_data %>% dplyr::select(-hb_tot) %>% filter(hb_oxy < 500 & hb_oxy > -100)
    df_hboxy10 %>% f_plot_hboxy(10)
    df_hboxy10$hb_oxy %>% var(na.rm = TRUE)
    
  # TM1032UM01
    # hb_tot
    df_raw_data %>% f_plot_hbtot(11) + theme_bw()
    df_hbtot11 <- df_raw_data %>% dplyr::select(-hb_oxy) %>% f_segment_signal(11, start = 800, end = 1150) %>% filter(hb_tot < 1000)
    df_hbtot11 %>% f_plot_hbtot(11)
    df_hbtot11$hb_tot %>% var(na.rm = TRUE)
    
    # hb_oxy
    df_raw_data %>% f_plot_hboxy(11)
    df_hboxy11 <- df_raw_data %>% dplyr::select(-hb_tot) %>% f_segment_signal(11, start = 800, end = 1150) %>% filter(hb_oxy < 200 & hb_oxy > -1110)
    df_hboxy11 %>% f_plot_hboxy(11)
    df_hboxy11$hb_oxy %>% var(na.rm = TRUE)
    
  # TM1033UM01
    # hb_tot
    df_raw_data %>% f_plot_hbtot(12) + theme_bw()
    df_hbtot12 <- df_raw_data %>% dplyr::select(-hb_oxy) %>% f_segment_signal(12, start = 250, end = 2000) %>% filter(hb_tot > -100)
    df_hbtot12 %>% f_plot_hbtot(12)
    df_hbtot12$hb_tot %>% mean(na.rm = TRUE)
    
    # hb_oxy
    df_raw_data %>% f_plot_hboxy(12)
    df_hboxy12 <- df_raw_data %>% dplyr::select(-hb_tot) %>% f_segment_signal(12, start = 250, end = 2000) %>% filter(hb_oxy < 200 & hb_oxy > -200)
    df_hboxy12 %>% f_plot_hboxy(12)
    df_hboxy12$hb_oxy %>% var(na.rm = TRUE)
    

    
  # combine
    df_hbtot_filt <- bind_rows(df_hbtot1, df_hbtot3, df_hbtot4, df_hbtot7, df_hbtot9, df_hbtot11, df_hbtot12)
    df_hboxy_filt <- bind_rows(df_hboxy1, df_hboxy3, df_hboxy4, df_hboxy7, df_hboxy9, df_hboxy10, df_hboxy11, df_hboxy12)
    
    

# take medians --------------------------------------------------------------


    df_medians <- df_hbtot_filt %>% 
      group_by(subject_id) %>% 
      summarise(cerebral_hb_tot = median(hb_tot)) %>% 
      # filter(number %in% c(1, 3, 4, 7, 9, 11, 12))
      left_join(df_hboxy_filt %>% 
                  group_by(subject_id) %>% 
                  summarise(cerebral_hb_oxy = median(hb_oxy)))
    
    
# format and export -------------------------------------------------------


  # add milliseconds in to Time
    
    # hb_tot
    df_hbtot_filt <- df_hbtot_filt %>% 
      group_by(number, Time) %>% 
      left_join(count(.)) %>% 
      mutate(Time = ifelse(n == 49, Time + 0.02*row_number(), Time + 0.02*(row_number() - 1))) %>% 
      dplyr::select(-n)
    
    # hb_oxy
    df_hboxy_filt <- df_hboxy_filt %>% 
      group_by(number, Time) %>% 
      left_join(count(.)) %>% 
      mutate(Time = ifelse(n == 49, Time + 0.02*row_number(), Time + 0.02*(row_number() - 1))) %>% 
      dplyr::select(-n)
    
  # export to .txt to filter in Matlab
    
    # Hb_tot
    l_hbtot_filt <- df_hbtot_filt %>% split(df_hbtot_filt$number)
    hbtot_names <- df_hbtot_filt$subject_id %>% unique()
    names(l_hbtot_filt) <- hbtot_names
    path <- "Data/signal_segments/cerebral_missing_UM/Hb_tot/"
    1:length(l_hbtot_filt) %>% map(~write.table(l_hbtot_filt[[.x]], file = paste0(path, names(l_hbtot_filt[.x]), ".txt")))
    
    # Hb_oxy
    l_hboxy_filt <- df_hboxy_filt %>% split(df_hboxy_filt$number)
    hboxy_names <- df_hboxy_filt$subject_id %>% unique()
    names(l_hboxy_filt) <- hboxy_names
    path <- "Data/signal_segments/cerebral_missing_UM/Hb_oxy/"
    1:length(l_hboxy_filt) %>% map(~write.table(l_hboxy_filt[[.x]], file = paste0(path, names(l_hboxy_filt[.x]), ".txt")))
    
    
    
  ### moving average & signal filtering are done in MatLab, read back into R for DFA   
    
 


# read back in data -------------------------------------------------------


    # Hb_tot
    hbtot_path <- "Data/filtered_signals/0.2s_mean_0.001_1_filt/cerebral_missing_UM/Hb_tot"
    hbtot_files <- list.files(hbtot_path)
    l_hbtot_filt <- 1:length(hbtot_files) %>% purrr::map(~read.table(paste0(hbtot_path, sep = "/", hbtot_files[.x])))
    names(l_hbtot_filt) <- sub("*.txt", "", hbtot_files)
    
    # Hb_oxy
    hboxy_path <- "Data/filtered_signals/0.2s_mean_0.001_1_filt/cerebral_missing_UM/Hb_oxy"
    hboxy_files <- list.files(hboxy_path)
    l_hboxy_filt <- 1:length(hboxy_files) %>% purrr::map(~read.table(paste0(hboxy_path, sep = "/", hboxy_files[.x])))
    names(l_hboxy_filt) <- sub("*.txt", "", hboxy_files)
    
    
    
    
# format data & run DFA -------------------------------------------------------------
    
    
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
    
    source("Scripts/dfa_functions.R")
    
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
    
    # ADD THESE FOR ANALYSIS/FIGURES!!!
    save(l_dfa_hboxy, l_dfa_hbtot, file = "cerebral_missing_UM_dfa_results.Rdata")
    


# second order detrending -------------------------------------------------

    
    library(segmented)
    

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
    
    

# plots and statistical analysis ------------------------------------------


    
    require(DescTools)
    require(gghalves)
    
# plot alphas for Hb_tot and Hb_oxy
    
    
  # overall
    
    # Hb_tot
    df_hbtot_all_alphas <- df_hbtot_all_alphas %>% 
      mutate(subject_id = sub("*.txt", "", names(l_hbtot_filt)),
             Status = substr(subject_id, start = 7, stop = 8),
             Status = ifelse(Status == "HV", "HC", Status),
             Status = factor(Status, levels = c("HC", "UM", "CM"))) %>% 
      filter(!is.na(Status))
    
    df_hbtot_all_alphas %>% 
      ggplot(aes(x = Status, y = overall_a, color = Status)) +  
      geom_half_point(shape = 1) +
      geom_half_boxplot() +
      ylim(0, 1.5) +
      theme_bw()
    
    kruskal.test(overall_a ~ Status, df_hbtot_all_alphas)
    # DunnTest(overall_a ~ Status, df_hbtot_all_alphas)
    
    # Hb_oxy
    df_hboxy_all_alphas <- df_hboxy_all_alphas %>% 
      mutate(subject_id = sub("*.txt", "", names(l_hboxy_filt)),
             Status = substr(subject_id, start = 7, stop = 8),
             Status = ifelse(Status == "HV", "HC", Status),
             Status = factor(Status, levels = c("HC", "UM", "CM"))) %>% 
      filter(!is.na(Status))
    
    df_hboxy_all_alphas %>% 
      ggplot(aes(x = Status, y = overall_a, color = Status)) +  
      geom_half_point(shape = 1) +
      geom_half_boxplot() +
      ylim(0, 1.5) +
      theme_bw()
    
    kruskal.test(overall_a ~ Status, df_hboxy_all_alphas)
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
    # DunnTest(short_a ~ Status, df_hbtot_all_alphas)
    
    # Hb_oxy
    df_hboxy_all_alphas %>% 
      ggplot(aes(x = Status, y = short_a, color = Status)) +  
      geom_half_point(shape = 1) +
      geom_half_boxplot() +
      ylim(0, 2) +
      theme_bw()
    
    kruskal.test(short_a ~ Status, df_hboxy_all_alphas)
    DunnTest(short_a ~ Status, df_hboxy_all_alphas)
    
  # long (use this: Chen et al 2002)
    
    # Hb_tot
    p_hbtot_long_a <- df_hbtot_all_alphas %>% 
      ggplot(aes(x = Status, y = long_a, color = Status)) +  
      geom_half_point(shape = 1) +
      geom_half_boxplot() +
      labs(y = "alpha", x = "") +
      ggtitle("Hb_tot") +
      ylim(0, 1.5) +
      theme_bw() +
      theme(legend.position = "none")
    
    kruskal.test(long_a ~ Status, df_hbtot_all_alphas)
    # DunnTest(long_a ~ Status, df_hbtot_all_alphas)
    
    # Hb_oxy
    p_hboxy_long_a <- df_hboxy_all_alphas %>% 
      ggplot(aes(x = Status, y = long_a, color = Status)) +  
      geom_half_point(shape = 1) +
      geom_half_boxplot() +
      labs(y = "alpha", x = "") +
      ggtitle("Hb_oxy") +
      ylim(0, 1.5) +
      theme_bw() +
      theme(legend.position = "none")
    
    kruskal.test(long_a ~ Status, df_hboxy_all_alphas)
    DunnTest(long_a ~ Status, df_hboxy_all_alphas)
    
    
  # medians
    df_hbtot_all_alphas %>% 
      summarise(median(overall_a, na.rm = TRUE),
                median(short_a, na.rm = TRUE),
                median(long_a, na.rm = TRUE))
    
    df_hboxy_all_alphas %>% 
      summarise(median(overall_a, na.rm = TRUE),
                median(short_a, na.rm = TRUE),
                median(long_a, na.rm = TRUE))

    
    
# export ----------------------------------------------
    
    
    df_bp_alphas <- df_hbtot_all_alphas %>% 
      rename_with(~paste0("cerebral_hbtot_", .), all_of(1:4)) %>% 
      left_join(df_hboxy_all_alphas %>% 
                  rename_with(~paste0("cerebral_hboxy_", .), all_of(1:4))) %>% 
      left_join(df_means, by = "subject_id") %>% 
      select(c(subject_id, Status, everything()))
    
    write.csv(df_bp_alphas, file = "Data/cerebral_missing_UM_alphas.csv")
    
    
    

    