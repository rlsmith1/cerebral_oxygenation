

#### original code for this data processing can be found in raw_data_processing.R ####


# read in data -------------------------------------------------------------------------

    library(tidyverse)
    library(purrr)
    
    
    # get names (use data from muscle folder, looks more complete)
    my_files <- list.files(path = "Data/NIRS files of brain and muscle/files_for_muscle_analysis/follow up/", pattern = "*.txt", full.names = TRUE, recursive = FALSE)
    
    # read in
    l_follow_up <- 1:length(my_files) %>% 
      purrr::map(~read.delim2(file = my_files[.x], header = FALSE, sep = "\t") %>% as_tibble()) 
    
    # name with Status
    l_names <- 1:length(my_files) %>% 
      purrr::map(~strsplit(my_files[.x], "//")[[1]][2] %>% 
                   strsplit("[.]") %>% .[[1]] %>% .[1])
    
    names(l_follow_up) <- l_names


    

# format raw data ---------------------------------------------------------

    
    
    library(data.table)
    library(zoo)
    
    # format to just THC and o2 sat in workable form
    l_follow_up <- 1:length(l_follow_up) %>% 
      
      purrr::map(~mutate(l_follow_up[[.x]], subject_id = l_names[[.x]])) %>% 
      purrr::map(~dplyr::filter(.x, !grepl("Patient|Data|[R|r]aw", V1))) %>% 
      purrr::map(~dplyr::select(.x, c(ncol(.x), 2:4))) %>% 
      purrr::map(~dplyr::rename(.x, "Time" = "V2", "O2_sat" = "V3", "THC" = "V4")) %>% 
      purrr::map(~mutate(.x, Time = as.numeric(as.character(Time)))) %>% 
      purrr::map(~mutate(.x, THC = as.numeric(as.character(THC)))) %>% 
      purrr::map(~mutate(.x, O2_sat = as.numeric(as.character(O2_sat))))
    
    # number patients
    l_follow_up <- 1:length(l_follow_up) %>% 
      purrr::map(~mutate(l_follow_up[[.x]], number = c(1:length(l_follow_up))[.x])) 
    
    # condense list to df
    df_follow_up <- rbindlist(l_follow_up) %>% as_tibble() 
    
    # add status
    df_follow_up <- df_follow_up %>% 
      mutate(Status = as.factor("CM")) %>% 
      dplyr::select(c(number, subject_id, Status, everything()))
    

    
# remove blips ------------------------------------------------------------

    
    # remove signals that aren't at least 200s  
    sampling_rate <- 1/.02
    df_follow_up %>% group_by(number) %>% count() %>% dplyr::filter(n < 200*sampling_rate) # all signals are at least 200s

    # plotting functions
    
      # THC
      f_plot_thc <- function(df, num) {
        
        df %>% 
          dplyr::filter(number == num) %>% 
          ggplot(aes(x = Time, y = THC)) +
          geom_line()
        
      }
      
      # O2 sat
      f_plot_o2sat <- function(df, num) {
        
        df %>% 
          dplyr::filter(number == num) %>% 
          ggplot(aes(x = Time, y = O2_sat)) +
          geom_line()
        
      }
      
    # segment signal function
      f_segment_signal <- function(df, num, start, end) {
        
        df %>% 
          dplyr::filter(number == num) %>% 
          dplyr::filter(Time > start & Time < end) # time is in SECONDS
      }
    
      
   ## variance needs to be < 10000 ##

  # 1: TM0005CMF
      
      # O2_sat
      df_follow_up %>% f_plot_o2sat(1)
      df_o2sat_fu1 <- df_follow_up %>% dplyr::select(-THC) %>% f_segment_signal(1, start = 0, end = 600)
      df_o2sat_fu1 %>% f_plot_o2sat(1) 
      df_o2sat_fu1$O2_sat %>% var(na.rm = TRUE)
      
      # THC
      df_follow_up %>% f_plot_thc(1) 
      df_thc_fu1 <- df_follow_up %>% dplyr::select(-O2_sat) %>% f_segment_signal(1, start = 1000, end = 2000)
      df_thc_fu1 %>% f_plot_thc(1) 
      df_thc_fu1$THC %>% var(na.rm = TRUE)
      
      
  # 2: TM0007CMF
      
      # O2_sat
      df_follow_up %>% f_plot_o2sat(2)
      df_o2sat_fu2 <- df_follow_up %>% dplyr::select(-THC) %>% f_segment_signal(2, start = 0, end = 2000)
      df_o2sat_fu2 %>% f_plot_o2sat(2) 
      df_o2sat_fu2$O2_sat %>% var(na.rm = TRUE)
      
      # THC
      df_follow_up %>% f_plot_thc(2) 
      df_thc_fu2 <- df_follow_up %>% dplyr::select(-O2_sat) %>% f_segment_signal(2, start = 0, end = 2000)
      df_thc_fu2 %>% f_plot_thc(2) 
      df_thc_fu2$THC %>% var(na.rm = TRUE)
      
      
  # 3: TM0009CMF
      
      # O2_sat
      df_follow_up %>% f_plot_o2sat(3)
      df_o2sat_fu3 <- df_follow_up %>% dplyr::select(-THC) %>% f_segment_signal(3, start = 1000, end = 2000)
      df_o2sat_fu3 %>% f_plot_o2sat(3) 
      df_o2sat_fu3$O2_sat %>% var(na.rm = TRUE)
      
      # THC
      df_follow_up %>% f_plot_thc(3) 
      df_thc_fu3 <- df_follow_up %>% dplyr::select(-O2_sat) %>% f_segment_signal(3, start = 1000, end = 2000) %>% dplyr::filter(THC < 200)
      df_thc_fu3 %>% f_plot_thc(3) 
      df_thc_fu3$THC %>% var(na.rm = TRUE)
      
      
  # 4: TM0018CMF
      
      # O2_sat
      df_follow_up %>% f_plot_o2sat(4)
      df_o2sat_fu4 <- df_follow_up %>% dplyr::select(-THC) %>% f_segment_signal(4, start = 0, end = 1150)
      df_o2sat_fu4 %>% f_plot_o2sat(4) 
      df_o2sat_fu4$O2_sat %>% var(na.rm = TRUE)
      
      # THC
      df_follow_up %>% f_plot_thc(4) 
      df_thc_fu4 <- df_follow_up %>% dplyr::select(-O2_sat) %>% f_segment_signal(4, start = 0, end = 1150)
      df_thc_fu4 %>% f_plot_thc(4) 
      df_thc_fu4$THC %>% var(na.rm = TRUE)
      
      
  # 5: TM0019CMF
      
      # O2_sat
      df_follow_up %>% f_plot_o2sat(5)
      df_o2sat_fu5 <- df_follow_up %>% dplyr::select(-THC) %>% f_segment_signal(5, start = 500, end = 1000)
      df_o2sat_fu5 %>% f_plot_o2sat(5) 
      df_o2sat_fu5$O2_sat %>% var(na.rm = TRUE)
      
      # THC
      df_follow_up %>% f_plot_thc(5) 
      df_thc_fu5 <- df_follow_up %>% dplyr::select(-O2_sat) %>% f_segment_signal(5, start = 0, end = 2000)
      df_thc_fu5 %>% f_plot_thc(5) 
      df_thc_fu5$THC %>% var(na.rm = TRUE)
      
      
  # 6: TM0020CMF
      
      # O2_sat
      df_follow_up %>% f_plot_o2sat(6)
      df_o2sat_fu6 <- df_follow_up %>% dplyr::select(-THC) %>% f_segment_signal(6, start = 0, end = 500) %>% dplyr::filter(O2_sat < 200)
      df_o2sat_fu6 %>% f_plot_o2sat(6) 
      df_o2sat_fu6$O2_sat %>% var(na.rm = TRUE)
      
      # THC
      df_follow_up %>% f_plot_thc(6) 
      df_thc_fu6 <- df_follow_up %>% dplyr::select(-O2_sat) %>% f_segment_signal(6, start = 0, end = 250)
      df_thc_fu6 %>% f_plot_thc(6) 
      df_thc_fu6$THC %>% var(na.rm = TRUE)
      
      
  # 7: TM0021CMF
      
      # O2_sat
      df_follow_up %>% f_plot_o2sat(7)
      df_o2sat_fu7 <- df_follow_up %>% dplyr::select(-THC) %>% dplyr::filter(O2_sat < 200 & O2_sat > -100)
      df_o2sat_fu7 %>% f_plot_o2sat(7) 
      df_o2sat_fu7$O2_sat %>% var(na.rm = TRUE)
      
      # THC
      df_follow_up %>% f_plot_thc(7) 
      df_thc_fu7 <- df_follow_up %>% dplyr::select(-O2_sat) %>% f_segment_signal(6, start = 0, end = 800)
      df_thc_fu7 %>% f_plot_thc(7) 
      df_thc_fu7$THC %>% var(na.rm = TRUE)
      
      
  # 8: TM0022CMF
      
      # O2_sat
      df_follow_up %>% f_plot_o2sat(8)
      df_o2sat_fu8 <- df_follow_up %>% dplyr::select(-THC) %>% dplyr::filter(O2_sat < 200 & O2_sat > -100)
      df_o2sat_fu8 %>% f_plot_o2sat(8) 
      df_o2sat_fu8$O2_sat %>% var(na.rm = TRUE)
      
      # THC - signal too noisy, variance too high

      
  # 9: TM0023CMF
      
      # O2_sat - signal too noisy

      # THC - signal too noisy

      
  # 10: TM0024CMF
      
      # O2_sat - signal too noisy

      # THC - signal too noisy, variance too high

      
  # 11: TM0026CMF
      
      # O2_sat - signal too noisy

      # THC - signal too noisy

      
  # 12: TM0027CMF
      
      # O2_sat
      df_follow_up %>% f_plot_o2sat(12)
      df_o2sat_fu12 <- df_follow_up %>% dplyr::select(-THC) %>% f_segment_signal(12, start = 250, end = 1300)
      df_o2sat_fu12 %>% f_plot_o2sat(12) 
      df_o2sat_fu12$O2_sat %>% var(na.rm = TRUE)
      
      # THC
      df_follow_up %>% f_plot_thc(12) 
      df_thc_fu12 <- df_follow_up %>% dplyr::select(-O2_sat) %>% f_segment_signal(12, start = 250, end = 1300)
      df_thc_fu12 %>% f_plot_thc(12) 
      df_thc_fu12$THC %>% var(na.rm = TRUE)
      
      
  # 13: TM0028CMF
      
      # O2_sat - signal too noisy

      # THC - signal too noisy

      
  # 14: TM0030CMF
      
      # O2_sat - signal too noisy

      # THC - signal too noisy

      
  # 15: TM0032CMF
      
      # O2_sat
      df_follow_up %>% f_plot_o2sat(15)
      df_o2sat_fu15 <- df_follow_up %>% dplyr::select(-THC) %>% f_segment_signal(15, start = 1250, end = 1500) %>% dplyr::filter(O2_sat < 3000)
      df_o2sat_fu15 %>% f_plot_o2sat(15) 
      df_o2sat_fu15$O2_sat %>% var(na.rm = TRUE)
      
      # THC
      df_follow_up %>% f_plot_thc(15) 
      df_thc_fu15 <- df_follow_up %>% dplyr::select(-O2_sat) %>% f_segment_signal(15, start = 1250, end = 1500)
      df_thc_fu15 %>% f_plot_thc(15) 
      df_thc_fu15$THC %>% var(na.rm = TRUE)
      
      
  # 16: TM0060CMF
      
      # O2_sat - signal too noisy

      # THC - signal too noisy

      
      # Once you have repeated this for all the patients, combine the filtered Hb_tot signals into one tibble and the filtered Hb_oxy signals into another
      
      df_o2sat_fu_filt <- bind_rows(df_o2sat_fu1, df_o2sat_fu2, df_o2sat_fu3, df_o2sat_fu4, df_o2sat_fu5,
                                    df_o2sat_fu6, df_o2sat_fu7, df_o2sat_fu8,
                                    df_o2sat_fu12, df_o2sat_fu15,
                                    )
      
      df_thc_fu_filt <- bind_rows(df_thc_fu1, df_thc_fu2, df_thc_fu3, df_thc_fu4, df_thc_fu5,
                                  df_thc_fu6, df_thc_fu7, 
                                  df_thc_fu12, df_thc_fu15,
                                  )
      

      

      

# determine avg THC & O2 sat ----------------------------------------------


            
      # total THC
      df_THC_fu_mean <- df_thc_fu_filt %>% 
        dplyr::filter(!is.na(THC)) %>% 
        group_by(number, Status) %>% 
        summarise(THC_mean = mean(THC))
      
      # plot
      df_THC_fu_mean  %>%  
        ggplot(aes(x = Status, y = THC_mean)) +
        geom_violin(aes(color = Status)) +
        geom_jitter(position = position_jitter(0.2), shape = 1) +
        stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
        
        theme_bw() +
        labs(y = "uM") +
        ggtitle("Cerebral tissue hemoglobin concentration in CM follow-up patients") +
        theme(legend.position = "none") 
      
      
      # O2 sat
      df_o2sat_fu_mean <- df_o2sat_fu_filt %>% 
        dplyr::filter(!is.na(O2_sat)) %>% 
        group_by(number, Status) %>% 
        summarise(o2sat_mean = mean(O2_sat))
      
      # plot
      df_o2sat_fu_mean %>%  
        ggplot(aes(x = Status, y = o2sat_mean)) +
        geom_violin(aes(color = Status)) +
        geom_jitter(position = position_jitter(0.2), shape = 1) +
        stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
        
        theme_bw() +
        labs(y = "percent") +
        ggtitle("Cerebral tissue hemoglobin oxygen saturation in CM follow-up patients") +
        theme(legend.position = "none")
      


# savitzky-golay signal filtering -----------------------------------------

      
      # add milliseconds in to Time
      
      # THC
      df_thc_fu_filt <- df_thc_fu_filt %>%
        group_by(number, Time) %>% 
        left_join(count(.)) %>% 
        mutate(Time = ifelse(n == 49, Time + 0.02*row_number(), Time + 0.02*(row_number() - 1))) %>% 
        dplyr::select(-n)
      
      # O2 sat
      df_o2sat_fu_filt <- df_o2sat_fu_filt %>% 
        group_by(number, Time) %>% 
        left_join(count(.)) %>% 
        mutate(Time = ifelse(n == 49, Time + 0.02*row_number(), Time + 0.02*(row_number() - 1))) %>% 
        dplyr::select(-n)
      
      
      # 4th order savitzy-golay filter
      library(signal)

  # THC
      df_thc_fu_filt <- df_thc_fu_filt %>% 
        group_by(number) %>% 
        mutate(sg_filt = sgolayfilt(THC, p = 3)) %>% 
        ungroup()
      
      # plot
      df_thc_fu_filt %>% dplyr::filter(number == 5) %>% 
        ggplot(aes(x = Time)) +
        geom_line(aes(y = THC)) +
        geom_line(aes(y = sg_filt), color = "red")
      
  # O2 sat
      df_o2sat_fu_filt <- df_o2sat_fu_filt %>% 
        group_by(number) %>% 
        mutate(sg_filt = sgolayfilt(O2_sat, p = 3)) %>% 
        ungroup()
      
      # plot
      df_o2sat_fu_filt %>% dplyr::filter(number == 5) %>% 
        ggplot(aes(x = Time)) +
        geom_line(aes(y = O2_sat)) +
        geom_line(aes(y = sg_filt), color = "red")
      

      
      

# high-pass signal filtering ----------------------------------------------


      # use a high-pass filter to remove anything below 0.01 Hz (associated with head displacements and motion noise) 
      
      library(seewave)
      library(dplR)
      library(data.table)
      
      # filter out anything lower than 0.01 Hz (noise)
      
      df_thc_fu_filt_pass %>% group_by(number) %>% mutate(mean = mean(THC))
      
      # THC
      l_thc_fu_filt <- df_thc_fu_filt %>% split(df_thc_fu_filt$number)
      sampling_rate <- 50
      butter_filter <- butter(4, W = 0.01/(sampling_rate/2), type = "high")
      l_thc_fu_filt_pass <- 1:length(l_thc_fu_filt) %>% purrr::map(~mutate(l_thc_fu_filt[[.x]], high_pass = filtfilt(butter_filter, sg_filt)))
      
      df_thc_fu_filt_pass <- rbindlist(l_thc_fu_filt_pass) %>% as_tibble()
      
      # O2 sat
      l_o2sat_fu_filt <- df_o2sat_fu_filt %>% split(df_o2sat_fu_filt$number)
      l_o2sat_fu_filt_pass <- 1:length(l_o2sat_fu_filt) %>% purrr::map(~mutate(l_o2sat_fu_filt[[.x]], high_pass = filtfilt(butter_filter, sg_filt)))
      
      df_o2sat_fu_filt_pass <- rbindlist(l_o2sat_fu_filt_pass) %>% as_tibble()
      
  # plot levels of signal processing
      df_thc_fu_filt_pass %>% dplyr::filter(number == 1) %>% 
        ggplot(aes(x = Time)) +
        geom_line(aes(y = THC)) +
        geom_line(aes(y = sg_filt), color = "red") +
        geom_line(aes(y = high_pass), color = "blue") +
        theme_bw()
      
      df_o2sat_fu_filt_pass %>% dplyr::filter(number == 1) %>% 
        ggplot(aes(x = Time)) +
        geom_line(aes(y = O2_sat)) +
        geom_line(aes(y = sg_filt), color = "red") +
        geom_line(aes(y = high_pass), color = "blue") +
        theme_bw()
      
      
      

# save filtered signal dfs -----------------------------------------------------------------------
      
      
      # save filtered signal
      save(df_thc_fu_filt, df_thc_fu_filt_pass, 
           df_o2sat_fu_filt, df_o2sat_fu_filt_pass,
           file = "follow_up_filtered_signals.Rdata")
      
      

# ### DFA -----------------------------------------------------------------

#### original code can be found in DFA.R

      
    load("follow_up_filtered_signals.Rdata")
      


# cumulative sum ----------------------------------------------------------

      
  ## THC
      # create list, where each object is a patient
      l_thc_fu_filt_pass <- df_thc_fu_filt_pass %>% split(f = df_thc_fu_filt_pass$number)
      
      # take the mean of each filtered signal
      l_fu_means <- 1:length(l_thc_fu_filt_pass) %>% purrr::map(~mean(l_thc_fu_filt_pass[[.x]]$sg_filt, na.rm = TRUE))
      
      # subtract the mean and take the cumulative sum
      l_thc_fu_cumsum <- 1:length(l_fu_means) %>% purrr::map(~mutate(l_thc_fu_filt_pass[[.x]], mean = l_fu_means[[.x]]))
      df_thc_fu_cumsum <- rbindlist(l_thc_fu_cumsum) %>% as_tibble()
      
      df_thc_fu_cumsum <- df_thc_fu_cumsum %>% 
        mutate(subtract_mean = sg_filt - mean) %>% 
        group_by(number) %>% 
        mutate(cumsum = cumsum(subtract_mean)) %>% 
        ungroup()
      
      # plot
      df_thc_fu_cumsum %>% dplyr::filter(number == 1) %>% 
        ggplot(aes(x = Time)) +
        geom_line(aes(y = cumsum)) +
        geom_line(aes(y = sg_filt), color = "red") +
        
        theme_bw() +
        ggtitle("Cumulative sum of Savitzy-Golar smoothed signal")
      
      
  ## O2 sat
      # create list, where each object is a patient
      l_o2sat_fu_filt_pass <- df_o2sat_fu_filt_pass %>% split(f = df_o2sat_fu_filt_pass$number)
      
      # take the mean of each filtered signal
      l_fu_means <- 1:length(l_o2sat_fu_filt_pass) %>% purrr::map(~mean(l_o2sat_fu_filt_pass[[.x]]$sg_filt, na.rm = TRUE))
      
      # subtract the mean and take the cumulative sum
      l_o2sat_fu_cumsum <- 1:length(l_fu_means) %>% purrr::map(~mutate(l_o2sat_fu_filt_pass[[.x]], mean = l_fu_means[[.x]]))
      df_o2sat_fu_cumsum <- rbindlist(l_o2sat_fu_cumsum) %>% as_tibble()
      
      df_o2sat_fu_cumsum <- df_o2sat_fu_cumsum %>% 
        mutate(subtract_mean = sg_filt - mean) %>% 
        group_by(number) %>% 
        mutate(cumsum = cumsum(subtract_mean)) %>% 
        ungroup()
      
      # plot
      df_o2sat_fu_cumsum %>% dplyr::filter(number == 1) %>% 
        ggplot(aes(x = Time)) +
        geom_line(aes(y = cumsum)) +
        geom_line(aes(y = sg_filt), color = "red") +
        
        theme_bw() +
        ggtitle("Cumulative sum of Savitzy-Golar smoothed signal")
      
      
      

# create DFA functions ----------------------------------------------------

     
      # determine window lengths
      f_window_sizes <- function(df, num_of_windows) {
        
        low_end <- 0.1
        high_end <- floor(nrow(df)/50*0.10) # take 10% of signal length & round down to nearest integer
        log_size_range <- log10(high_end) - log10(low_end)
        log_step_size <- log_size_range/num_of_windows # divide by number of windows desired
        
        10^(seq(log10(low_end), log10(high_end), by = log_step_size))
        
        
      }
      
      # split by window sizes (window length in s)
      f_split_windows <- function(df, window_length) {
        
        window_size <- window_length*50
        signal_length <- nrow(df)
        reps <- rep(1:ceiling(signal_length/window_size), each = window_size)[1:signal_length]
        list <- df %>% split(reps) 
        1:length(list) %>% purrr::map(~mutate(list[[.x]], window_size = window_length))
        
      }
      
      # compute linear least squares regression line for all windows
      f_lm_least_squares <- function(l_splits) {
        
        1:length(l_splits) %>% purrr::map(~lm(cumsum ~ Time, data = l_splits[[.x]]))
        
      }
      
      # detrend by calculating residuals
      f_calc_residuals <- function(l_splits, l_lms) {
        
        1:length(l_lms) %>% purrr::map(~tibble(window_size = l_splits[[.x]]$window_size, 
                                               Time = l_splits[[.x]]$Time, 
                                               residuals = l_lms[[.x]]$residuals))  
        
      }

      # plot window function
      f_plot_window <- function(l_splits, window_num) {
        
        l_splits[[window_num]] %>% 
          ggplot(aes(x = Time)) +
          geom_line(aes(y = sg_filt), color = "red", alpha = 0.5) +
          geom_line(aes(y = cumsum), color = "black", alpha = 0.65) +
          geom_smooth(aes(y = cumsum), method = "lm", color = "blue", se = FALSE) +
          theme_bw()
        
      }
      
      # plot detrended window
      f_plot_resids <- function(l_resids, window_num) {
        
        l_resids[[window_num]] %>% 
          ggplot(aes(x = Time, y = residuals)) +
          geom_line(color = "black") +
          geom_hline(yintercept = 0, color = "blue") +
          theme_bw()
        
      }
      
      # return window size, number, and fluctuation
      f_calc_fluct <- function(l_resids) {
        
        1:length(l_resids) %>% 
          purrr::map(~tibble(window_size = l_resids[[.x]]$window_size, 
                             window_num = .x, 
                             fluctuation = sd(l_resids[[.x]]$residuals))) %>% 
          rbindlist() %>% 
          as_tibble()  
        
      }
      
      # calculate average fluctuation for each window size  
      f_calc_avg_fluct <- function(l_fluct) {
        
        1:length(l_fluct) %>% purrr::map(~tibble(window_size = l_fluct[[.x]]$window_size, 
                                                 avg_fluctuation = mean(l_fluct[[.x]]$fluctuation, na.rm = TRUE)) %>% 
                                           unique()) %>% 
          rbindlist() %>% as_tibble()
        
      }
      

      ## combine all into one DFA function
      f_dfa <- function(df, num_of_windows) {
        
        
        # determine the window lengths needed for n number of windows equally spaced on a log scale across the signal
        window_lengths <- f_window_sizes(df, num_of_windows)
        
        # split df into a list of lists, where length of the overall list is the number of window lengths, and each list item contains windows of that length
        l_splits <- window_lengths %>% purrr::map(~f_split_windows(df, .x))
        names(l_splits) <- window_lengths
        
        # compute linear least squares regression line for all windows
        l_lms <- 1:length(l_splits) %>% purrr::map(~f_lm_least_squares(l_splits[[.x]]))
        names(l_lms) <- window_lengths
        
        # detrend by calculating and subtractin out residuals
        l_resids <- 1:length(l_splits) %>% purrr::map(~f_calc_residuals(l_splits[[.x]], l_lms[[.x]]))
        names(l_resids) <- window_lengths
        
        # calculate fluctuations (standard deviation) for each window per window size 
        l_fluct <- 1:length(l_resids) %>% purrr::map(~f_calc_fluct(l_resids[[.x]]) %>% unique())
        
        # calculate average fluctuation for each window size  
        df_avg_fluct <- f_calc_avg_fluct(l_fluct)
        
        # find line of best fit and extract alpha
        alpha <- lm(log10(avg_fluctuation) ~ log10(window_size), data = df_avg_fluct) %>% .$coefficients %>% .[2]
        
        # return variables of interest
        return(list(df_avg_fluct, alpha))
        
        
        
      }
      
      


# perform DFA on follow-up data -------------------------------------------

      
      # create list of df cumsum, where each object is a patient
      l_thc_fu_cumsum <- df_thc_fu_cumsum %>% split(f = df_thc_fu_cumsum$number)
      l_o2sat_fu_cumsum <- df_o2sat_fu_cumsum %>% split(f = df_o2sat_fu_cumsum$number)
      
      
      # repeat for all patients!!
      
      library(doParallel)
      registerDoParallel()
      
      
      l_dfa_thc_fu <- 1:length(l_thc_fu_cumsum) %>% purrr::map(~f_dfa(l_thc_fu_cumsum[[.x]], num_of_windows = 10))
      l_dfa_o2sat_fu <- 1:length(l_o2sat_fu_cumsum) %>% purrr::map(~f_dfa(l_o2sat_fu_cumsum[[.x]], num_of_windows = 10))
      
      
      # save DFA analysis results
      save(l_dfa_thc_fu, l_dfa_o2sat_fu, file = "dfa_follow_up_results.Rdata")   
      
      


# second order detrending -------------------------------------------------

      
      library(segmented)
      
      # extract the second slope after the "bend" in the alpha slope 
      
      f_segment_slope <- function(l_dfa) {
        
        df <- l_dfa[[1]] %>% mutate(log10_avg_fluct = log10(avg_fluctuation), log10_window_size = log10(window_size))
        fit <- lm(log10_avg_fluct ~ log10_window_size, data = df) 
        segfit <- segmented(fit)
        
        if(length(segfit$coefficients) == 2){
          
          a2 <- segfit$coefficients[2]
          
        } else if (length(segfit$coefficients) == 4) {
          
          a2 <- slope(segfit)$log10_window_size[2,1]
          
        }
        
        return(a2)
        
      }
      
      
      # THC
      thc_fu_alphas <- 1:length(l_dfa_thc_fu) %>% purrr::map(~l_dfa_thc_fu[[.x]][[2]]) %>% unlist()
      
      l_thc_fu_second_alpha <- 1:length(l_dfa_thc_fu) %>% purrr::map(~f_segment_slope(l_dfa_thc_fu[[.x]]))
      thc_fu_second_alpha <- 1:length(l_dfa_thc_fu) %>% purrr::map(~f_segment_slope(l_dfa_thc_fu[[.x]])) %>% unlist()
      
      # combine dfa results into df
      df_thc_fu_dfa <- df_thc_fu_cumsum %>% 
        count(number, Status) %>% 
        dplyr::select(-n) %>% 
        mutate(alpha = thc_fu_alphas, 
               second_alpha = thc_fu_second_alpha,
               Status = as.factor(Status))
      
      
      # O2 sat
      o2sat_fu_alphas <- 1:length(l_dfa_o2sat_fu) %>% purrr::map(~l_dfa_o2sat_fu[[.x]][[2]]) %>% unlist()
      
      l_o2sat_fu_second_alpha <- 1:length(l_dfa_o2sat_fu) %>% purrr::map(~f_segment_slope(l_dfa_o2sat_fu[[.x]]))
      o2sat_fu_second_alpha <- 1:length(l_dfa_o2sat_fu) %>% purrr::map(~f_segment_slope(l_dfa_o2sat_fu[[.x]])) %>% unlist()
      
      # combine dfa results into df
      df_o2sat_fu_dfa <- df_o2sat_fu_cumsum %>% 
        count(number, Status) %>% 
        dplyr::select(-n) %>% 
        mutate(alpha = o2sat_fu_alphas, 
               second_alpha = o2sat_fu_second_alpha,
               Status = as.factor(Status))
      
      
      

# plot alphas -------------------------------------------------------------

      
      
      # THC
      df_thc_fu_dfa %>% 
        ggplot(aes(x = Status, y = second_alpha)) +
        geom_violin(aes(color = Status)) +
        geom_jitter(position = position_jitter(0.2), shape = 1) +
        stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
        
        labs(y = "alpha") +
        ggtitle("Cerebral tissue hemoglobin concentration alpha follow-up") +
        
        theme_bw() +
        theme(legend.position = "none")
      
      # O2 sat
      df_o2sat_fu_dfa %>% 
        ggplot(aes(x = Status, y = second_alpha)) +
        geom_violin(aes(color = Status)) +
        geom_jitter(position = position_jitter(0.2), shape = 1) +
        stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
        
        labs(y = "alpha") +
        ggtitle("Cerebral tissue hemoglobin oxygen saturation alpha follow-up") +
        
        theme_bw() +
        theme(legend.position = "none")
      
      
      


# combine numbers with patient ID and save results ------------------------

      
      # pair number with subject id
      df_num_id <- df_thc_fu_filt %>% distinct(number, subject_id)
      
      # combine all results
      df_fu_results <- df_num_id %>% 
        left_join(df_thc_fu_dfa, by = "number") %>% 
        dplyr::rename("Hb_conc_fu_alpha" = "alpha", "Hb_conc_fu_second_alpha" = "second_alpha") %>% 
        left_join(df_o2sat_fu_dfa, by = c("number", "Status")) %>% 
        rename("Hb_o2sat_fu_alpha" = "alpha", "Hb_o2sat_fu_second_alpha" = "second_alpha")
        
      # add means
      df_fu_results <- df_fu_results %>% 
        left_join(df_THC_fu_mean) %>% 
        left_join(df_o2sat_fu_mean) %>% 
        dplyr::rename("avg_Hb_conc_fu" = "THC_mean", "avg_Hb_o2sat_fu" = "o2sat_mean") %>% 
        mutate(visit = factor("02")) %>% 
        
        # reorder
        dplyr::select(number, subject_id, visit, Status, 
                      avg_Hb_conc_fu, Hb_conc_fu_alpha, Hb_conc_fu_second_alpha,
                      avg_Hb_o2sat_fu, Hb_o2sat_fu_alpha, Hb_o2sat_fu_second_alpha)
      
      # save for figure plotting
      save(df_fu_results, file = "follow_up_results.Rdata")
      
      # write csv
      write.csv(df_fu_results, file = "Data/follow_up_results.csv")
      
      
      
      

# plot patient trends -----------------------------------------------------
 

      
      load("results.Rdata") # from DFA.Rdata
      
      # combine visit 1 with follow-ups
      df_results_all <- df_results %>% 
        mutate(subject_id = substr(subject_id, start = 1, stop = 6)) %>% 
        left_join(df_fu_results %>% 
                    mutate(subject_id = substr(subject_id, start = 1, stop = 6)) %>% 
                    select(-c(number, Status)), by = "subject_id") 
        
      # reorganize for plotting
      df_results_all_plot <- df_results_all %>% 
        dplyr::select(-visit) %>% 
        dplyr::filter(!is.na(avg_Hb_o2sat_fu)) %>% 
        pivot_longer(cols = 4:ncol(.), names_to = "variable", values_to = "value") %>% 
        mutate(follow_up = ifelse(grepl("fu", variable), "follow-up", "initial visit")) %>% 
        mutate(follow_up = factor(follow_up, levels = c("initial visit", "follow-up")))
        
      # plot THC
      df_results_all_plot %>% dplyr::filter(variable %in% c("Hb_conc_overall_alpha", "Hb_conc_fu_overall_alpha")) %>% 
        
        ggplot(aes(x = follow_up, y = value)) +
        geom_violin(aes(color = follow_up)) +
        geom_point(shape = 1) +
        stat_summary(fun = "median", geom = "crossbar", aes(x = follow_up, color = follow_up), size = 0.2, width = 0.5) +
        geom_path(aes(group = subject_id), lty = 2, alpha = 0.7) +
        
        geom_segment(aes(x = 1, y = 1.4, xend = 2, yend = 1.4), size = 0.3) +
        geom_segment(aes(x = 1, y = 1.4, xend = 1, yend = 1.35), size = 0.3) +
        geom_segment(aes(x = 2, y = 1.4, xend = 2, yend = 1.35), size = 0.3) +
        annotate(geom = "text", x = 1.5, y = 1.43, label = paste0("p = ", round(thc_fu_pval, 3))) +

        ylim(0.5, 1.5) +
        labs(y = "alpha", x = "") +
        ggtitle("Cerebral tissue hemoglobin concentration alpha follow-up visits") +
        
        theme_bw() +
        theme(legend.position = "none")
      
      # plot O2 sat
      df_results_all_plot %>% dplyr::filter(variable %in% c("Hb_o2sat_overall_alpha", "Hb_o2sat_fu_alpha")) %>% 
        
        ggplot(aes(x = follow_up, y = value)) +
        geom_violin(aes(color = follow_up)) +
        geom_point(shape = 1) +
        stat_summary(fun = "median", geom = "crossbar", aes(x = follow_up, color = follow_up), size = 0.2, width = 0.5) +
        geom_path(aes(group = subject_id), lty = 2, alpha = 0.7) +
        
        geom_segment(aes(x = 1, y = 1.4, xend = 2, yend = 1.4), size = 0.3) +
        geom_segment(aes(x = 1, y = 1.4, xend = 1, yend = 1.35), size = 0.3) +
        geom_segment(aes(x = 2, y = 1.4, xend = 2, yend = 1.35), size = 0.3) +
        annotate(geom = "text", x = 1.5, y = 1.43, label = paste0("p = ", round(o2sat_fu_pval, 3))) +
        
        ylim(0.5, 1.5) +
        labs(y = "alpha", x = "") +
        ggtitle("Cerebral tissue hemoglobin oxygen saturation alpha follow-up visits") +
        
        theme_bw() +
        theme(legend.position = "none")
      
      
      # descriptive statistics
      
        # THC
        df_thc_compare <- df_results_all_plot %>% 
          dplyr::filter(variable %in% c("Hb_conc_overall_alpha", "Hb_conc_fu_alpha"))
          
        thc_fu_pval <- wilcox.test(value ~ follow_up, data = df_thc_compare) %>% .$p.value
      
        # O2 sat
        df_o2sat_compare <- df_results_all_plot %>% 
          dplyr::filter(variable %in% c("Hb_o2sat_overall_alpha", "Hb_o2sat_fu_alpha"))
        
        o2sat_fu_pval <- wilcox.test(value ~ follow_up, data = df_o2sat_compare) %>% .$p.value
        
      
         
      
      
      
      
      
      