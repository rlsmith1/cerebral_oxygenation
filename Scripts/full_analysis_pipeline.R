


# Below is the code for reading in, processing, and performing detrended fluctuation analysis on hemoglobin concentration and
# hemoglobin oxygen saturation signals from OxiplexTS. Code is commented to describe the purpose of each step.




# read in & format data ---------------------------------------------------


### READ IN DATA

    # libraries needed for this step (and future steps)
    library(tidyverse)
    library(purrr)
    
    # get names of patient files, mine are in the directory Data/Raw data/ within my R project
    my_files <- list.files(path = "Data/Raw data/", pattern = "*.txt", full.names = TRUE, recursive = FALSE)
    
    # remove duplicate patient files 
    # check through patient files, make sure there are no duplicates of patients & session numbers
  
    # read in .txt files - output is a list of tibbles, where each tibble is a patient reading
    l_raw_data <- 1:length(my_files) %>% 
      purrr::map(~read.delim2(file = my_files[.x], header = FALSE, sep = "\t") %>% as_tibble()) 
    
    # name each tibble in list with the patient ID 
    # (as per the file name, may need to adjust start & stop if using Malawi as opposed to NIAID ID)
    l_names <- 1:length(my_files) %>% 
      purrr::map(~substr(my_files[.x], start = 16, stop = 25)) 
    
    names(l_raw_data) <- l_names

### FORMAT DATA
  
    # load required libraries
    library(data.table)
    library(zoo)
    
    # format each tibble in list: 
    # output is a list of tibbles each containing subject ID (from .txt file name), Time, Hb_oxy (cerebral Hb O2 sat), & Hb_tot (cerebal Hb conc)
    l_raw_data <- 1:length(l_raw_data) %>% 
      
      purrr::map(~mutate(l_raw_data[[.x]], subject_id = l_names[[.x]])) %>% 
      purrr::map(~dplyr::filter(.x, !grepl("Patient|Data|[R|r]aw", V1))) %>% 
      purrr::map(~dplyr::select(.x, c(ncol(.x), 2:4))) %>% 
      purrr::map(~dplyr::rename(.x, "Time" = "V2", "Hb_oxy" = "V3", "Hb_tot" = "V4")) %>% 
      purrr::map(~mutate(.x, Time = as.numeric(Time))) %>% 
      purrr::map(~mutate(.x, Hb_tot = as.numeric(Hb_tot))) %>% 
      purrr::map(~mutate(.x, Hb_oxy = as.numeric(Hb_oxy)))
    
    # renumber patients for use in processing steps
    l_raw_data <- 1:length(l_raw_data) %>% 
      purrr::map(~mutate(l_raw_data[[.x]], number = c(1:length(l_raw_data))[.x])) 
    
    # condense list to df
    df_raw_data <- rbindlist(l_raw_data) %>% as_tibble() 
    
    # add patient status as a factor column, extracted from patient id
    df_raw_data <- df_raw_data %>% 
      mutate(Status = substr(subject_id, start = 7, stop = 8)) %>% 
      dplyr::select(c(number, subject_id, Status, everything())) %>% 
      mutate(Status = factor(Status, levels = c("HC", "UM", "CM")))
    
    # make sure patient status is one of HC, UM, or CM. adjust those that aren't
    df_raw_data %>% count(Status)

    # final output is one tibble combining all patient readings with the following columns:
        # number (int - patient number), subject_id (chr - patient ID  extracted from .txt file name),
        # Status (fct - extracted from subject_id), Time (dbl - from Oxiplex reading), Hb_oxy (dbl -
        # Hb oxygen saturation measured by Oxiplex), Hb_tot (dbl - Hb concentration measured by Oxiplex)
    
    # may be useful to save df_raw_data at this point using the save() function so you don't have to read everything in again

  

    
    


# signal preprocessing ------------------------------------------------------

 
    
### SEGMENTING SIGNALS
    
    # remove signals that aren't at least 200s  
    sampling_rate <- 1/.02 # may need to adjust the sampling rate if it's not 50 Hz
    df_raw_data %>% group_by(number) %>% count() %>% dplyr::filter(n < 200*sampling_rate) # will output numbers of signals that aren't at least 200s
    df_raw_data <- df_raw_data %>% filter(!(number %in% c())) # in c(), include numbers of signals that aren't long enough
    
    # I have found it best to manually filter signals to remove "blips", or segments where the machine clearly malfunctioned.
    # A function that will automatically do this can be found in /Functions/f_filter_signal.R, 
    # but this function is very slow and not realistic for long signals & many patients
    
    # First, write functions to visualize the signals and easily segment to sections without blips
    
        # plotting functions
    
            # input: df = tibble containing columns number, Time, Hb_tot, and Hb_oxy (i.e. df_raw_data), num = number of patient to plot
            # output: plot of respective signal for the individual patient indicated
        
            # Hb_tot
            f_plot_hbtot <- function(df, num) {
              
              df %>% 
                dplyr::filter(number == num) %>% 
                ggplot(aes(x = Time, y = Hb_tot)) +
                geom_line()
              
            }
            
            # Hb_oxy
            f_plot_hboxy <- function(df, num) {
              
              df %>% 
                dplyr::filter(number == num) %>% 
                ggplot(aes(x = Time, y = Hb_oxy)) +
                geom_line()
              
            }
        
        # segment signal function
            
            # input: df = same df as plotting functions (df_raw_data), num = patient number to segment, 
            #       start = time in seconds indicating beginning of desired signal, end = time inn seconds indicating the end of desired signal
            # output: df containing same columns as input df, but only includes rows within the time frame indicated
            
            f_segment_signal <- function(df, num, start, end) {
          
          df %>% 
            dplyr::filter(number == num) %>% 
            dplyr::filter(Time > start & Time < end) # time is in SECONDS
            }
            
    # Then, use these functions to visually assess the signal and manually determine the longest segment that doesn't contain any "blips"
    # check the variance to make sure it isn't astronomically high (I typically used around 10,000 as a cut-off - 
            # a much higher variance may indicate that the signal needs to be cut shorter or is potentially too noisy to include in analysis)
            
    ### Example
            
        # patient 1
            
            # Hb_tot
            df_raw_data %>% f_plot_hbtot(1) # eyeball where blips are
            df_hbtot1 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(1, start = 0, end = 2000) # shorten signal to desired length, must be > 200s
            df_hbtot1 %>% f_plot_hbtot(1) # reassess visually
            df_hbtot1$Hb_tot %>% var(na.rm = TRUE) # check var less than 10000
            
            # Hb_oxy - repeat for Hb_oxy signal
            df_raw_data %>% f_plot_hboxy(1)
            df_hboxy1 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(1, start = 0, end = 2000)
            df_hboxy1 %>% f_plot_hboxy(1) 
            df_hboxy1$Hb_oxy %>% var(na.rm = TRUE)
            
            
    # Once you have repeated this for all the patients, combine the filtered Hb_tot signals into one tibble and the filtered Hb_oxy signals into another
            
        df_hbtot_filt <- bind_rows(df_hbtot1, ...)
        df_hboxy_filt <- bind_rows(df_hboxy1, ...)
        
    
### CALCULATE AVERAGE HB_TOT AND HB_OXY FOR EACH SIGNAL
        
    
        # Hb_tot
        df_hbtot_mean <- df_hbtot_filt %>% 
          dplyr::filter(!is.na(Hb_tot)) %>% 
          group_by(number, Status) %>% 
          summarise(hbtot_mean = mean(Hb_tot))
        
            # plot
            df_hbtot_mean  %>%  
              ggplot(aes(x = Status, y = hbtot_mean)) +
              geom_violin(aes(color = Status)) +
              geom_jitter(position = position_jitter(0.2), shape = 1) +
              stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
              
              theme_bw() +
              labs(y = "uM") +
              ggtitle("Cerebral tissue hemoglobin concentration") +
              theme(legend.position = "none") 
        
        
        # Hb_oxy
        df_hboxy_mean <- df_hboxy_filt %>% 
          dplyr::filter(!is.na(Hb_oxy)) %>% 
          group_by(number, Status) %>% 
          summarise(hboxy_mean = mean(Hb_oxy))
        
            # plot
            df_hboxy_mean %>%  
              ggplot(aes(x = Status, y = hboxy_mean)) +
              geom_violin(aes(color = Status)) +
              geom_jitter(position = position_jitter(0.2), shape = 1) +
              stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
              
              theme_bw() +
              labs(y = "percent") +
              ggtitle("Cerebral tissue hemoglobin oxygen saturation") +
              theme(legend.position = "none")
        
      # Explore any nonsensical values, especially any negatives in the Hb_oxy. you may want to remove those patients from the analysis due to unreliable readings
    
    
            
### COMPLETE SIGNAL PREPROCESSING IN MATLAB
            
      # add milliseconds in to Time
            
          # Hb_tot
            df_hbtot_filt <- df_hbtot_filt %>% 
              group_by(number, Time) %>% 
              left_join(count(.)) %>% 
              mutate(Time = ifelse(n == 49, Time + 0.02*row_number(), Time + 0.02*(row_number() - 1))) %>% 
              dplyr::select(-n)
            
          # O2 sat
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
            path <- "Data/signal_segments/Hb_tot/"
            1:length(l_hbtot_filt) %>% map(~write.table(l_hbtot_filt[[.x]], file = paste0(path, names(l_hbtot_filt[.x]), ".txt")))
            
          # Hb_oxy
            l_hboxy_filt <- df_hboxy_filt %>% split(df_hboxy_filt$number)
            hboxy_names <- df_hboxy_filt$subject_id %>% unique()
            names(l_hboxy_filt) <- hboxy_names
            path <- "Data/signal_segments/Hb_oxy/"
            1:length(l_hboxy_filt) %>% map(~write.table(l_hboxy_filt[[.x]], file = paste0(path, names(l_hboxy_filt[.x]), ".txt")))
            
            
      # MatLab scripts calculate moving average (across 10 points) and perform bandpass (or highpass) filtering
            
            
            
    
# Detrended fluctuation analysis ------------------------------------------


### LOAD IN DATA
            
      # If not already in your environment, load the processed signals
            
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
            
           
             
### COMPUTE CUMULATIVE SUM OF EACH SIGNAL
            
            
      # DFA is performed on the cumulative sum of the signal, not the signal itself
            
      ## Hb_tot
            
            # create list, where each object is a patient
            l_hbtot_filt_pass <- df_hbtot_filt_pass %>% split(f = df_hbtot_filt_pass$number)
            
            # take the mean of each filtered signal
            l_means <- 1:length(l_hbtot_filt_pass) %>% purrr::map(~mean(l_hbtot_filt_pass[[.x]]$sg_filt, na.rm = TRUE))
            
            # subtract the mean and take the cumulative sum
            l_hbtot_cumsum <- 1:length(l_means) %>% purrr::map(~mutate(l_hbtot_filt_pass[[.x]], mean = l_means[[.x]]))
            df_hbtot_cumsum <- rbindlist(l_hbtot_cumsum) %>% as_tibble()
            
            df_hbtot_cumsum <- df_hbtot_cumsum %>% 
              mutate(subtract_mean = sg_filt - mean) %>% 
              group_by(number) %>% 
              mutate(cumsum = cumsum(subtract_mean)) %>% 
              ungroup()
            
                # plot
                df_hbtot_cumsum %>% dplyr::filter(number == 1) %>% 
                  ggplot(aes(x = Time)) +
                  geom_line(aes(y = cumsum)) +
                  geom_line(aes(y = sg_filt), color = "red") +
                  
                  theme_bw() +
                  ggtitle("Cumulative sum of Savitzy-Golar smoothed signal")
            
            
        ## Hb_oxy
                
            # create list, where each object is a patient
            l_hboxy_filt_pass <- df_hboxy_filt_pass %>% split(f = df_hboxy_filt_pass$number)
            
            # take the mean of each filtered signal
            l_means <- 1:length(l_hboxy_filt_pass) %>% purrr::map(~mean(l_hboxy_filt_pass[[.x]]$sg_filt, na.rm = TRUE))
            
            # subtract the mean and take the cumulative sum
            l_hboxy_cumsum <- 1:length(l_means) %>% purrr::map(~mutate(l_hboxy_filt_pass[[.x]], mean = l_means[[.x]]))
            df_hboxy_cumsum <- rbindlist(l_hboxy_cumsum) %>% as_tibble()
            
            df_hboxy_cumsum <- df_hboxy_cumsum %>% 
              mutate(subtract_mean = sg_filt - mean) %>% 
              group_by(number) %>% 
              mutate(cumsum = cumsum(subtract_mean)) %>% 
              ungroup()
            
                # plot
                df_hboxy_cumsum %>% dplyr::filter(number == 1) %>% 
                  ggplot(aes(x = Time)) +
                  geom_line(aes(y = cumsum)) +
                  geom_line(aes(y = sg_filt), color = "red") +
                  
                  theme_bw() +
                  ggtitle("Cumulative sum of Savitzy-Golar smoothed signal")
            
            

### WRITE DFA FUNCTIONS
                
                
      # write individual DFA functions, then combine into one big function to run at once
                
            
            ## Split into windows
                
                # determine window lengths
                
                # criteria (Hardstone et al 2012):
                # lower end = at least 4 samples (linear detrending will perform poorly with less points)
                # high end = <10% of signal length (anything higher is more noisy bc less than 10 windows to average)
                # need to be equally spaced on a log10 scale (gives equal weight to time scales when fitting  a line)
                # 50% overlap between windows? (used to increase number of windows; Hardstone et al 2012), or non-overlapping (Tagliazucchi et al 2016)
                
                # function: will determine window lengths based on above criteria for a given df
                    # input: df = df_hboxy_cumsum or df_hbtot_cumsum, num_of_windows = number of windows desired (default = 10)
                    # output: vector of numbers indicating window lengths (window_length)
                
                f_window_sizes <- function(df, num_of_windows = 10) {
                  
                  low_end <- 0.1
                  high_end <- floor(nrow(df)/50*0.10) # take 10% of signal length & round down to nearest integer
                  log_size_range <- log10(high_end) - log10(low_end)
                  log_step_size <- log_size_range/num_of_windows # divide by number of windows desired
                  
                  10^(seq(log10(low_end), log10(high_end), by = log_step_size))
                  
                  
                }
                
                # function: split signal by window sizes (window length in s)
                    # input: df = df_hboxy_cumsum or df_hbtot_cumsum, window_length = vector of window sizes (output from above function: window_length)
                    # output: list of lists of tibbles - outer list = window length, inner list = signal broken up into tibbles of that window length (l_splits)
                
                f_split_windows <- function(df, window_length) {
                  
                  window_size <- window_length*50
                  signal_length <- nrow(df)
                  reps <- rep(1:ceiling(signal_length/window_size), each = window_size)[1:signal_length]
                  list <- df %>% split(reps) 
                  1:length(list) %>% purrr::map(~mutate(list[[.x]], window_size = window_length))
                  
                }
                
                
            ## Detrending: remove the linear trend using a least squares fit
                
                # function: compute linear least squares regression line for all windows
                    # input: output from f_split_windows (l_splits)
                    # output: list of lists of linear models for each window (l_lms)
                
                f_lm_least_squares <- function(l_splits) {
                  
                  1:length(l_splits) %>% purrr::map(~lm(cumsum ~ Time, data = l_splits[[.x]]))
                  
                }
                
                # function: detrend by calculating residuals
                    # input: output from f_lm_least_squares (l_lms)
                    # output: list of lists of tibbles containing the residuals at each time point from the linear models (l_resids)
                
                f_calc_residuals <- function(l_splits, l_lms) {
                  
                  1:length(l_lms) %>% purrr::map(~tibble(window_size = l_splits[[.x]]$window_size, 
                                                         Time = l_splits[[.x]]$Time, 
                                                         residuals = l_lms[[.x]]$residuals))  
                  
                }
                
                
            ## Calculate the standard deviation ("fluctuation") of the detrended line
                
                # function: calculate the standard deviation of each line
                    # input: output from f_calc_residuals (l_resids)
                    # output: list of tibbles including window size, number, and fluctuation (l_fluct)
                
                f_calc_fluct <- function(l_resids) {
                  
                  1:length(l_resids) %>% 
                    purrr::map(~tibble(window_size = l_resids[[.x]]$window_size, 
                                       window_num = .x, 
                                       fluctuation = sd(l_resids[[.x]]$residuals))) %>% 
                    rbindlist() %>% 
                    as_tibble()  
                  
                }
                
                # calculate average fluctuation for each window size  
                  # input: output from f_calc_fluct (l_fluct)
                  # output: tibble containing window size with its corresponding average standard deviation (df_avg_fluct)
                
                f_calc_avg_fluct <- function(l_fluct) {
                  
                  1:length(l_fluct) %>% purrr::map(~tibble(window_size = l_fluct[[.x]]$window_size, 
                                                           avg_fluctuation = mean(l_fluct[[.x]]$fluctuation, na.rm = TRUE)) %>% 
                                                     unique()) %>% 
                    rbindlist() %>% as_tibble()
                  
                }
                
              
            ##  Find line of best fit and extract alpha
                
                alpha <- lm(log10(avg_fluctuation) ~ log10(window_size), data = df_avg_fluct) %>% .$coefficients %>% .[2]
                
               
                 
      ### Once you've run/loaded all these functions, combine into one DFA function!
            # input: df = df_hbtot_cumsum or df_hboxy_cumsum, num_of_windows = number of desired windows
            # output: list with two items: the first is a df_avg_fluct, the second is the associated alpha value
                
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
                
                
        # to run for all patients at once: create list of df cumsum, where each object is a patient
            l_hbtot_cumsum <- df_hbtot_cumsum %>% split(f = df_hbtot_cumsum$number)
            l_hboxy_cumsum <- df_hboxy_cumsum %>% split(f = df_hboxy_cumsum$number)
                
                
        # repeat for all patients!!
                
            library(doParallel) # runs faster in parallel
            registerDoParallel()
                
                
            l_dfa_hbtot <- 1:length(l_hbtot_cumsum) %>% purrr::map(~f_dfa(l_hbtot_cumsum[[.x]], num_of_windows = 10))
            l_dfa_hboxy <- 1:length(l_hboxy_cumsum) %>% purrr::map(~f_dfa(l_hboxy_cumsum[[.x]], num_of_windows = 10))
                
                
       # would be wise to save DFA results here because it takes a while to run
            save(l_dfa_hbtot, l_dfa_hboxy, file = "dfa_results.Rdata")   
            
                
                
                
### SECOND ORDER DETRENDING
            
      # in some instances, it is necessary to perform second-order detrending on the log10(avg_fluctuation) vs log10(window_size)
      # because there are differences between the long- vs short-range correlations. (Seleznov et al 2019) 
      # To determine if this is necessary, inspect individual signals. signals that require second-order detrending have a break point or a curve.
                
                
      # Write function to determine the breakpoint in the curve and calculate the short- vs long-range alpha values.
      # If no breakpoint is detected, the function will return the original alpha value (no second-order detrending required)
            
            # load required library
            library(segmented)
            
            # extract the second slope after the "bend" in the alpha slope 
                # input: l_dfa (output from f_dfa function across patients)
                # output: a numeric containing the alpha value representing long-range correlations
            
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
       
            
      # Perform function on DFA output for Hb_tot and Hb_oxy, combine all results into one DF
            
            # Hb_tot
                hbtot_alphas <- 1:length(l_dfa_hbtot) %>% purrr::map(~l_dfa_hbtot[[.x]][[2]]) %>% unlist()
                
                l_hbtot_second_alpha <- 1:length(l_dfa_hbtot) %>% purrr::map(~f_segment_slope(l_dfa_hbtot[[.x]]))
                hbtot_second_alpha <- 1:length(l_dfa_hbtot) %>% purrr::map(~f_segment_slope(l_dfa_hbtot[[.x]])) %>% unlist()
                
                # combine dfa results into df
                df_hbtot_dfa <- df_hbtot_cumsum %>% 
                  count(number, Status) %>% 
                  dplyr::select(-n) %>% 
                  mutate(alpha = hbtot_alphas, 
                         second_alpha = hbtot_second_alpha,
                         Status = as.factor(Status))
            
            # Hb_oxy
                hboxy_alphas <- 1:length(l_dfa_hboxy) %>% purrr::map(~l_dfa_hboxy[[.x]][[2]]) %>% unlist()
                
                l_hboxy_second_alpha <- 1:length(l_dfa_hboxy) %>% purrr::map(~f_segment_slope(l_dfa_hboxy[[.x]]))
                hboxy_second_alpha <- 1:length(l_dfa_hboxy) %>% purrr::map(~f_segment_slope(l_dfa_hboxy[[.x]])) %>% unlist()
                
                # combine dfa results into df
                df_hboxy_dfa <- df_hboxy_cumsum %>% 
                  count(number, Status) %>% 
                  dplyr::select(-n) %>% 
                  mutate(alpha = hboxy_alphas, 
                         second_alpha = hboxy_second_alpha,
                         Status = as.factor(Status))
            
            
                

# plots and statistical analyses ------------------------------------------

                
       # plot alphas for Hb_tot and Hb_oxy
                
                
            # Hb_tot
                
                # plot
                df_hbtot_dfa %>% 
                  ggplot(aes(x = Status, y = second_alpha)) +
                  geom_violin(aes(color = Status)) +
                  geom_jitter(position = position_jitter(0.2), shape = 1) +
                  stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
                  
                  labs(y = "alpha") +
                  ggtitle("Cerebral tissue hemoglobin concentration alpha") +
                  
                  theme_bw() +
                  theme(legend.position = "none")
                
                # test for differences
                df_hbtot_dfa %>% group_by(Status) %>% summarise(median = median(second_alpha))
                kruskal.test(second_alpha ~ Status, data = df_hbtot_dfa)
                DunnTest(second_alpha ~ Status, data = df_hbtot_dfa)
                
                
                
             # Hb_oxy
                
                # plot
                df_hboxy_dfa %>% 
                  ggplot(aes(x = Status, y = second_alpha)) +
                  geom_violin(aes(color = Status)) +
                  geom_jitter(position = position_jitter(0.2), shape = 1) +
                  stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
                  
                  labs(y = "alpha") +
                  ggtitle("Cerebral tissue hemoglobin oxygen saturation alpha") +
                  
                  theme_bw() +
                  theme(legend.position = "none")
                
                
                # test for differences
                df_hboxy_dfa %>% group_by(Status) %>% summarise(median = median(second_alpha))
                kruskal.test(second_alpha ~ Status, data = df_hboxy_dfa)
                DunnTest(second_alpha ~ Status, data = df_hboxy_dfa)
                
                
                
                
                
                
                
            
            
            
            
            
            
            
  
  
  
  
  
  
  
  