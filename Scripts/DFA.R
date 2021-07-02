



# load in processed signal data -------------------------------------------

  load("filtered_signals.Rdata")





# calculate amplitude envelope using Hilbert transform --------------------

  # library(seewave)
  # library(zoo)
  # library(DescTools)

  # # smooth the amplitude envelope
  # hil_amp <- env(df_thc_filt_pass1$THC, f = 50, msmooth = c(50, 0), plot = FALSE) %>% 
  #   as_tibble(.name_repair = "universal") %>% 
  #   dplyr::rename("hil_amp" = "...1") %>% 
  #   slice(rep(1:n(), each = 50)) # smoothed hilbert transform 1 sec
  # 
  # # check with plotting
  # df_thc_filt_pass1[1:58450,] %>% bind_cols(hil_amp) %>% 
  #   
  #   ggplot(aes(x = Time)) +
  #   geom_line(aes(y = THC)) +
  #   geom_line(aes(y = hil_amp), color = "red")
  # 
  # # use the hilbert amplitude envelope for all signals
  # l_hil_amp <- unique(df_thc_filt_pass$number) %>% 
  #   purrr::map(~dplyr::filter(df_thc_filt_pass, number == .x, !is.na(THC))) %>% 
  #   purrr::map(~env(.x$THC, f = 50, msmooth = c(50, 0), plot = FALSE) %>% 
  #                        as_tibble(.name_repair = "universal") %>% 
  #                        dplyr::rename("hil_amp" = "...1") %>% 
  #                        slice(rep(1:n(), each = 50))) # nrow rounds down to the nearest 50
  # names(l_hil_amp) <-  unique(df_thc_filt_pass$number)
  # 
  # # combine with df_thc_filt_pass
  # l_signal_abbr <- unique(df_thc_filt_pass$number) %>% 
  #   purrr::map(~dplyr::filter(df_thc_filt_pass, number == .x)) %>% 
  #   purrr::map(~.x[1:RoundTo(nrow(.x), multiple = 50, FUN = "floor"),])
  # 
  # l_signal_hil_amp <- 1:length(l_signal_abbr) %>% purrr::map(~bind_cols(l_signal_abbr[[.x]], l_hil_amp[[.x]]))
  # 
  # # create df
  # df_thc_hil_amp <- l_signal_hil_amp %>% rbindlist() %>% as_tibble()
  
  
  

    

# compute cumulative sum --------------------------------------------------

    # create list, where each object is a patient
    l_thc_filt_pass <- df_thc_filt_pass %>% split(f = df_thc_filt_pass$number)

    # # take the mean of each signal
    # l_means <- 1:length(l_thc_filt_pass) %>% purrr::map(~mean(l_thc_filt_pass[[.x]]$roll_mean, na.rm = TRUE))
    # 
    # # subtract the mean and take the cumulative sum
    # l_thc_cumsum <- 1:length(l_thc_filt_pass) %>% purrr::map(~filter(l_thc_filt_pass[[.x]], !is.na(roll_mean)))
    # l_thc_cumsum <- 1:length(l_means) %>% purrr::map(~mutate(l_thc_filt_pass[[.x]], mean = l_means[[.x]]))
    # df_thc_cumsum <- rbindlist(l_thc_cumsum) %>% as_tibble()
    # 
    # df_thc_cumsum <- df_thc_cumsum %>% 
    #   mutate(subtract_mean = roll_mean - mean) %>% 
    #   group_by(number) %>% 
    #   mutate(cumsum = cumsum(subtract_mean)) %>% 
    #   ungroup()
    
    # take the mean of each filtered signal
    l_means <- 1:length(l_thc_filt_pass) %>% purrr::map(~mean(l_thc_filt_pass[[.x]]$band_pass, na.rm = TRUE))
    
    # subtract the mean and take the cumulative sum
    l_thc_cumsum <- 1:length(l_thc_filt_pass) %>% purrr::map(~filter(l_thc_filt_pass[[.x]], !is.na(band_pass)))
    l_thc_cumsum <- 1:length(l_means) %>% purrr::map(~mutate(l_thc_filt_pass[[.x]], mean = l_means[[.x]]))
    df_thc_cumsum <- rbindlist(l_thc_cumsum) %>% as_tibble()
    
    df_thc_cumsum <- df_thc_cumsum %>% 
      mutate(subtract_mean = band_pass - mean) %>% 
      group_by(number) %>% 
      mutate(cumsum = cumsum(subtract_mean)) %>% 
      ungroup()
    
    # plot
    df_thc_cumsum %>% filter(number == 1) %>% 
      
      ggplot(aes(x = Time)) +
      geom_line(aes(y = band_pass)) +
      geom_line(aes(y = cumsum), color = "red")
    
    
    
    
# split into windows --------------------------------------------------------------------
    
    
    
    # start with signal from patient 1
    df_thc_cumsum1 <- df_thc_cumsum %>% filter(number == 1)
    df_thc_cumsum1 %>% nrow()/50 # signal is 1169 seconds long
    
    # determine window lengths
    
      # criteria:
      # lower end = at least 4 samples (linear detrending will perform poorly with less points)
      # high end = <10% of signal length (anything higher is more noisy bc less than 10 windows to average)
      # need to be equally spaced on a log10 scale (gives equal weight to time scales when fitting  a line)
      # 50% overlap between windows? (used to increase number of windows; Hardstone et al 2012), or non-overlapping (Tagliazucchi et al 2016)
    
    # function: will determine window lengths based on above criteria for a given df
    f_window_sizes <- function(df, num_of_windows) {
      
      low_end <- 0.1
      high_end <- floor(nrow(df)/50*0.10) # take 10% of signal length & round down to nearest integer
      log_size_range <- log10(high_end) - log10(low_end)
      log_step_size <- log_size_range/num_of_windows # divide by number of windows desired
      
      10^(seq(log10(low_end), log10(high_end), by = log_step_size))
      

    }
    
    f_window_sizes(df_thc_cumsum1, 10)
    
    # function to split by window sizes (window length in s)
    f_split_windows <- function(df, window_length) {
      
      window_size <- window_length*50
      signal_length <- nrow(df)
      reps <- rep(1:ceiling(signal_length/window_size), each = window_size)[1:signal_length]
      list <- df %>% split(reps) 
      1:length(list) %>% purrr::map(~mutate(list[[.x]], window_size = window_length))

    }
    
    # start simple with 4 window sizes - 0.1 s, 1 s, 10s, and 100s
    window_lengths <- f_window_sizes(df_thc_cumsum1, 10)
    
    l_splits <- window_lengths %>% purrr::map(~f_split_windows(df_thc_cumsum1, .x))
    names(l_splits) <- window_lengths
    

    
    
# remove the linear trend using a least squares fit -----------------------


    # compute linear least squares regression line for all windows

    f_lm_least_squares <- function(l_splits) {
      
      1:length(l_splits) %>% purrr::map(~lm(cumsum ~ Time, data = l_splits[[.x]]))
      
    }
    
    l_lms <- 1:length(l_splits) %>% purrr::map(~f_lm_least_squares(l_splits[[.x]]))
    names(l_lms) <- window_lengths
    
    
    # detrend by calculating residuals
    
    f_calc_residuals <- function(l_splits, l_lms) {
      
      1:length(l_lms) %>% purrr::map(~tibble(window_size = l_splits[[.x]]$window_size, 
                                             Time = l_splits[[.x]]$Time, 
                                             residuals = l_lms[[.x]]$residuals))  
      
    }
    
    l_resids <- 1:length(l_splits) %>% purrr::map(~f_calc_residuals(l_splits[[.x]], l_lms[[.x]]))
    names(l_resids) <- window_lengths
    
    
    # plot window function
    f_plot_window <- function(l_splits, window_num) {
      
      l_splits[[window_num]] %>% 
        ggplot(aes(x = Time)) +
        geom_line(aes(y = THC), color = "blue") +
        geom_line(aes(y = cumsum), color = "black") +
        geom_smooth(aes(y = cumsum), method = "lm", color = "red") +
        theme_bw()
      
    }
    
    f_plot_window(l_splits[[6]], 1)
    
    
    # plot detrended window
    p_plot_resids <- function(l_resids, window_num) {
      
      l_resids[[window_num]] %>% 
        ggplot(aes(x = Time, y = residuals)) +
        geom_line() +
        geom_hline(yintercept = 0, color = "red") +
        theme_bw()
      

    }
    
    p_plot_resids(l_resids[[6]], 1)
    
    
    

# calculate standard deviation of detrended line --------------------------

    
    # write function to return window size, number, and fluctuation

    f_calc_fluct <- function(l_resids) {
      
      1:length(l_resids) %>% 
        purrr::map(~tibble(window_size = l_resids[[.x]]$window_size, 
                           window_num = .x, 
                           fluctuation = sd(l_resids[[.x]]$residuals))) %>% 
        rbindlist() %>% 
        as_tibble()  
      
    }

    l_fluct <- 1:length(l_resids) %>% purrr::map(~f_calc_fluct(l_resids[[.x]]) %>% unique())
    

    # calculate average fluctuation for each window size  
    
    f_calc_avg_fluct <- function(l_fluct) {
      
      1:length(l_fluct) %>% purrr::map(~tibble(window_size = l_fluct[[.x]]$window_size, 
                                               avg_fluctuation = mean(l_fluct[[.x]]$fluctuation, na.rm = TRUE)) %>% 
                                         unique()) %>% 
        rbindlist() %>% as_tibble()
      
    }
    
    df_avg_fluct <- f_calc_avg_fluct(l_fluct)
    
    
    # find line of best fit and extract alpha
    alpha <- lm(log10(avg_fluctuation) ~ log10(window_size), data = df_avg_fluct) %>% .$coefficients %>% .[2]
    
    
    # plot
    df_avg_fluct %>% 
      ggplot(aes(x = log10(window_size), y = log10(avg_fluctuation))) +
      geom_point() +
      geom_smooth(method = "lm")
    
    
    
    
# DFA function!! ----------------------------------------------------------

    

    
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
    
    
    # create list of df cumsum, where each object is a patient
    l_thc_cumsum <- df_thc_cumsum %>% split(f = df_thc_cumsum$number)
    
    # example
    dfa1 <- f_dfa(l_thc_cumsum[[1]], num_of_windows = 10)
    dfa2 <- f_dfa(l_thc_cumsum[[2]], num_of_windows = 10)
    dfa99 <- f_dfa(l_thc_cumsum[[93]], num_of_windows = 10)
    dfa98 <- f_dfa(l_thc_cumsum[[92]], num_of_windows = 10)
    
    # repeat for all patients!!
    
    library(doParallel)
    registerDoParallel()
    
      # on rolling mean, not filtered signal
        l_dfa <- 1:length(l_thc_cumsum) %>% purrr::map(~f_dfa(l_thc_cumsum[[.x]], num_of_windows = 10))
        
        # combine into df
        df_alphas <- 1:length(l_dfa) %>% 
          purrr::map(~tibble(alpha = l_dfa[[.x]][[2]])) %>% 
          rbindlist() %>% 
          as_tibble() %>% 
          bind_cols(df_thc_cumsum %>% count(number, Status) %>% dplyr::select(c(-n))) %>% 
          mutate(Status = as.factor(Status), Hurst = ifelse(alpha > 1, alpha - 1, alpha)) 
        
        # plot
        df_alphas %>% 
          ggplot(aes(x = Status, y = alpha)) +
          geom_violin(aes(color = Status)) +
          geom_jitter(position = position_jitter(0.2), shape = 1) +
          stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
          theme_bw()
        
        # test for differences
        kruskal.test(alpha ~ Status, data = df_alphas)
        DunnTest(alpha ~ Status, data = df_alphas)
        
      # on filtered signal, not rolling mean
        l_dfa_filt <- 1:length(l_thc_cumsum) %>% purrr::map(~f_dfa(l_thc_cumsum[[.x]], num_of_windows = 10))
        
        # combine into df
        df_alphas_filt <- 1:length(l_dfa_filt) %>% 
          purrr::map(~tibble(alpha = l_dfa_filt[[.x]][[2]])) %>% 
          rbindlist() %>% 
          as_tibble() %>% 
          bind_cols(df_thc_cumsum %>% count(number, Status) %>% dplyr::select(c(-n))) %>% 
          mutate(Status = as.factor(Status), Hurst = ifelse(alpha > 1, alpha - 1, alpha)) 
        
        # plot
        df_alphas_filt %>% 
          ggplot(aes(x = Status, y = alpha)) +
          geom_violin(aes(color = Status)) +
          geom_jitter(position = position_jitter(0.2), shape = 1) +
          stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
          theme_bw()
        
        # test for differences
        kruskal.test(alpha ~ Status, data = df_alphas_filt)
        DunnTest(alpha ~ Status, data = df_alphas_filt)
        

    
    

# are signals segmented?? --------------------------------------------------


        library(segmented)
        
        # band pass
        df_dfa1 <- dfa1[[1]] %>% mutate(log10_avg_fluct = log10(avg_fluctuation), log10_window_size = log10(window_size))
        fit1 <- lm(log10_avg_fluct ~ log10_window_size, data = df_dfa1) 
        segfit1 <- segmented(fit1)
        
        segfit1 %>% summary()
        segfit1 %>% slope()
        segfit1 %>% fitted()
        
        df_dfa1 %>% 
          ggplot(aes(x = log10_window_size, y = log10_avg_fluct)) +
          geom_point() +
          geom_line(aes(y = segfit1 %>% fitted()), color = "red")
        
        # roll mean
          # 1
          df_dfa1 <- l_dfa[[1]][[1]] %>% mutate(log10_avg_fluct = log10(avg_fluctuation), log10_window_size = log10(window_size))
          fit1 <- lm(log10_avg_fluct ~ log10_window_size, data = df_dfa1) 
          segfit1 <- segmented(fit1)
          
          segfit1 %>% summary()
          segfit1 %>% slope()
          segfit1 %>% fitted()
  
          df_dfa1 %>% 
            ggplot(aes(x = log10_window_size, y = log10_avg_fluct)) +
            geom_point() +
            geom_line(aes(y = segfit1 %>% fitted()), color = "red")
          
          # 99
          df_dfa99 <- l_dfa[[93]][[1]] %>% mutate(log10_avg_fluct = log10(avg_fluctuation), log10_window_size = log10(window_size))
          fit99 <- lm(log10_avg_fluct ~ log10_window_size, data = df_dfa99) 
          segfit99 <- segmented(fit99)
          
          segfit99 %>% summary()
          segfit99 %>% slope()
          segfit99 %>% fitted()
          
          df_dfa99 %>% 
            ggplot(aes(x = log10_window_size, y = log10_avg_fluct)) +
            geom_point() +
            geom_line(aes(y = segfit99 %>% fitted()), color = "red") ## maybe??
          
        # write function to extract second slope 
          
          f_segment_slope <- function(l_dfa) {
            
            df <- l_dfa[[1]] %>% mutate(log10_avg_fluct = log10(avg_fluctuation), log10_window_size = log10(window_size))
            fit <- lm(log10_avg_fluct ~ log10_window_size, data = df) 
            segfit <- segmented(fit)
            
            slope(segfit)$log10_window_size[2,1]
            
          }
          
          l_second_alpha <- 1:length(l_dfa) %>% purrr::map(~f_segment_slope(l_dfa[[.x]]))
          
          # plot
          df_second_alphas <- df_alphas %>% mutate(second_alpha = l_second_alpha %>% unlist()) 
            
          df_second_alphas %>% 
            ggplot(aes(x = Status, y = second_alpha)) +
            geom_violin(aes(color = Status)) +
            geom_jitter(position = position_jitter(0.2), shape = 1) +
            stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
            theme_bw()
          
          # test for differences
          kruskal.test(second_alpha ~ Status, data = df_second_alphas)
          DunnTest(second_alpha ~ Status, data = df_second_alphas)
          
          
          

# try savitzky golay & new band pass --------------------------------------

          # take the mean of each filtered signal
          l_means <- 1:length(l_thc_filt_pass) %>% purrr::map(~mean(l_thc_filt_pass[[.x]]$band_pass, na.rm = TRUE))
          
          # subtract the mean and take the cumulative sum
          l_thc_cumsum <- 1:length(l_thc_filt_pass) %>% purrr::map(~filter(l_thc_filt_pass[[.x]], !is.na(band_pass)))
          l_thc_cumsum <- 1:length(l_means) %>% purrr::map(~mutate(l_thc_filt_pass[[.x]], mean = l_means[[.x]]))
          df_thc_cumsum <- rbindlist(l_thc_cumsum) %>% as_tibble()
          
          df_thc_cumsum <- df_thc_cumsum %>% 
            mutate(subtract_mean = band_pass - mean) %>% 
            group_by(number) %>% 
            mutate(cumsum = cumsum(subtract_mean)) %>% 
            ungroup()
          
          # create list of df cumsum, where each object is a patient
          l_thc_cumsum <- df_thc_cumsum %>% split(f = df_thc_cumsum$number)
          
          # examples
          df_thc_cumsum %>% count(number, Status) %>% filter(Status == "HV")
            # CM
            dfa1 <- f_dfa(l_thc_cumsum[[1]], num_of_windows = 10)
            dfa15 <- f_dfa(l_thc_cumsum[[13]], num_of_windows = 10)
            
            # UM
            dfa58 <- f_dfa(l_thc_cumsum[[54]], num_of_windows = 10)
            
            # HV
            dfa99 <- f_dfa(l_thc_cumsum[[93]], num_of_windows = 10)
            dfa86 <- f_dfa(l_thc_cumsum[[80]], num_of_windows = 10)
            
            
            
            # take the mean of each filtered signal
            l_means <- 1:length(l_thc_filt_pass) %>% purrr::map(~mean(l_thc_filt_pass[[.x]]$sg_filt, na.rm = TRUE))
            
            # subtract the mean and take the cumulative sum
            l_thc_cumsum <- 1:length(l_thc_filt_pass) %>% purrr::map(~filter(l_thc_filt_pass[[.x]], !is.na(sg_filt)))
            l_thc_cumsum <- 1:length(l_means) %>% purrr::map(~mutate(l_thc_filt_pass[[.x]], mean = l_means[[.x]]))
            df_thc_cumsum <- rbindlist(l_thc_cumsum) %>% as_tibble()
            
            df_thc_cumsum <- df_thc_cumsum %>% 
              mutate(subtract_mean = sg_filt - mean) %>% 
              group_by(number) %>% 
              mutate(cumsum = cumsum(subtract_mean)) %>% 
              ungroup()
            
            # create list of df cumsum, where each object is a patient
            l_thc_cumsum <- df_thc_cumsum %>% split(f = df_thc_cumsum$number)
            
            # examples
            # CM
            dfa1 <- f_dfa(l_thc_cumsum[[1]], num_of_windows = 10)
            dfa15 <- f_dfa(l_thc_cumsum[[13]], num_of_windows = 10)
            dfa26 <- f_dfa(l_thc_cumsum[[23]], num_of_windows = 10)
            dfa42 <- f_dfa(l_thc_cumsum[[37]], num_of_windows = 10)
            dfa45 <- f_dfa(l_thc_cumsum[[40]], num_of_windows = 10)
            
            # UM
            dfa55 <- f_dfa(l_thc_cumsum[[50]], num_of_windows = 10)
            dfa58 <- f_dfa(l_thc_cumsum[[54]], num_of_windows = 10)
            dfa62 <- f_dfa(l_thc_cumsum[[57]], num_of_windows = 10)
            dfa65 <- f_dfa(l_thc_cumsum[[60]], num_of_windows = 10)
            dfa71 <- f_dfa(l_thc_cumsum[[65]], num_of_windows = 10)

            # HV            
            dfa73 <- f_dfa(l_thc_cumsum[[67]], num_of_windows = 10)
            dfa77 <- f_dfa(l_thc_cumsum[[71]], num_of_windows = 10)
            dfa81 <- f_dfa(l_thc_cumsum[[75]], num_of_windows = 10)
            dfa86 <- f_dfa(l_thc_cumsum[[80]], num_of_windows = 10)
            dfa99 <- f_dfa(l_thc_cumsum[[93]], num_of_windows = 10)
            
          
            dfa1[[1]] %>% mutate(log10_avg_fluct = log10(avg_fluctuation), log10_window_size = log10(window_size)) %>% 
              
              ggplot(aes(x = log10_window_size, y = log10_avg_fluct)) +
              geom_point() +
              geom_smooth(method = "lm") +
              theme_bw()

            # second alphas?
            f_segment_slope(dfa99)
            l_examples <- list(dfa1, dfa15, dfa26, dfa42, dfa45,
                               dfa55, dfa58, dfa62, dfa65, dfa71,
                               dfa73, dfa77, dfa81, dfa86, dfa99)
          
          alpha2 <- 1:length(l_examples) %>% 
            purrr::map(~f_segment_slope(l_examples[[.x]])) %>% 
            unlist()
   
         df_alpha2_ex <- tibble(Status = as.factor(c(rep("CM", 5), rep("UM", 5), rep("HV", 5))), second_alpha = alpha2)
         df_alpha2_ex %>% 
           ggplot(aes(x = Status, y = second_alpha)) +
           geom_violin(aes(color = Status)) +
           geom_jitter(position = position_jitter(0.2), shape = 1) +
           stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
           theme_bw()
         
         # test for differences
         kruskal.test(second_alpha ~ Status, data = df_alpha2_ex)
         library(DescTools)
         DunnTest(second_alpha ~ Status, data = df_alpha2_ex)
   
   
         
         ### perform on all patients!!
         1:length(l_thc_cumsum) %>% purrr::map(f_dfa(~l_thc_cumsum[[.x]], num_of_windows = 10))
   
   
   
   
   
   
   
    
    
    