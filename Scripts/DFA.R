



# load in processed signal data -------------------------------------------

  load("filtered_signals.Rdata")





# calculate amplitude envelope using Hilbert transform --------------------

  library(seewave)
  library(zoo)
  library(DescTools)

  # add milliseconds in to Time
  df_thc_filt_pass <- df_thc_filt_pass %>% 
    group_by(number, Time) %>% 
    left_join(count(.)) %>% 
    mutate(Time = ifelse(n == 49, Time + 0.02*row_number(), Time + 0.02*(row_number() - 1)))
  
  # 1 second rolling mean?
  df_thc_filt_pass <- df_thc_filt_pass %>% group_by(number) %>% 
    mutate(roll_mean = zoo::rollmean(THC, k = 50, fill = NA)) %>% 
    ungroup()

  # or we can just smooth the amplitude envelope
  hil_amp <- env(df_thc_filt_pass1$THC, f = 50, msmooth = c(50, 0), plot = FALSE) %>% 
    as_tibble(.name_repair = "universal") %>% 
    dplyr::rename("hil_amp" = "...1") %>% 
    slice(rep(1:n(), each = 50)) # smoothed hilbert transform 1 sec
  
  # check with plotting
  df_thc_filt_pass1[1:58450,] %>% bind_cols(hil_amp) %>% 
    
    ggplot(aes(x = Time)) +
    geom_line(aes(y = THC)) +
    geom_line(aes(y = hil_amp), color = "red")
  
  # use the hilbert amplitude envelope for all signals
  l_hil_amp <- unique(df_thc_filt_pass$number) %>% 
    purrr::map(~dplyr::filter(df_thc_filt_pass, number == .x, !is.na(THC))) %>% 
    purrr::map(~env(.x$THC, f = 50, msmooth = c(50, 0), plot = FALSE) %>% 
                         as_tibble(.name_repair = "universal") %>% 
                         dplyr::rename("hil_amp" = "...1") %>% 
                         slice(rep(1:n(), each = 50))) # nrow rounds down to the nearest 50
  names(l_hil_amp) <-  unique(df_thc_filt_pass$number)

  # combine with df_thc_filt_pass
  l_signal_abbr <- unique(df_thc_filt_pass$number) %>% 
    purrr::map(~dplyr::filter(df_thc_filt_pass, number == .x)) %>% 
    purrr::map(~.x[1:RoundTo(nrow(.x), multiple = 50, FUN = "floor"),])
  
  l_signal_hil_amp <- 1:length(l_signal_abbr) %>% purrr::map(~bind_cols(l_signal_abbr[[.x]], l_hil_amp[[.x]]))
  
  # create df
  df_thc_hil_amp <- l_signal_hil_amp %>% rbindlist() %>% as_tibble()
  
  
  

    

# compute cumulative sum --------------------------------------------------

    
    # take the mean of each signal
    l_means <- 1:length(l_signal_hil_amp) %>% purrr::map(~mean(l_signal_hil_amp[[.x]]$roll_mean, na.rm = TRUE))

    # subtract the mean and take the cumulative sum
    l_thc_cumsum <- 1:length(l_signal_hil_amp) %>% purrr::map(~filter(l_signal_hil_amp[[.x]], !is.na(roll_mean)))
    
    l_thc_cumsum <- 1:length(l_means) %>% purrr::map(~mutate(l_signal_hil_amp[[.x]], mean = l_means[[.x]]))
    
    df_thc_cumsum <- rbindlist(l_thc_cumsum) %>% as_tibble() %>% filter(!is.na(roll_mean))
    
    df_thc_cumsum <- df_thc_cumsum %>% 
      mutate(subtract_mean = roll_mean - mean) %>% 
      group_by(number) %>% 
      mutate(cumsum = cumsum(subtract_mean)) %>% 
      ungroup()
    
    # plot
    df_thc_cumsum %>% filter(number == 7) %>% 
      
      ggplot(aes(x = Time)) +
      geom_line(aes(y = roll_mean)) +
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
    window_lengths <- c(0.1, 1, 10, 100)
    
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
    f_plot_window <- function(l_split, window_num) {
      
      l_split[[window_num]] %>% 
        ggplot(aes(x = Time)) +
        geom_line(aes(y = THC), color = "blue") +
        geom_line(aes(y = cumsum), color = "black") +
        geom_smooth(aes(y = cumsum), method = "lm", color = "red") +
        theme_bw()
      
    }
    
    f_plot_window(l_splits[[1]], 1)
    
    
    # plot detrended window
    p_plot_resid <- function(l_resid, window_num) {
      
      l_resid[[window_num]] %>% 
        ggplot(aes(x = Time, y = residuals)) +
        geom_line() +
        geom_hline(yintercept = 0, color = "red") +
        theme_bw()
      

    }
    
    p_plot_resid(l_resids[[1]], 1)
    
    
    

# calculate standard deviation of detrended line --------------------------

    
    # write function to return window size, number, and fluctuation

    f_calc_fluct <- function(l_resids) {
      
      1:length(l_resids) %>% 
        purrr::map(~tibble(window_size = l_resids[[.x]]$window_size, window_num = .x, fluctuation = sd(l_resids[[.x]]$residuals))) %>% 
        rbindlist() %>% 
        as_tibble()  
      
    }

    l_fluct <- 1:length(l_resids) %>% purrr::map(~f_calc_fluct(l_resids[[.x]]) %>% unique())
    

    # calculate average fluctuation for each window size  
    
    f_calc_avg_fluct <- function(l_fluct) {
      
      1:length(l_fluct) %>% purrr::map(~tibble(window_size = l_fluct[[.x]]$window_size, avg_fluctuation = mean(l_fluct[[.x]]$fluctuation, na.rm = TRUE)) %>% unique()) %>% 
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
      window_lengths <- f_window_sizes(df_thc_cumsum, num_of_windows)
      
      # split df into a list of lists, where length of the overall list is the number of window lengths, and each list item contains windows of that length
      l_splits <- window_lengths %>% purrr::map(~f_split_windows(df_thc_cumsum, .x))
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
    
    
    dfa1 <- f_dfa(df_thc_cumsum1, num_of_windows = 10)
    
    # plot
    dfa1[[1]] %>% 
      ggplot(aes(x = log10(window_size), y = log10(avg_fluctuation))) +
      geom_point() +
      geom_smooth(method = "lm")
    
    
    
    l_dfa <- unique(df_thc_cumsum$number) %>% purrr::map(~f_dfa(df_thc_cumsum %>% filter(number == .x), num_of_windows = 10))
    
    
    l_thc_cumsum <- df_thc_cumsum %>% split(f = df_thc_cumsum$number)
    
    length(l_thc_cumsum)

    dfa2 <- l_thc_cumsum[[2]] %>% f_dfa(num_of_windows = 10)
    
    
    
    
    
    
    
    
    