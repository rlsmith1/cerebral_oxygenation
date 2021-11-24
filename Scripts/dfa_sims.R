


# libraries ---------------------------------------------------------------

    library(tidyverse)
    library(purrr)
    library(data.table)



# split into windows --------------------------------------------------------------------



  # determine window lengths
  
  
`  # function: will determine window lengths based on above criteria for a given df
  
  sampling_freq <- 50
  
  f_window_sizes <- function(df, num_of_windows) {
    
    low_end <- 0.1 # need at least 4 points to regress
    high_end <- floor(nrow(df)/sampling_freq*0.10) # take 10% of signal length & round down to nearest integer
    log_size_range <- log10(high_end) - log10(low_end)
    log_step_size <- log_size_range/num_of_windows # divide by number of windows desired
    
    10^(seq(log10(low_end), log10(high_end), by = log_step_size))
    
    
  }

  # function to split by window sizes (window length in s)
  f_split_windows <- function(df, window_length) {
    
    window_size <- window_length*sampling_freq
    signal_length <- nrow(df)
    reps <- rep(1:ceiling(signal_length/window_size), each = window_size)[1:signal_length]
    list <- df %>% split(reps) 
    1:length(list) %>% purrr::map(~mutate(list[[.x]], window_size = window_length))
    
  }





# remove the linear trend using a least squares fit -----------------------


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
  
  
  # calculate average fluctuation for each window size  
  
  f_calc_avg_fluct <- function(l_fluct) {
    
    1:length(l_fluct) %>% purrr::map(~tibble(window_size = l_fluct[[.x]]$window_size, 
                                             avg_fluctuation = mean(l_fluct[[.x]]$fluctuation, na.rm = TRUE)) %>% 
                                       unique()) %>% 
      rbindlist() %>% as_tibble()
    
  }
  
  
  # find line of best fit and extract alpha
  alpha <- lm(log10(avg_fluctuation) ~ log10(window_size), data = df_avg_fluct) %>% .$coefficients %>% .[2]
  
  # plot
  df_avg_fluct %>% 
    ggplot(aes(x = log10(window_size), y = log10(avg_fluctuation))) +
    geom_point() +
    geom_smooth(method = "lm")

  
  
# full DFA function!! ----------------------------------------------------------
  
  
  
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
  

  
# generate signals ---------------------------------------------------------
  
  
  # white noise
  
    white_noise <- arima.sim(model = list(order = c(0, 0, 0)), n = 1000) %>% 
      as_tibble_col(column_name = "signal") %>% 
      mutate(Time = 1:nrow(.))
    
    # take the mean
    mean <- mean(white_noise$signal)
    
    # subtract the mean and take the cumulative sum
    df_white_noise <- white_noise %>% 
      mutate(subtract_mean = signal - mean) %>% 
      mutate(cumsum = cumsum(subtract_mean))
    
    # plot
    df_white_noise %>% 
      ggplot(aes(x = Time)) +
      geom_line(aes(y = cumsum)) +
      geom_line(aes(y = signal), color = "red") +
      theme_bw()
    
    
  # pink noise
    require(R.matlab)
    pinknoise <- readMat('Data/pinknoise.mat')
    
    require(matlabr)
    run_matlab_code("pinknoise(1000);", paths_to_add = "~/")
    
    setwd("~")
    setwd("/Applications/MATLAB_R2013a.app/bin/matlab")
    system('matlab -nodisplay -r "a=2; b=1; display(a+b); exit"')
    
    library(R.matlab)
    Matlab$startServer()
    matlab <- Matlab()
    isOpen <- open(matlab)
    
    pink_noise <- c(pinknoise$pinknoise) %>% 
      as_tibble_col(column_name = "signal") %>% 
      mutate(Time = 1:nrow(.))
    
    # take the mean
    mean <- mean(pink_noise$signal)
    
    # subtract the mean and take the cumulative sum
    df_pink_noise <- pink_noise %>% 
      mutate(subtract_mean = signal - mean) %>% 
      mutate(cumsum = cumsum(subtract_mean))
  
    df_pink_noise %>% f_dfa(10)

  # brown noise
    
    brown_noise <- arima.sim(model = list(order = c(0, 1, 0)), n = 1000) %>% 
      as_tibble_col(column_name = "signal") %>% 
      mutate(Time = 1:nrow(.))
    
    # take the mean
    mean <- mean(brown_noise$signal)
    
    # subtract the mean and take the cumulative sum
    df_brown_noise <- brown_noise %>% 
      mutate(subtract_mean = signal - mean) %>% 
      mutate(cumsum = cumsum(subtract_mean))
    
  

# run function on white,  pink,  and brown noise --------------------------

  

  # noise simulations function
  
  noise_simulations <- function(noise, n_sims = 100) {
    
    noise_model <- case_when(
      
      noise == "white" ~ list(order = c(0, 0, 0)),
      noise == "brown" ~ list(order = c(0, 1, 0)),
      noise == "pink" ~ list(order = c(1, 0, 1))
      
      )
    
    names(noise_model) <- "order"
    
    l_noise <- list()
    for (i in 1:n_sims) {
      
      # generate i instances of white noise
      tmp <- arima.sim(model = noise_model, n = 1000) %>% 
        as_tibble_col(column_name = "signal") %>% 
        mutate(Time = 1:nrow(.))
      
      l_noise[[i]] <- tmp
      
    }
    
    # take the mean
    l_means <- 1:length(l_noise) %>% purrr::map(~mean(l_noise[[.x]]$signal))
    
    # subtract mean and take the cumulative sum
    l_cumsum <- 1:length(l_means) %>% purrr::map(~mutate(l_noise[[.x]], mean = l_means[[.x]], number = .x))
    l_cumsum <- 1:length(l_cumsum) %>% purrr::map(~mutate(l_cumsum[[.x]], subtract_mean = signal - mean))
    l_cumsum <- 1:length(l_cumsum) %>% purrr::map(~mutate(l_cumsum[[.x]], cumsum = cumsum(subtract_mean)))
    

    return(l_cumsum)

  }
  
  # run simulations
  l_white_noise <- noise_simulations("white")
  l_brown_noise <- noise_simulations("brown")
  
  # run DFA on simulations
  library(doParallel)
  registerDoParallel()
  
  l_wn_dfa <- 1:length(l_white_noise) %>% purrr::map(~f_dfa(l_white_noise[[.x]],10))
  l_bn_dfa <- 1:length(l_brown_noise) %>% purrr::map(~f_dfa(l_brown_noise[[.x]],10))
  
  # plot distribution of simulated alphas
  
    # combine noise alphas
    df_wn_alpha <- 1:length(l_wn_dfa) %>% 
      purrr::map(~l_wn_dfa[[.x]][2]) %>% 
      unlist() %>% 
      as_tibble_col(column_name = "alpha") %>% 
      mutate(noise = "white")
    
    df_bn_alpha <- 1:length(l_bn_dfa) %>% 
      purrr::map(~l_bn_dfa[[.x]][2]) %>% 
      unlist() %>% 
      as_tibble_col(column_name = "alpha") %>% 
      mutate(noise = "brown")
    
    # calculate 95% CIs
    df_alpha <- df_bn_alpha %>% 
      bind_rows(df_wn_alpha) %>% 
      
      group_by(noise) %>% 
      mutate(median = median(alpha),
             stdev = sd(alpha),
             error = qnorm(0.975)*stdev/sqrt(100),
             low = median - error,
             high = median + error) %>% 
      ungroup()
    
    # plot
    library(gghalves)
    
    df_alpha %>% 
      
      ggplot(aes(x = noise, y = alpha)) +
      geom_half_violin(aes(color = noise)) +
      geom_half_boxplot(aes(color = noise), side = "r") +
      geom_half_point(shape = 1, size = 1, alpha = 0.5, side = "r") +
      
      ylim(0,2) +
      theme_bw() +
      theme(legend.position = "none")
  
  
  
  
  
  

  
  
  
  
  
