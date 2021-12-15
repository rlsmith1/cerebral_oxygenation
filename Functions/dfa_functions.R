

library(tidyverse)

## source this script to run DFA

# write individual DFA functions, then combine into one big function to run at once
# assumes sampling freq of 50 Hz

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

# alpha <- lm(log10(avg_fluctuation) ~ log10(window_size), data = df_avg_fluct) %>% .$coefficients %>% .[2]



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
