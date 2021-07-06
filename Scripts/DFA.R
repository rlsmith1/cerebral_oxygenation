



# load in processed signal data -------------------------------------------

  load("filtered_signals.Rdata")



# compute cumulative sum --------------------------------------------------

    # create list, where each object is a patient
    l_thc_filt_pass <- df_thc_filt_pass %>% split(f = df_thc_filt_pass$number)

    # take the mean of each filtered signal
    l_means <- 1:length(l_thc_filt_pass) %>% purrr::map(~mean(l_thc_filt_pass[[.x]]$sg_filt, na.rm = TRUE))
    
    # subtract the mean and take the cumulative sum
    l_thc_cumsum <- 1:length(l_means) %>% purrr::map(~mutate(l_thc_filt_pass[[.x]], mean = l_means[[.x]]))
    df_thc_cumsum <- rbindlist(l_thc_cumsum) %>% as_tibble()
    
    df_thc_cumsum <- df_thc_cumsum %>% 
      mutate(subtract_mean = sg_filt - mean) %>% 
      group_by(number) %>% 
      mutate(cumsum = cumsum(subtract_mean)) %>% 
      ungroup()
    

    # plot
    df_thc_cumsum %>% dplyr::filter(number == 1) %>% 
      ggplot(aes(x = Time)) +
      geom_line(aes(y = cumsum)) +
      geom_line(aes(y = sg_filt), color = "red") +
      
      theme_bw() +
      ggtitle("Cumulative sum of Savitzy-Golar smoothed signal")
        

        
# split into windows --------------------------------------------------------------------
    
    
    
    # start with signal from patient 1
    df_thc_cumsum1 %>% nrow()/50 # signal is 1169 seconds long
    
    # determine window lengths
    
      # criteria (Hardstone et al 2012):
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
        geom_line(aes(y = sg_filt), color = "red", alpha = 0.5) +
        geom_line(aes(y = cumsum), color = "black", alpha = 0.65) +
        geom_smooth(aes(y = cumsum), method = "lm", color = "blue", se = FALSE) +
        theme_bw()
      
    }
    
    f_plot_window(l_splits[[11]], 1)

    # plot detrended window
    p_plot_resids <- function(l_resids, window_num) {
      
      l_resids[[window_num]] %>% 
        ggplot(aes(x = Time, y = residuals)) +
        geom_line(color = "black") +
        geom_hline(yintercept = 0, color = "blue") +
        theme_bw()
      
    }
    
    p_plot_resids(l_resids[[11]], 1)
    

    

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
    
    
    # create list of df cumsum, where each object is a patient
    l_thc_cumsum <- df_thc_cumsum %>% split(f = df_thc_cumsum$number)
    
    # repeat for all patients!!
    
    library(doParallel)
    registerDoParallel()
    
    
    l_dfa <- 1:length(l_thc_cumsum) %>% purrr::map(~f_dfa(l_thc_cumsum[[.x]], num_of_windows = 10))

    
    # save DFA analysis results
    save(l_dfa, l_dfa2, file = "dfa_results.Rdata")   
    
    load("dfa_results.Rdata")
    

# second order detrending --------------------------------------------------

## separately compute scaling exponents for short-length vs long-length windws, for short-range and long-range correlations of the time series (Seleznov et al 2019)
    
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
    
    alphas <- 1:length(l_dfa) %>% purrr::map(~l_dfa[[.x]][[2]]) %>% unlist()
    
    l_second_alpha <- 1:length(l_dfa) %>% purrr::map(~f_segment_slope(l_dfa[[.x]]))
    second_alpha <- 1:length(l_dfa) %>% purrr::map(~f_segment_slope(l_dfa[[.x]])) %>% unlist()
    
    # combine dfa results into df
    df_dfa_res <- df_thc_cumsum %>% 
      count(number, Status) %>% 
      dplyr::select(-n) %>% 
      mutate(alpha = alphas, 
             second_alpha = second_alpha,
             Status = as.factor(Status))
    


          

# plot & explore results --------------------------------------

    # plot example of detrending cumulative sum
    df_thc_cumsum1 <- df_thc_cumsum %>% dplyr::filter(number == 1)
    
    window_size <- window_lengths[11]*50
    signal_length <- nrow(df_thc_cumsum1)
    reps <- rep(1:ceiling(signal_length/window_size), each = window_size)[1:signal_length]
    df_thc_cumsum1$window_num <- reps
    
    ggplot(data = df_thc_cumsum1, aes(x = Time)) +
      geom_line(aes(y = sg_filt), color = "red", alpha = 0.5) +
      geom_line(aes(y = cumsum), alpha = 0.5) +
      
      # draw in windows
      geom_vline(xintercept = df_thc_cumsum1[1,]$Time + 0*window_lengths[11], lty = 2) +
      geom_vline(xintercept = df_thc_cumsum1[1,]$Time + 1*window_lengths[11], lty = 2) +
      geom_vline(xintercept = df_thc_cumsum1[1,]$Time + 2*window_lengths[11], lty = 2) +
      geom_vline(xintercept = df_thc_cumsum1[1,]$Time + 3*window_lengths[11], lty = 2) +
      geom_vline(xintercept = df_thc_cumsum1[1,]$Time + 4*window_lengths[11], lty = 2) +
      geom_vline(xintercept = df_thc_cumsum1[1,]$Time + 5*window_lengths[11], lty = 2) +
      geom_vline(xintercept = df_thc_cumsum1[1,]$Time + 6*window_lengths[11], lty = 2) +
      geom_vline(xintercept = df_thc_cumsum1[1,]$Time + 7*window_lengths[11], lty = 2) +
      geom_vline(xintercept = df_thc_cumsum1[1,]$Time + 8*window_lengths[11], lty = 2) +
      geom_vline(xintercept = df_thc_cumsum1[1,]$Time + 9*window_lengths[11], lty = 2) +
      geom_vline(xintercept = df_thc_cumsum1[1,]$Time + 10*window_lengths[11], lty = 2) +
      
      geom_smooth(data = df_thc_cumsum1 %>% dplyr::filter(window_num == 1), aes(y = cumsum), method = "lm", se = FALSE) +
      geom_smooth(data = df_thc_cumsum1 %>% dplyr::filter(window_num == 2), aes(y = cumsum), method = "lm", se = FALSE) +
      geom_smooth(data = df_thc_cumsum1 %>% dplyr::filter(window_num == 3), aes(y = cumsum), method = "lm", se = FALSE) +
      geom_smooth(data = df_thc_cumsum1 %>% dplyr::filter(window_num == 4), aes(y = cumsum), method = "lm", se = FALSE) +
      geom_smooth(data = df_thc_cumsum1 %>% dplyr::filter(window_num == 5), aes(y = cumsum), method = "lm", se = FALSE) +
      geom_smooth(data = df_thc_cumsum1 %>% dplyr::filter(window_num == 6), aes(y = cumsum), method = "lm", se = FALSE) +
      geom_smooth(data = df_thc_cumsum1 %>% dplyr::filter(window_num == 7), aes(y = cumsum), method = "lm", se = FALSE) +
      geom_smooth(data = df_thc_cumsum1 %>% dplyr::filter(window_num == 8), aes(y = cumsum), method = "lm", se = FALSE) +
      geom_smooth(data = df_thc_cumsum1 %>% dplyr::filter(window_num == 9), aes(y = cumsum), method = "lm", se = FALSE) +
      geom_smooth(data = df_thc_cumsum1 %>% dplyr::filter(window_num == 10), aes(y = cumsum), method = "lm", se = FALSE) +
      
      labs(y = "THC") +
      theme_bw() +
      ggtitle("Cumulative sum of Savitzy-Golay smoothed signal")
    
    
    
library(DescTools)
    
    # plot DFA results
    
    names(l_dfa) <- names(l_thc_cumsum)  
    
      # segmented
      df1 <- l_dfa[[1]][[1]] %>% mutate(log10_avg_fluct = log10(avg_fluctuation), log10_window_size = log10(window_size))
      fit1 <- lm(log10_avg_fluct ~ log10_window_size, data = df1) 
      segfit1 <- segmented(fit1)
      slopes1 <- coef(segfit1) # the coefficients are the differences in slope in comparison to the previous slope

      
      # first line: 
      #y = b0 + b1*x; y = intercept1 + slope1 * x
      
      # second line:
      # y = c0 + c1*x; y = intercept2 + slope2 * x

      # At the breakpoint (break1), the segments b and c intersect: b0 + b1*x = c0 + c1*x
      
      b0 <- slopes1[[1]]
      b1 <- slopes1[[2]]
      
      c1 <- slopes1[[2]] + slopes1[[3]]
      break1 <- segfit1$psi[[2]]
      
      # Solve for c0 (intercept of second segment):
      c0 <- b0 + b1 * break1 - c1 * break1
    
      l_dfa[[1]][[1]] %>% 
        ggplot(aes(x = log10(window_size), y = log10(avg_fluctuation))) +
        geom_point(shape = 1, size = 3) +
        geom_vline(xintercept = break1, lty = 2) + # add line showing breakpoint
        geom_abline(intercept = b0, slope = b1, color = "#F8766D") +
        geom_abline(intercept = c0, slope = c1, color = "#00BFC4") +
        
        geom_text(x = -0.5, y = 2, label = paste0("alpha = ", round(b1, 2)), color = "#F8766D") +
        geom_text(x = 1.5, y = 1.25, label = paste0("alpha = ", round(c1, 2)), color = "#00BFC4") +
        
        theme_bw() +
        ggtitle("Patient 1 (CM) DFA results")
      
      f_segment_slope(l_dfa[[1]])
      
      # not segmented
      df_combine <-  l_dfa[[6]][[1]] %>% mutate(number = 8, Status = "CM") %>% 
        bind_rows(l_dfa[[60]][[1]] %>% mutate(number = 65, Status = "UM"), 
                  l_dfa[[74]][[1]] %>% mutate(number = 80, Status = "HV"))
      
      df80 <- l_dfa[[74]][[1]] %>% mutate(log10_avg_fluct = log10(avg_fluctuation), log10_window_size = log10(window_size))
      fit80 <- lm(log10_avg_fluct ~ log10_window_size, df80)
      
      df_combine %>% dplyr::filter(number == 80) %>% 
        ggplot(aes(x = log10(window_size), y = log10(avg_fluctuation))) +
        geom_point(shape = 1, size = 3) +
        geom_abline(intercept = coef(fit80)[1], slope = coef(fit80)[2], color = "#F8766D") +
        geom_text(x = -0.5, y = 2, label = paste0("alpha = ", round(coef(fit80)[2], 2))) +
        
        theme_bw() +
        ggtitle("Patient 80 (HV) DFA results")

      
    # plot alpha by group
    df_dfa_res %>% 
      ggplot(aes(x = Status, y = second_alpha)) +
      geom_violin(aes(color = Status)) +
      geom_jitter(position = position_jitter(0.2), shape = 1) +
      stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
      
      labs(y = "alpha") +
      ggtitle("Cerebral tissue hemoglobin concentration alpha") +
    
      theme_bw() +
      theme(legend.position = "none")
    
    # test for differences
    kruskal.test(second_alpha ~ Status, data = df_dfa_res)
    DunnTest(second_alpha ~ Status, data = df_dfa_res)
    
    
    # looks like there could be two distinct distributions in the CM group...
    df_dfa_res2 <- df_dfa_res %>% mutate(Status_split = as.factor(case_when(
      
      Status == "UM" ~ "UM",
      Status == "HV" ~ "HV",
      Status == "CM" & second_alpha > 0.75 ~ "CM_norm",
      Status == "CM" & second_alpha < 0.75 ~ "CM_disrupt"
      
    )) )
    
    # plot split
    df_dfa_res2 %>% 
      ggplot(aes(x = Status, y = second_alpha)) +
      geom_violin(aes(color = Status_split), position = position_dodge(0.0)) +
      geom_jitter(position = position_jitter(0.1), shape = 1) +
      stat_summary(fun = "median", geom = "crossbar", aes(color = Status_split), size = 0.2, width = 0.5) +
      
      labs(y = "alpha") +
      ggtitle("Cerebral tissue hemoglobin concentration alpha") +
      
      theme_bw() +
      theme(legend.position = "none")
    
    # test for differences
    kruskal.test(second_alpha ~ Status_split, data = df_dfa_res2)
    DunnTest(second_alpha ~ Status_split, data = df_dfa_res2)

    
    


# combine results with average THC for export -----------------------------

    df_results <- df_thc_cumsum %>% 
      group_by(number) %>% 
      summarise(mean(THC)) %>% 
      left_join(df_dfa_res, by = "number") %>% 
      rename("Hb_total" = `mean(THC)`) %>% 
      dplyr::select(c(number, Status, everything()))

    save(df_results, file = "results.Rdata")   
    
    df_results_binomial <- df_results %>% 
      dplyr::filter(Status != "HV") %>% 
      mutate(Status = as.factor(ifelse(Status == "CM", 1, 0))) 
    
    glm(Status ~ Hb_total + second_alpha, family = "binomial", data = df_results_binomial) %>% summary()
    
    
    
    
        
    
    
    
    
    
    
    
    



    