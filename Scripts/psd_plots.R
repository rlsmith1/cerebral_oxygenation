


# libraries ---------------------------------------------------------------


  library(tidyverse)
  library(reticulate)
  library(purrr)
  library(data.table)


  
# data --------------------------------------------------------------------

  
  load("filtered_signals.Rdata")
  load("filtered_muscle_signals.Rdata")
  

# python script to calculate freq & power for each patient ---------------------------------


  np <- import("numpy")
  
  sampling_rate <- 50 # Hz
  
  f_psd <- function(df, num){
    
    df_id <- df %>% count(number, Status) %>% dplyr::select(-n)
    
    df <- df %>% filter(number == num)
    fourier_transform <- np$fft$rfft(df$sg_filt)  
    abs_fourier_transform <- np$abs(fourier_transform)  
    power_spectrum <- np$square(abs_fourier_transform)  
    frequency <- np$linspace(0, sampling_rate/2, length(power_spectrum))  
    df_psd <- tibble(number = num, freq = frequency, power = power_spectrum)
    
    df_psd <- df_psd %>% left_join(df_id, by = "number")
    
    df_psd <- df_psd %>% 
      dplyr::filter(freq > 0.001) %>% 
      mutate(freq_range = cut(freq, seq(from = 0, to = max(freq), by = 0.001)))
    
    df_psd
    
  }

  

# run function on all signals ---------------------------------------------



  # Brain Hb_tot
  df_brain_hbtot_psd <- unique(df_thc_filt_pass$number) %>% 
    purrr::map(~f_psd(df_thc_filt_pass, .x)) %>% 
    rbindlist() %>% 
    as_tibble() %>% 
    group_by(Status, freq_range) %>% 
    summarise(med_power = median(power)) %>% 
    mutate(tissue = "brain", variable = "Hb_tot")
  
  # Brain Hb_oxy
  df_brain_hboxy_psd <- unique(df_o2sat_filt_pass$number) %>% 
    purrr::map(~f_psd(df_o2sat_filt_pass, .x)) %>% 
    rbindlist() %>% 
    as_tibble() %>% 
    group_by(Status, freq_range) %>% 
    summarise(med_power = median(power)) %>% 
    mutate(tissue = "brain", variable = "Hb_oxy")
  
  # Muscle Hb_tot
  df_muscle_hbtot_psd <- unique(df_hbtot_filt_pass_muscle$number) %>% 
    purrr::map(~f_psd(df_hbtot_filt_pass_muscle, .x)) %>% 
    rbindlist() %>% 
    as_tibble() %>% 
    group_by(Status, freq_range) %>% 
    summarise(med_power = median(power)) %>% 
    mutate(tissue = "muscle", variable = "Hb_tot")
  
  # Muscle Hb_oxy
  df_muscle_hboxy_psd <- unique(df_hboxy_filt_pass_muscle$number) %>% 
    purrr::map(~f_psd(df_hboxy_filt_pass_muscle, .x)) %>% 
    rbindlist() %>% 
    as_tibble() %>% 
    group_by(Status, freq_range) %>% 
    summarise(med_power = median(power)) %>% 
    mutate(tissue = "muscle", variable = "Hb_oxy")
  
  
  # combine
  df_psd <- bind_rows(df_brain_hbtot_psd, df_brain_hboxy_psd, df_muscle_hbtot_psd, df_muscle_hboxy_psd)
  df_psd <- df_psd %>% 
    mutate(Status = replace(Status, Status == "HV", "HC")) %>% 
    dplyr::filter(Status != "?")


# plot --------------------------------------------------------------------

  
  # convert freq range to numeric
    freq <- df_psd$freq_range %>% 
      as.character() %>% 
      purrr::map(~strsplit(.x, split = c(","))[[1]][1]) %>% 
      purrr::map(~substr(.x, start = 2, stop = nchar(.x))) %>% 
      as.numeric()
  
    df_psd <- df_psd %>% ungroup() %>% mutate(freq = freq)
    
    # save object for plotting
    save(df_psd, file = "psd_plot_data.Rdata")
    
    
    # plot
    f_plot_psd <- function(df, status, upper_lim) {
      
      df %>% 
        dplyr::filter(Status == status) %>% 
        dplyr::filter(freq > 0.01 & freq < upper_lim) %>% 
        ggplot(aes(x = freq, y = med_power)) +
        geom_line() +
        geom_vline(xintercept = 0.1, color = "red", lty = 2) +
        facet_grid(tissue ~ variable, scales = "free_y") +
        labs(x = "Frequency", y = "Median power across patients") +
        theme_bw()  
      
    }
    

    ## 0-1 Hz
    
        # HC
        f_plot_psd(df_psd, "HC", 1)
    
        # UM
        f_plot_psd(df_psd, "UM", 1)
    
        # CM
        f_plot_psd(df_psd, "CM", 1)
    
    ## 0-3 Hz
        
        # HC
        f_plot_psd(df_psd, "HC", 3)
        
        # UM
        f_plot_psd(df_psd, "UM", 3)
        
        # CM
        f_plot_psd(df_psd, "CM", 3)
        
    

    