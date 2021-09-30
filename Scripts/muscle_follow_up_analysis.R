


### Muscle data preprocessing and DFA analysis from CM follow-up patients




# read in & format data ---------------------------------------------------



### READ IN DATA

    # libraries needed for this step (and future steps)
    library(tidyverse)
    library(purrr)
    
    # get names of patient files, mine are in the directory Data/NIRS files of brain and muscle/ within my R project
    my_files <- list.files(path = "Data/NIRS files of brain and muscle/files_for_muscle_analysis/follow up/", 
                           pattern = "*.txt", full.names = TRUE, recursive = FALSE) # TM0006CMF is an empty file
    
    # remove duplicate patient files 
    # check through patient files, make sure there are no duplicates of patients & session numbers
    
    # read in .txt files - output is a list of tibbles, where each tibble is a patient reading
    l_raw_data <- 1:length(my_files) %>% 
      purrr::map(~read.delim2(file = my_files[.x], header = FALSE, sep = "\t") %>% as_tibble()) 
    
    # name each tibble in list with the patient ID 
    # (as per the file name)
    l_names <- 1:length(my_files) %>% 
      purrr::map(~strsplit(my_files[.x], "//")[[1]][2] %>% 
                   strsplit("[.]") %>% .[[1]] %>% .[1])
    
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
      purrr::map(~dplyr::select(.x, c(ncol(.x), 2, 7:8))) %>% 
      purrr::map(~dplyr::rename(.x, "Time" = "V2", "Hb_oxy" = "V7", "Hb_tot" = "V8")) %>% # muscle data columns
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
      dplyr::select(c(number, subject_id, Status, everything())) 
    
    # make sure patient status is CM
    df_raw_data %>% count(Status)







# signal visualization and segmentation ------------------------------------------------------
    
    
    
### SEGMENTING SIGNALS
    
    # remove signals that aren't at least 200s  
    sampling_rate <- 1/.02 # may need to adjust the sampling rate if it's not 50 Hz
    df_raw_data %>% group_by(number, subject_id) %>% count() %>% dplyr::filter(n < 200*sampling_rate) # will output numbers of signals that aren't at least 200s
    # df_raw_data <- df_raw_data %>% filter(!(number %in% c())) # in c(), include numbers of signals that aren't long enough
    
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
    # For the muscle data, it's important to cut off the signal prior to occlusion
    # check the variance to make sure it isn't astronomically high (I typically used around 10,000 as a cut-off - 
    # a much higher variance may indicate that the signal needs to be cut shorter or is potentially too noisy to include in analysis)
    
    
  # 1: TM0005CMF
    
      # Hb_oxy
      df_raw_data %>% f_plot_hboxy(1)
      df_hboxy_fu1 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(1, start = 1500, end = 2000)
      df_hboxy_fu1 %>% f_plot_hboxy(1) 
      df_hboxy_fu1$Hb_oxy %>% var(na.rm = TRUE)
    
      # Hb_tot
      df_raw_data %>% f_plot_hbtot(1) 
      df_hbtot_fu1 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(1, start = 1500, end = 2000)
      df_hbtot_fu1 %>% f_plot_hbtot(1) 
      df_hbtot_fu1$Hb_tot %>% var(na.rm = TRUE)
      
    
  # 2: TM0007CMF
      
      # Hb_oxy
      df_raw_data %>% f_plot_hboxy(2)
      df_hboxy_fu2 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(2, start = 0, end = 1100)
      df_hboxy_fu2 %>% f_plot_hboxy(2) 
      df_hboxy_fu2$Hb_oxy %>% var(na.rm = TRUE)
      
      # Hb_tot
      df_raw_data %>% f_plot_hbtot(2) 
      df_hbtot_fu2 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(2, start = 0, end = 800)
      df_hbtot_fu2 %>% f_plot_hbtot(2) 
      df_hbtot_fu2$Hb_tot %>% var(na.rm = TRUE)
      
      
  # 3: TM0009CMF
      
      # Hb_oxy
      df_raw_data %>% f_plot_hboxy(3)
      df_hboxy_fu3 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(3, start = 250, end = 500)
      df_hboxy_fu3 %>% f_plot_hboxy(3) 
      df_hboxy_fu3$Hb_oxy %>% var(na.rm = TRUE)
      
      # Hb_tot
      df_raw_data %>% f_plot_hbtot(3) 
      df_hbtot_fu3 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(3, start = 250, end = 600)
      df_hbtot_fu3 %>% f_plot_hbtot(3) 
      df_hbtot_fu3$Hb_tot %>% var(na.rm = TRUE)
      
      
  # 4: TM0018CMF
      
      # Hb_oxy
      df_raw_data %>% f_plot_hboxy(4)
      df_hboxy_fu4 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(4, start = 0, end = 500)
      df_hboxy_fu4 %>% f_plot_hboxy(4) 
      df_hboxy_fu4$Hb_oxy %>% var(na.rm = TRUE)
      
      # Hb_tot
      df_raw_data %>% f_plot_hbtot(4) 
      df_hbtot_fu4 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(4, start = 0, end = 350)
      df_hbtot_fu4 %>% f_plot_hbtot(4) 
      df_hbtot_fu4$Hb_tot %>% var(na.rm = TRUE)
      
      
  # 5: TM0019CMF
      
      # Hb_oxy
      df_raw_data %>% f_plot_hboxy(5)
      df_hboxy_fu5 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(5, start = 150, end = 800)
      df_hboxy_fu5 %>% f_plot_hboxy(5) 
      df_hboxy_fu5$Hb_oxy %>% var(na.rm = TRUE)
      
      # Hb_tot
      df_raw_data %>% f_plot_hbtot(5) 
      df_hbtot_fu5 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(5, start = 150, end = 800)
      df_hbtot_fu5 %>% f_plot_hbtot(5) 
      df_hbtot_fu5$Hb_tot %>% var(na.rm = TRUE)
      
      
  # 6: TM0020CMF
      
      # Hb_oxy
      df_raw_data %>% f_plot_hboxy(6)
      df_hboxy_fu6 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(6, start = 0, end = 800)
      df_hboxy_fu6 %>% f_plot_hboxy(6) 
      df_hboxy_fu6$Hb_oxy %>% var(na.rm = TRUE)
      
      # Hb_tot
      df_raw_data %>% f_plot_hbtot(6) 
      df_hbtot_fu6 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(6, start = 0, end = 800)
      df_hbtot_fu6 %>% f_plot_hbtot(6) 
      df_hbtot_fu6$Hb_tot %>% var(na.rm = TRUE)
      
      
  # 7: TM0021CMF
      
      # Hb_oxy
      df_raw_data %>% f_plot_hboxy(7)
      df_hboxy_fu7 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(7, start = 0, end = 800)
      df_hboxy_fu7 %>% f_plot_hboxy(7) 
      df_hboxy_fu7$Hb_oxy %>% var(na.rm = TRUE)
      
      # Hb_tot
      df_raw_data %>% f_plot_hbtot(7) 
      df_hbtot_fu7 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(7, start = 0, end = 800)
      df_hbtot_fu7 %>% f_plot_hbtot(7) 
      df_hbtot_fu7$Hb_tot %>% var(na.rm = TRUE)
      
      
  # 8: TM0022CMF
      
      # Hb_oxy
      df_raw_data %>% f_plot_hboxy(8)
      df_hboxy_fu8 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(8, start = 0, end = 600)
      df_hboxy_fu8 %>% f_plot_hboxy(8) 
      df_hboxy_fu8$Hb_oxy %>% var(na.rm = TRUE)
      
      # Hb_tot
      df_raw_data %>% f_plot_hbtot(8) 
      df_hbtot_fu8 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(8, start = 0, end = 600)
      df_hbtot_fu8 %>% f_plot_hbtot(8) 
      df_hbtot_fu8$Hb_tot %>% var(na.rm = TRUE)
      
      
  # 9: TM0023CMF
      
      # Hb_oxy
      df_raw_data %>% f_plot_hboxy(9)
      df_hboxy_fu9 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(9, start = 250, end = 750)
      df_hboxy_fu9 %>% f_plot_hboxy(9) 
      df_hboxy_fu9$Hb_oxy %>% var(na.rm = TRUE)
      
      # Hb_tot
      df_raw_data %>% f_plot_hbtot(9) 
      df_hbtot_fu9 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(9, start = 0, end = 800)
      df_hbtot_fu9 %>% f_plot_hbtot(9) 
      df_hbtot_fu9$Hb_tot %>% var(na.rm = TRUE)
      
      
  # 10: TM0024CMF
      
      # Hb_oxy
      df_raw_data %>% f_plot_hboxy(10)
      df_hboxy_fu10 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(10, start = 300, end = 800)
      df_hboxy_fu10 %>% f_plot_hboxy(10) 
      df_hboxy_fu10$Hb_oxy %>% var(na.rm = TRUE)
      
      # Hb_tot
      df_raw_data %>% f_plot_hbtot(10) 
      df_hbtot_fu10 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(10, start = 0, end = 800)
      df_hbtot_fu10 %>% f_plot_hbtot(10) 
      df_hbtot_fu10$Hb_tot %>% var(na.rm = TRUE)
      
      
  # 11: TM0026CMF
      
      # Hb_oxy
      df_raw_data %>% f_plot_hboxy(11)
      df_hboxy_fu11 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(11, start = 100, end = 800)
      df_hboxy_fu11 %>% f_plot_hboxy(11) 
      df_hboxy_fu11$Hb_oxy %>% var(na.rm = TRUE)
      
      # Hb_tot
      df_raw_data %>% f_plot_hbtot(11) 
      df_hbtot_fu11 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(11, start = 0, end = 500)
      df_hbtot_fu11 %>% f_plot_hbtot(11) 
      df_hbtot_fu11$Hb_tot %>% var(na.rm = TRUE)
      
      
  # 12: TM0027CMF
      
      # Hb_oxy
      df_raw_data %>% f_plot_hboxy(12)
      df_hboxy_fu12 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(12, start = 100, end = 800)
      df_hboxy_fu12 %>% f_plot_hboxy(12) 
      df_hboxy_fu12$Hb_oxy %>% var(na.rm = TRUE)
      
      # Hb_tot
      df_raw_data %>% f_plot_hbtot(12) 
      df_hbtot_fu12 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(12, start = 0, end = 800)
      df_hbtot_fu12 %>% f_plot_hbtot(12) 
      df_hbtot_fu12$Hb_tot %>% var(na.rm = TRUE)
      
      
  # 13: TM0028CMF
      
      # Hb_oxy
      df_raw_data %>% f_plot_hboxy(13)
      df_hboxy_fu13 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(13, start = 125, end = 400)
      df_hboxy_fu13 %>% f_plot_hboxy(13) 
      df_hboxy_fu13$Hb_oxy %>% var(na.rm = TRUE)
      
      # Hb_tot
      df_raw_data %>% f_plot_hbtot(13) 
      df_hbtot_fu13 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(13, start = 125, end = 400)
      df_hbtot_fu13 %>% f_plot_hbtot(13) 
      df_hbtot_fu13$Hb_tot %>% var(na.rm = TRUE)
      
      
  # 14: TM0030CMF
      
      # Hb_oxy
      df_raw_data %>% f_plot_hboxy(14)
      df_hboxy_fu14 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(14, start = 0, end = 800)
      df_hboxy_fu14 %>% f_plot_hboxy(14) 
      df_hboxy_fu14$Hb_oxy %>% var(na.rm = TRUE)
      
      # Hb_tot
      df_raw_data %>% f_plot_hbtot(14) 
      df_hbtot_fu14 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(14, start = 0, end = 450)
      df_hbtot_fu14 %>% f_plot_hbtot(14) 
      df_hbtot_fu14$Hb_tot %>% var(na.rm = TRUE)
      
      
  # 15: TM0032CMF
      
      # Hb_oxy
      df_raw_data %>% f_plot_hboxy(15)
      df_hboxy_fu15 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(15, start = 0, end = 800)
      df_hboxy_fu15 %>% f_plot_hboxy(15) 
      df_hboxy_fu15$Hb_oxy %>% var(na.rm = TRUE)
      
      # Hb_tot
      df_raw_data %>% f_plot_hbtot(15) 
      df_hbtot_fu15 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(15, start = 100, end = 500)
      df_hbtot_fu15 %>% f_plot_hbtot(15) 
      df_hbtot_fu15$Hb_tot %>% var(na.rm = TRUE)
      
      
  # 16: TM0060CMF
      
      # Hb_oxy
      df_raw_data %>% f_plot_hboxy(16)
      df_hboxy_fu16 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(16, start = 250, end = 1100)
      df_hboxy_fu16 %>% f_plot_hboxy(16) 
      df_hboxy_fu16$Hb_oxy %>% var(na.rm = TRUE)
      
      # Hb_tot
      df_raw_data %>% f_plot_hbtot(16) 
      df_hbtot_fu16 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(16, start = 775, end = 1100)
      df_hbtot_fu16 %>% f_plot_hbtot(16) 
      df_hbtot_fu16$Hb_tot %>% var(na.rm = TRUE)
      

    # Once you have repeated this for all the patients, combine the filtered Hb_tot signals into one tibble and the filtered Hb_oxy signals into another
    
    df_hboxy_filt_muscle_fu <- bind_rows(df_hboxy_fu1, df_hboxy_fu2, df_hboxy_fu3, df_hboxy_fu4, df_hboxy_fu5,
                                         df_hboxy_fu6, df_hboxy_fu7, df_hboxy_fu8, df_hboxy_fu9, df_hboxy_fu10,
                                         df_hboxy_fu11, df_hboxy_fu12, df_hboxy_fu13, df_hboxy_fu14, df_hboxy_fu15,
                                         df_hboxy_fu16)
    
    df_hbtot_filt_muscle_fu <- bind_rows(df_hbtot_fu1, df_hbtot_fu2, df_hbtot_fu3, df_hbtot_fu4, df_hbtot_fu5,
                                         df_hbtot_fu6, df_hbtot_fu7, df_hbtot_fu8, df_hbtot_fu9, df_hbtot_fu10,
                                         df_hbtot_fu11, df_hbtot_fu12, df_hbtot_fu13, df_hbtot_fu14, df_hbtot_fu15,
                                         df_hbtot_fu16)
    
    
    
    
    
    
    
# Pre-process segmented signals -------------------------------------------
    
    
    
### CALCULATE AVERAGE HB_TOT AND HB_OXY FOR EACH SIGNAL
    
    
    # Hb_tot
    df_hbtot_mean <- df_hbtot_filt_muscle_fu %>% 
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
      ggtitle("muscle follow-up hemoglobin concentration") +
      theme(legend.position = "none") 
    
    
    # Hb_oxy
    df_hboxy_mean <- df_hboxy_filt_muscle_fu %>% 
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
      ggtitle("Muscle follow-up hemoglobin oxygen saturation") +
      theme(legend.position = "none")
    
    # Explore any nonsensical values, especially any negatives in the Hb_oxy. you may want to remove those patients from the analysis due to unreliable readings
    
    
    
### SAVITZY-GOLAY SIGNAL SMOOTHING
    
    
    # For some reason, the Oxiplex removes milliseconds after a few seconds. Add milliseconds back in to Time
    
    # Hb_tot
    df_hbtot_filt_muscle_fu <- df_hbtot_filt_muscle_fu %>% 
      group_by(number, Time) %>% 
      left_join(count(.)) %>% 
      mutate(Time = ifelse(n == 49, Time + 0.02*row_number(), Time + 0.02*(row_number() - 1))) %>% 
      dplyr::select(-n)
    
    # Hb_oxy
    df_hboxy_filt_muscle_fu <- df_hboxy_filt_muscle_fu %>% 
      group_by(number, Time) %>% 
      left_join(count(.)) %>% 
      mutate(Time = ifelse(n == 49, Time + 0.02*row_number(), Time + 0.02*(row_number() - 1))) %>% 
      dplyr::select(-n)
    
    
    # 4th order savitzy-golay filter
    
    library(signal)
    
    # Hb_tot
    df_hbtot_filt_muscle_fu <- df_hbtot_filt_muscle_fu %>% 
      group_by(number) %>% 
      mutate(sg_filt = sgolayfilt(Hb_tot, p = 4)) %>% 
      ungroup()
    
    # plot to visualize smoothing
    df_hbtot_filt_muscle_fu %>% dplyr::filter(number == 10) %>% # can pick any arbitrary time segment to look at
      ggplot(aes(x = Time)) +
      geom_line(aes(y = Hb_tot)) +
      geom_line(aes(y = sg_filt), color = "red")
    
    # Hb_oxy 
    df_hboxy_filt_muscle_fu <- df_hboxy_filt_muscle_fu %>% 
      group_by(number) %>% 
      mutate(sg_filt = sgolayfilt(Hb_oxy, p = 4)) %>% 
      ungroup()
    
    # plot
    df_hboxy_filt_muscle_fu %>% dplyr::filter(number == 10 & Time < 1500 & Time > 1490) %>% 
      ggplot(aes(x = Time)) +
      geom_line(aes(y = Hb_oxy)) +
      geom_line(aes(y = sg_filt), color = "red")
    
    
    
### HIGH-PASS SIGNAL FILTERING
    
    
    # use a high-pass filter to remove anything below 0.01 Hz (associated with head displacements and motion noise) 
    # NOTE: we did NOT use this step in the 2021 paper
    
    
    # load required libraries
    library(seewave)
    library(data.table)
    
    # Hb_tot
    l_hbtot_filt_muscle_fu <- df_hbtot_filt_muscle_fu %>% split(df_hbtot_filt_muscle_fu$number) # split filtered signals into list of tibbles
    sampling_rate <- 50 # again, may need to adjust with different sampling rates
    butter_filter <- butter(4, W = 0.01/(sampling_rate/2), type = "high") # create butter filter
    l_hbtot_filt_pass_muscle_fu <- 1:length(l_hbtot_filt_muscle_fu) %>% 
      purrr::map(~mutate(l_hbtot_filt_muscle_fu[[.x]], high_pass = filtfilt(butter_filter, sg_filt))) # apply to SG filtered signal
    
    df_hbtot_filt_pass_muscle_fu <- rbindlist(l_hbtot_filt_pass_muscle_fu) %>% as_tibble() # combing again into tibble
    
    # repeat for Hb_oxy
    l_hboxy_filt_muscle_fu <- df_hboxy_filt_muscle_fu %>% split(df_hboxy_filt_muscle_fu$number)
    l_hboxy_filt_pass_muscle_fu <- 1:length(l_hboxy_filt_muscle_fu) %>% 
      purrr::map(~mutate(l_hboxy_filt_muscle_fu[[.x]], high_pass = filtfilt(butter_filter, sg_filt)))
    
    df_hboxy_filt_pass_muscle_fu <- rbindlist(l_hboxy_filt_pass_muscle_fu) %>% as_tibble()
    
    
    
### FAST-FOURIER TRANSFORM
    
    
    # perform and visualize the FFT for each signal. we didn't really use this for the 2021 paper either
    
    
    ## function to plot fast fourier transform of signal
    # df = df_hboxy_filt_pass or df_hbtot_filt_pass (NIRS signal + sg filter + high pass)
    # num = patient number to plot
    # zoom = look at overall FFT (FALSE) or freq range of interest (0-1.5 Hz; TRUE)
    
    
    f_plot_fft <- function(df, num, zoom = TRUE) {
      
      # compute fft
      
      fft <- df %>% dplyr::filter(number == num) %>% .$sg_filt %>% fft() 
      
      # determine power spectra
      
      freq <- 50  #sample frequency in Hz 
      duration <- df %>% dplyr::filter(number == num) %>% nrow()/freq # length of signal in seconds
      amo <- Mod(fft) # frequency "amounts" (power)
      freqvec <- 1:length(amo) # associated frequency
      
      freqvec <- freqvec/duration # normalize to signal length to get frequency ranges
      df <- tibble(freq = freqvec, power = amo) # create df assigning power to each frequency
      df <- df[(1:as.integer(0.5*freq*duration)),] # select rows within Nyquist frequency
      
      # plot power spectra
      p <- df %>% dplyr::filter(power < 4000 & freq < ifelse(zoom == TRUE, 1.5, max(freq))) %>% 
        
        ggplot(aes(x = freq, y = power)) + 
        geom_line(stat = "identity") +
        geom_vline(xintercept = 0.01, lty = 2, color = "red") +
        
        theme_bw() +
        
        theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
      
      ifelse(zoom == TRUE, p <- p + scale_x_continuous(breaks = seq(0, 1.5, 0.1)), p <- p)
      
      p
      
    }
    
    
    f_plot_fft(df_hboxy_filt_pass_muscle_fu, 10, zoom = FALSE)  
    f_plot_fft(df_hbtot_filt_pass_muscle_fu, 2, zoom = TRUE)  
    
    
    
### SAVE THE FILTERED SIGNALS
    
    # This is important so you don't have to repeat all the signal pre-processing before running the DFA
    

    # e.g.
    save(df_hbtot_filt_muscle_fu, df_hbtot_filt_pass_muscle_fu, 
         df_hboxy_filt_muscle_fu, df_hboxy_filt_pass_muscle_fu,
         file = "filtered_muscle_fu_signals.Rdata")
    
    
    
    
    
    
    
# Detrended fluctuation analysis ------------------------------------------
    
    
### LOAD IN DATA
    
    # If not already in your environment, load the processed signals
    
    load("filtered_muscle_fu_signals.Rdata")
    
    
    
### COMPUTE CUMULATIVE SUM OF EACH SIGNAL
    
    
    # DFA is performed on the cumulative sum of the signal, not the signal itself
    
  ## Hb_tot
    
    # create list, where each object is a patient
    l_hbtot_filt_pass_muscle_fu <- df_hbtot_filt_pass_muscle_fu %>% split(f = df_hbtot_filt_pass_muscle_fu$number)
    
    # take the mean of each filtered signal
    l_means_muscle_fu <- 1:length(l_hbtot_filt_pass_muscle_fu) %>% 
      purrr::map(~mean(l_hbtot_filt_pass_muscle_fu[[.x]]$sg_filt, na.rm = TRUE))
    
    # subtract the mean and take the cumulative sum
    l_hbtot_cumsum_muscle_fu <- 1:length(l_means_muscle_fu) %>% purrr::map(~mutate(l_hbtot_filt_pass_muscle_fu[[.x]], 
                                                                             mean = l_means_muscle_fu[[.x]]))
    df_hbtot_cumsum_muscle_fu <- rbindlist(l_hbtot_cumsum_muscle_fu) %>% as_tibble()
    
    df_hbtot_cumsum_muscle_fu <- df_hbtot_cumsum_muscle_fu %>% 
      mutate(subtract_mean = sg_filt - mean) %>% 
      group_by(number) %>% 
      mutate(cumsum = cumsum(subtract_mean)) %>% 
      ungroup()
    
    # plot
    df_hbtot_cumsum_muscle_fu %>% dplyr::filter(number == 1) %>% 
      ggplot(aes(x = Time)) +
      geom_line(aes(y = cumsum)) +
      geom_line(aes(y = sg_filt), color = "red") +
      
      theme_bw() +
      ggtitle("Cumulative sum of Savitzy-Golay smoothed signal")
    
    
  ## Hb_oxy
    
    # create list, where each object is a patient
    l_hboxy_filt_pass_muscle_fu <- df_hboxy_filt_pass_muscle_fu %>% split(f = df_hboxy_filt_pass_muscle_fu$number)
    
    # take the mean of each filtered signal
    l_means_muscle_fu <- 1:length(l_hboxy_filt_pass_muscle_fu) %>% 
      purrr::map(~mean(l_hboxy_filt_pass_muscle_fu[[.x]]$sg_filt, na.rm = TRUE))
    
    # subtract the mean and take the cumulative sum
    l_hboxy_cumsum_muscle_fu <- 1:length(l_means_muscle_fu) %>% purrr::map(~mutate(l_hboxy_filt_pass_muscle_fu[[.x]], 
                                                                             mean = l_means_muscle_fu[[.x]]))
    df_hboxy_cumsum_muscle_fu <- rbindlist(l_hboxy_cumsum_muscle_fu) %>% as_tibble()
    
    df_hboxy_cumsum_muscle_fu <- df_hboxy_cumsum_muscle_fu %>% 
      mutate(subtract_mean = sg_filt - mean) %>% 
      group_by(number) %>% 
      mutate(cumsum = cumsum(subtract_mean)) %>% 
      ungroup()
    
    # plot
    df_hboxy_cumsum_muscle_fu %>% dplyr::filter(number == 1) %>% 
      ggplot(aes(x = Time)) +
      geom_line(aes(y = cumsum)) +
      geom_line(aes(y = sg_filt), color = "red") +
      
      theme_bw() +
      ggtitle("Cumulative sum of Savitzy-Golay smoothed signal")
    
    
    
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
    l_hbtot_cumsum_muscle_fu <- df_hbtot_cumsum_muscle_fu %>% split(f = df_hbtot_cumsum_muscle_fu$number)
    l_hboxy_cumsum_muscle_fu <- df_hboxy_cumsum_muscle_fu %>% split(f = df_hboxy_cumsum_muscle_fu$number)
    
    
  # repeat for all patients!!
    
    library(doParallel) # runs faster in parallel
    registerDoParallel()
    
    
    l_dfa_hbtot_muscle_fu <- 1:length(l_hbtot_cumsum_muscle_fu) %>% 
      purrr::map(~f_dfa(l_hbtot_cumsum_muscle_fu[[.x]], num_of_windows = 10))
    
    l_dfa_hboxy_muscle_fu <- 1:length(l_hboxy_cumsum_muscle_fu) %>% 
      purrr::map(~f_dfa(l_hboxy_cumsum_muscle_fu[[.x]], num_of_windows = 10))
    
    
    # would be wise to save DFA results here because it takes a while to run
    save(l_dfa_hbtot_muscle_fu, l_dfa_hboxy_muscle_fu, file = "muscle_fu_dfa_results.Rdata")   
    
    
    
    
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
    hbtot_alphas_muscle_fu <- 1:length(l_dfa_hbtot_muscle_fu) %>% purrr::map(~l_dfa_hbtot_muscle_fu[[.x]][[2]]) %>% unlist()
    
    l_hbtot_second_alpha_muscle_fu <- 1:length(l_dfa_hbtot_muscle_fu) %>% 
      purrr::map(~f_segment_slope(l_dfa_hbtot_muscle_fu[[.x]]))
    
    hbtot_second_alpha_muscle_fu <- 1:length(l_dfa_hbtot_muscle_fu) %>% 
      purrr::map(~f_segment_slope(l_dfa_hbtot_muscle_fu[[.x]])) %>% unlist()
    
    # combine dfa results into df
    df_hbtot_dfa_muscle_fu <- df_hbtot_cumsum_muscle_fu %>% 
      count(number, Status) %>% 
      dplyr::select(-n) %>% 
      mutate(alpha = hbtot_alphas_muscle_fu, 
             second_alpha = hbtot_second_alpha_muscle_fu,
             Status = as.factor(Status))
    
    # Hb_oxy
    hboxy_alphas_muscle_fu <- 1:length(l_dfa_hboxy_muscle_fu) %>% purrr::map(~l_dfa_hboxy_muscle_fu[[.x]][[2]]) %>% unlist()
    
    l_hboxy_second_alpha_muscle_fu <- 1:length(l_dfa_hboxy_muscle_fu) %>% 
      purrr::map(~f_segment_slope(l_dfa_hboxy_muscle_fu[[.x]]))
    
    hboxy_second_alpha_muscle_fu <- 1:length(l_dfa_hboxy_muscle_fu) %>% 
      purrr::map(~f_segment_slope(l_dfa_hboxy_muscle_fu[[.x]])) %>% unlist()
    
    # combine dfa results into df
    df_hboxy_dfa_muscle_fu <- df_hboxy_cumsum_muscle_fu %>% 
      count(number, Status) %>% 
      dplyr::select(-n) %>% 
      mutate(alpha = hboxy_alphas_muscle_fu, 
             second_alpha = hboxy_second_alpha_muscle_fu,
             Status = as.factor(Status))
    
    
    
    
    
# plots and statistical analyses ------------------------------------------
    
    
    library(DescTools)
    
    # plot alphas for Hb_tot and Hb_oxy
    
    
    # Hb_tot
    
    # plot
    df_hbtot_dfa_muscle_fu %>% 
      ggplot(aes(x = Status, y = alpha)) +
      geom_violin(aes(color = Status)) +
      geom_jitter(position = position_jitter(0.2), shape = 1) +
      stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
      
      labs(y = "alpha") +
      ggtitle("Muscle follow-up hemoglobin concentration alpha") +
      
      theme_bw() +
      theme(legend.position = "none")
    

    
    
    # Hb_oxy
    
    # plot
    df_hboxy_dfa_muscle_fu %>% 
      ggplot(aes(x = Status, y = alpha)) +
      geom_violin(aes(color = Status)) +
      geom_jitter(position = position_jitter(0.2), shape = 1) +
      stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
      
      labs(y = "alpha") +
      ggtitle("Muscle follow-up hemoglobin oxygen saturation alpha") +
      
      theme_bw() +
      theme(legend.position = "none")
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
# combine results for export ----------------------------------------------
    
    
    # pair number with subject_id, remove duplicates (ending in 2)
    df_num_id_muscle_fu <- df_hboxy_cumsum_muscle_fu %>% 
      count(number, subject_id) %>% 
      dplyr::select(-n)
    
    # Take means
    df_hbtot_muscle_fu_mean <- df_hbtot_filt_muscle_fu %>% 
      dplyr::filter(!is.na(sg_filt)) %>% 
      group_by(number, Status) %>% 
      summarise(avg_Hb_conc_muscle_fu = mean(sg_filt))
    
    df_hboxy_muscle_fu_mean <- df_hboxy_filt_muscle_fu %>% 
      dplyr::filter(!is.na(sg_filt)) %>% 
      group_by(number, Status) %>% 
      summarise(avg_Hb_o2sat_muscle_fu = mean(sg_filt))
    
    # combine into df
    
    df_muscle_fu_results <- df_num_id_muscle_fu %>% 
      left_join(df_hbtot_dfa_muscle_fu, by = "number") %>% 
      dplyr::rename("muscle_Hb_conc_fu_alpha" = "alpha", "muscle_Hb_conc_fu_second_alpha" = "second_alpha") %>% 
      left_join(df_hboxy_dfa_muscle_fu, by = c("number", "Status")) %>% 
      rename("muscle_Hb_o2sat_fu_alpha" = "alpha", "muscle_Hb_o2sat_fu_second_alpha" = "second_alpha") %>% 
      
      # add hbtot and hboxy means
      left_join(df_hbtot_muscle_fu_mean, by = "number") %>% 
      left_join(df_hboxy_muscle_fu_mean, by = "number") %>% 
      
      # reorder
      dplyr::select(number, subject_id, Status, 
                    avg_Hb_conc_muscle_fu, muscle_Hb_conc_fu_alpha, muscle_Hb_conc_fu_second_alpha,
                    avg_Hb_o2sat_muscle_fu, muscle_Hb_o2sat_fu_alpha, muscle_Hb_o2sat_fu_second_alpha)
    
    # save!!!
    save(df_muscle_fu_results, file = "muscle_follow_up_results.Rdata")
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    