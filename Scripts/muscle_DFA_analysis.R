


### Muscle data preprocessing and DFA analysis


# read in & format data ---------------------------------------------------



### READ IN DATA
    
    # libraries needed for this step (and future steps)
    library(tidyverse)
    library(purrr)
    
    # get names of patient files, mine are in the directory Data/NIRS files of brain and muscle/ within my R project
    my_files <- list.files(path = "Data/NIRS files of brain and muscle/files_for_muscle_analysis/initial hospital visit/", 
                             pattern = "*.txt|*.text", full.names = TRUE, recursive = FALSE)

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
      dplyr::select(c(number, subject_id, Status, everything())) %>% 
      mutate(Status = replace(Status, Status == "HV", "HC"),
             Status = replace(Status, Status == "-0", "?"),
             Status = factor(Status, levels = c("HC", "UM", "CM", "?")))
    
    # make sure patient status is one of HC, UM, or CM. adjust those that aren't
    df_raw_data %>% count(Status) # need to figure out what these files are from


    
    
    
# signal visualization and segmentation ------------------------------------------------------
    
    
    
### SEGMENTING SIGNALS
    
    # remove signals that aren't at least 200s  
    sampling_rate <- 1/.02 # may need to adjust the sampling rate if it's not 50 Hz
    df_raw_data %>% group_by(number, subject_id) %>% count() %>% dplyr::filter(n < 200*sampling_rate) # will output numbers of signals that aren't at least 200s
    df_raw_data <- df_raw_data %>% filter(!(number %in% c(19, 24, 25, 82))) # in c(), include numbers of signals that aren't long enough
    
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
    
    
# 1: TM0001CM01
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(1) # eyeball where blips are
    df_hbtot1 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(1, start = 0, end = 1400) # shorten signal to desired length, must be > 200s
    df_hbtot1 %>% f_plot_hbtot(1) # reassess visually
    df_hbtot1$Hb_tot %>% var(na.rm = TRUE) # check var less than 10000
    
    # Hb_oxy - repeat for Hb_oxy signal
    df_raw_data %>% f_plot_hboxy(1)
    df_hboxy1 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(1, start = 0, end = 1400)
    df_hboxy1 %>% f_plot_hboxy(1) 
    df_hboxy1$Hb_oxy %>% var(na.rm = TRUE)
   
     
# 2: TM0002CM01
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(2) 
    df_hbtot2 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(2, start = 0, end = 1200)
    df_hbtot2 %>% f_plot_hbtot(2) 
    df_hbtot2$Hb_tot %>% var(na.rm = TRUE) 
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(2)
    df_hboxy2 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(2, start = 0, end = 1200)
    df_hboxy2 %>% f_plot_hboxy(2) 
    df_hboxy2$Hb_oxy %>% var(na.rm = TRUE)
    

# 3: TM0003CM01
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(3) 
    df_hbtot3 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(3, start = 0, end = 1150)
    df_hbtot3 %>% f_plot_hbtot(3) 
    df_hbtot3$Hb_tot %>% var(na.rm = TRUE) 
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(3)
    df_hboxy3 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(3, start = 0, end = 1150)
    df_hboxy3 %>% f_plot_hboxy(3) 
    df_hboxy3$Hb_oxy %>% var(na.rm = TRUE)
    
  
# 4: TM0005CM01
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(4) 
    df_hbtot4 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(4, start = 1500, end = 2000)
    df_hbtot4 %>% f_plot_hbtot(4) 
    df_hbtot4$Hb_tot %>% var(na.rm = TRUE) 
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(4)
    df_hboxy4 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(4, start = 1500, end = 2000)
    df_hboxy4 %>% f_plot_hboxy(4) 
    df_hboxy4$Hb_oxy %>% var(na.rm = TRUE)
    
    
# 5: TM0006CM01
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(5) 
    df_hbtot5 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(5, start = 0, end = 1150)
    df_hbtot5 %>% f_plot_hbtot(5) 
    df_hbtot5$Hb_tot %>% var(na.rm = TRUE) 
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(5)
    df_hboxy5 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(5, start = 0, end = 1150)
    df_hboxy5 %>% f_plot_hboxy(5) 
    df_hboxy5$Hb_oxy %>% var(na.rm = TRUE)
    
    
# 6: TM0007CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(6)
    df_hboxy6 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(6, start = 0, end = 1150)
    df_hboxy6 %>% f_plot_hboxy(6) 
    df_hboxy6$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(6) 
    df_hbtot6 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(6, start = 0, end = 1150)
    df_hbtot6 %>% f_plot_hbtot(6) 
    df_hbtot6$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 7: TM0010CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(7)
    df_hboxy7 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(7, start = 0, end = 1150)
    df_hboxy7 %>% f_plot_hboxy(7) 
    df_hboxy7$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(7) 
    df_hbtot7 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(7, start = 0, end = 1150)
    df_hbtot7 %>% f_plot_hbtot(7) 
    df_hbtot7$Hb_tot %>% var(na.rm = TRUE) 

    
# 8: TM0012CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(8)
    df_hboxy8 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(8, start = 475, end = 1200)
    df_hboxy8 %>% f_plot_hboxy(8) 
    df_hboxy8$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(8) 
    df_hbtot8 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(8, start = 475, end = 1200)
    df_hbtot8 %>% f_plot_hbtot(8) 
    df_hbtot8$Hb_tot %>% var(na.rm = TRUE) 

    
# 9: TM0013CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(9)
    df_hboxy9 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(9, start = 0, end = 1150)
    df_hboxy9 %>% f_plot_hboxy(9) 
    df_hboxy9$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(9) 
    df_hbtot9 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(9, start = 0, end = 1150)
    df_hbtot9 %>% f_plot_hbtot(9) 
    df_hbtot9$Hb_tot %>% var(na.rm = TRUE) 

    
# 10: TM0014CM01 2
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(10)
    df_hboxy10 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(10, start = 0, end = 1150)
    df_hboxy10 %>% f_plot_hboxy(10) 
    df_hboxy10$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(10) 
    df_hbtot10 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(10, start = 0, end = 1150)
    df_hbtot10 %>% f_plot_hbtot(10) 
    df_hbtot10$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 11: TM0014CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(11)
    df_hboxy11 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(11, start = 0, end = 1150)
    df_hboxy11 %>% f_plot_hboxy(11) 
    df_hboxy11$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(11) 
    df_hbtot11 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(11, start = 0, end = 1150)
    df_hbtot11 %>% f_plot_hbtot(11) 
    df_hbtot11$Hb_tot %>% var(na.rm = TRUE) 

    
# 12: TM0015CM01 2
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(12)
    df_hboxy12 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(12, start = 650, end = 2000)
    df_hboxy12 %>% f_plot_hboxy(12) 
    df_hboxy12$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(12) 
    df_hbtot12 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(12, start = 650, end = 2000)
    df_hbtot12 %>% f_plot_hbtot(12) 
    df_hbtot12$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 13: TM0015CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(13)
    df_hboxy13 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(13, start = 650, end = 2000)
    df_hboxy13 %>% f_plot_hboxy(13) 
    df_hboxy13$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(13) 
    df_hbtot13 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(13, start = 650, end = 2000)
    df_hbtot13 %>% f_plot_hbtot(13) 
    df_hbtot13$Hb_tot %>% var(na.rm = TRUE) 
  
    
# 14: TM0016CM01 2
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(14)
    df_hboxy14 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(14, start = 1600, end = 2000)
    df_hboxy14 %>% f_plot_hboxy(14) 
    df_hboxy14$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(14) 
    df_hbtot14 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(14, start = 1600, end = 2000)
    df_hbtot14 %>% f_plot_hbtot(14) 
    df_hbtot14$Hb_tot %>% var(na.rm = TRUE) 

    
# 15: TM0016CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(15)
    df_hboxy15 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(15, start = 1600, end = 2000)
    df_hboxy15 %>% f_plot_hboxy(15) 
    df_hboxy15$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(15) 
    df_hbtot15 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(15, start = 1600, end = 2000)
    df_hbtot15 %>% f_plot_hbtot(15) 
    df_hbtot15$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 16: TM0018CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(16)
    df_hboxy16 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(16, start = 0, end = 600)
    df_hboxy16 %>% f_plot_hboxy(16) 
    df_hboxy16$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(16) 
    df_hbtot16 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(16, start = 0, end = 1150)
    df_hbtot16 %>% f_plot_hbtot(16) 
    df_hbtot16$Hb_tot %>% var(na.rm = TRUE) 

    
# 17: TM0019CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(17)
    df_hboxy17 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(17, start = 0, end = 750)
    df_hboxy17 %>% f_plot_hboxy(17) 
    df_hboxy17$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(17) 
    df_hbtot17 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(17, start = 300, end = 650)
    df_hbtot17 %>% f_plot_hbtot(17) 
    df_hbtot17$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 18: TM0020CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(18)
    df_hboxy18 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(18, start = 0, end = 750)
    df_hboxy18 %>% f_plot_hboxy(18) 
    df_hboxy18$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(18) 
    df_hbtot18 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(18, start = 0, end = 750)
    df_hbtot18 %>% f_plot_hbtot(18) 
    df_hbtot18$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 20: TM0022CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(20)
    df_hboxy20 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(20, start = 0, end = 800)
    df_hboxy20 %>% f_plot_hboxy(20) 
    df_hboxy20$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(20) 
    df_hbtot20 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(20, start = 0, end = 800)
    df_hbtot20 %>% f_plot_hbtot(20) 
    df_hbtot20$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 21: TM0023CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(21)
    df_hboxy21 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(21, start = 0, end = 800)
    df_hboxy21 %>% f_plot_hboxy(21) 
    df_hboxy21$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(21) 
    df_hbtot21 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(21, start = 0, end = 800)
    df_hbtot21 %>% f_plot_hbtot(21) 
    df_hbtot21$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 22: TM0024CM01 part 2
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(22)
    df_hboxy22 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(22, start = 400, end = 2000)
    df_hboxy22 %>% f_plot_hboxy(22) 
    df_hboxy22$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(22) 
    df_hbtot22 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(22, start = 400, end = 2000)
    df_hbtot22 %>% f_plot_hbtot(22) 
    df_hbtot22$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 23: TM0024CM01
    
    ## too noisy, only use part 2 for this patient
    
    
# 26: TM0027CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(26)
    df_hboxy26 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(26, start = 0, end = 800)
    df_hboxy26 %>% f_plot_hboxy(26) 
    df_hboxy26$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(26) 
    df_hbtot26 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(26, start = 0, end = 800)
    df_hbtot26 %>% f_plot_hbtot(26) 
    df_hbtot26$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 27: TM0028CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(27)
    df_hboxy27 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(27, start = 300, end = 800)
    df_hboxy27 %>% f_plot_hboxy(27) 
    df_hboxy27$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(27) 
    df_hbtot27 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(27, start = 300, end = 800)
    df_hbtot27 %>% f_plot_hbtot(27) 
    df_hbtot27$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 28: TM0029CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(28)
    df_hboxy28 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(28, start = 300, end = 800)
    df_hboxy28 %>% f_plot_hboxy(28) 
    df_hboxy28$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(28) 
    df_hbtot28 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(28, start = 300, end = 800)
    df_hbtot28 %>% f_plot_hbtot(28) 
    df_hbtot28$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 29: TM0030CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(29)
    df_hboxy29 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(29, start = 0, end = 450)
    df_hboxy29 %>% f_plot_hboxy(29) 
    df_hboxy29$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(29) 
    df_hbtot29 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(29, start = 0, end = 450)
    df_hbtot29 %>% f_plot_hbtot(29) 
    df_hbtot29$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 30: TM0031CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(30)
    df_hboxy30 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(30, start = 0, end = 500)
    df_hboxy30 %>% f_plot_hboxy(30) 
    df_hboxy30$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(30) 
    df_hbtot30 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(30, start = 0, end = 500)
    df_hbtot30 %>% f_plot_hbtot(30) 
    df_hbtot30$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 31: TM0032CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(31)
    df_hboxy31 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(31, start = 1250, end = 2000)
    df_hboxy31 %>% f_plot_hboxy(31) 
    df_hboxy31$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(31) 
    df_hbtot31 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(31, start = 1250, end = 2000)
    df_hbtot31 %>% f_plot_hbtot(31) 
    df_hbtot31$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 32: TM0033CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(32)
    df_hboxy32 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(32, start = 0, end = 700)
    df_hboxy32 %>% f_plot_hboxy(32) 
    df_hboxy32$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(32) 
    df_hbtot32 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(32, start = 0, end = 700)
    df_hbtot32 %>% f_plot_hbtot(32) 
    df_hbtot32$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 33: TM0034CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(33)
    df_hboxy33 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(33, start = 0, end = 800)
    df_hboxy33 %>% f_plot_hboxy(33) 
    df_hboxy33$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(33) 
    df_hbtot33 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(33, start = 50, end = 800)
    df_hbtot33 %>% f_plot_hbtot(33) 
    df_hbtot33$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 34: TM0035CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(34)
    df_hboxy34 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(34, start = 0, end = 800)
    df_hboxy34 %>% f_plot_hboxy(34) 
    df_hboxy34$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(34) 
    df_hbtot34 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(34, start = 100, end = 700)
    df_hbtot34 %>% f_plot_hbtot(34) 
    df_hbtot34$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 35: TM0036CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(35)
    df_hboxy35 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(35, start = 0, end = 800)
    df_hboxy35 %>% f_plot_hboxy(35) 
    df_hboxy35$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(35) 
    df_hbtot35 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(35, start = 0, end = 800)
    df_hbtot35 %>% f_plot_hbtot(35) 
    df_hbtot35$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 36: TM0037CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(36)
    df_hboxy36 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(36, start = 0, end = 800)
    df_hboxy36 %>% f_plot_hboxy(36) 
    df_hboxy36$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(36) 
    df_hbtot36 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(36, start = 400, end = 750)
    df_hbtot36 %>% f_plot_hbtot(36) 
    df_hbtot36$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 37: TM0038CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(37)
    df_hboxy37 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(37, start = 1200, end = 2000)
    df_hboxy37 %>% f_plot_hboxy(37) 
    df_hboxy37$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(37) 
    df_hbtot37 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(37, start = 1200, end = 2000)
    df_hbtot37 %>% f_plot_hbtot(37) 
    df_hbtot37$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 38: TM0039CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(38)
    df_hboxy38 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(38, start = 0, end = 800)
    df_hboxy38 %>% f_plot_hboxy(38) 
    df_hboxy38$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(38) 
    df_hbtot38 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(38, start = 0, end = 800)
    df_hbtot38 %>% f_plot_hbtot(38) 
    df_hbtot38$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 39: TM0040CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(39)
    df_hboxy39 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(39, start = 0, end = 800)
    df_hboxy39 %>% f_plot_hboxy(39) 
    df_hboxy39$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(39) 
    df_hbtot39 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(39, start = 0, end = 800)
    df_hbtot39 %>% f_plot_hbtot(39) 
    df_hbtot39$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 40: TM0041CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(40)
    df_hboxy40 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(40, start = 0, end = 800)
    df_hboxy40 %>% f_plot_hboxy(40) 
    df_hboxy40$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(40) 
    df_hbtot40 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(40, start = 0, end = 800)
    df_hbtot40 %>% f_plot_hbtot(40) 
    df_hbtot40$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 41: TM0042CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(41)
    df_hboxy41 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(41, start = 0, end = 800)
    df_hboxy41 %>% f_plot_hboxy(41) 
    df_hboxy41$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(41) 
    df_hbtot41 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(41, start = 0, end = 800)
    df_hbtot41 %>% f_plot_hbtot(41) 
    df_hbtot41$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 42: TM0043CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(42)
    df_hboxy42 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(42, start = 0, end = 800)
    df_hboxy42 %>% f_plot_hboxy(42) 
    df_hboxy42$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(42) 
    df_hbtot42 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(42, start = 0, end = 800)
    df_hbtot42 %>% f_plot_hbtot(42) 
    df_hbtot42$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 43: TM0044CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(43)
    df_hboxy43 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(43, start = 0, end = 800)
    df_hboxy43 %>% f_plot_hboxy(43) 
    df_hboxy43$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(43) 
    df_hbtot43 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(43, start = 0, end = 800)
    df_hbtot43 %>% f_plot_hbtot(43) 
    df_hbtot43$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 44: TM0045CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(44)
    df_hboxy44 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(44, start = 1250, end = 2000)
    df_hboxy44 %>% f_plot_hboxy(44) 
    df_hboxy44$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(44) 
    df_hbtot44 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(44, start = 300, end = 600)
    df_hbtot44 %>% f_plot_hbtot(44) 
    df_hbtot44$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 45: TM0046CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(45)
    df_hboxy45 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(45, start = 1250, end = 2000)
    df_hboxy45 %>% f_plot_hboxy(45) 
    df_hboxy45$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(45) 
    df_hbtot45 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(45, start = 1250, end = 2000)
    df_hbtot45 %>% f_plot_hbtot(45) 
    df_hbtot45$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 46: TM0047CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(46)
    df_hboxy46 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(46, start = 0, end = 600)
    df_hboxy46 %>% f_plot_hboxy(46) 
    df_hboxy46$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(46) 
    df_hbtot46 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(46, start = 0, end = 600)
    df_hbtot46 %>% f_plot_hbtot(46) 
    df_hbtot46$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 47: TM0048CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(47)
    df_hboxy47 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(47, start = 0, end = 800)
    df_hboxy47 %>% f_plot_hboxy(47) 
    df_hboxy47$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(47) 
    df_hbtot47 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(47, start = 0, end = 800)
    df_hbtot47 %>% f_plot_hbtot(47) 
    df_hbtot47$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 48: TM0060CM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(48)
    df_hboxy48 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(48, start = 0, end = 1150)
    df_hboxy48 %>% f_plot_hboxy(48) 
    df_hboxy48$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(48) 
    df_hbtot48 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(48, start = 550, end = 1150)
    df_hbtot48 %>% f_plot_hbtot(48) 
    df_hbtot48$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 49: TM1003UM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(49)
    df_hboxy49 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(49, start = 0, end = 1150)
    df_hboxy49 %>% f_plot_hboxy(49) 
    df_hboxy49$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(49) 
    df_hbtot49 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(49, start = 0, end = 1000)
    df_hbtot49 %>% f_plot_hbtot(49) 
    df_hbtot49$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 50: TM1004UM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(50)
    df_hboxy50 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(50, start = 50, end = 800)
    df_hboxy50 %>% f_plot_hboxy(50) 
    df_hboxy50$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(50) 
    df_hbtot50 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(50, start = 500, end = 1000)
    df_hbtot50 %>% f_plot_hbtot(50) 
    df_hbtot50$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 51: TM1005UM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(51)
    df_hboxy51 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(51, start = 0, end = 1200)
    df_hboxy51 %>% f_plot_hboxy(51) 
    df_hboxy51$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(51) 
    df_hbtot51 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(51, start = 0, end = 800)
    df_hbtot51 %>% f_plot_hbtot(51) 
    df_hbtot51$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 52: TM1006UM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(52)
    df_hboxy52 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(52, start = 250, end = 1200)
    df_hboxy52 %>% f_plot_hboxy(52) 
    df_hboxy52$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(52) 
    df_hbtot52 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(52, start = 300, end = 800)
    df_hbtot52 %>% f_plot_hbtot(52) 
    df_hbtot52$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 53: TM1007UM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(53)
    df_hboxy53 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(53, start = 100, end = 1200)
    df_hboxy53 %>% f_plot_hboxy(53) 
    df_hboxy53$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(53) 
    df_hbtot53 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(53, start = 100, end = 1200)
    df_hbtot53 %>% f_plot_hbtot(53) 
    df_hbtot53$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 54: TM1008UM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(54)
    df_hboxy54 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(54, start = 0, end = 1100)
    df_hboxy54 %>% f_plot_hboxy(54) 
    df_hboxy54$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(54) 
    df_hbtot54 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(54, start = 0, end = 1100)
    df_hbtot54 %>% f_plot_hbtot(54) 
    df_hbtot54$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 55: TM1009UM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(55)
    df_hboxy55 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(55, start = 100, end = 1000)
    df_hboxy55 %>% f_plot_hboxy(55) 
    df_hboxy55$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(55) 
    df_hbtot55 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(55, start = 0, end = 1100)
    df_hbtot55 %>% f_plot_hbtot(55) 
    df_hbtot55$Hb_tot %>% var(na.rm = TRUE) 
    
# 56: TM1010UM01
    
    ## signal too noisy, no segment long enough to analyze
    

# 57: TM1011UM01
    
    ## signal too noisy
    
    
# 58: TM1012UM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(58)
    df_hboxy58 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(58, start = 0, end = 900)
    df_hboxy58 %>% f_plot_hboxy(58) 
    df_hboxy58$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(58) 
    df_hbtot58 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(58, start = 0, end = 900)
    df_hbtot58 %>% f_plot_hbtot(58) 
    df_hbtot58$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 59: TM1013UM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(59)
    df_hboxy59 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(59, start = 150, end = 800)
    df_hboxy59 %>% f_plot_hboxy(59) 
    df_hboxy59$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(59) 
    df_hbtot59 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(59, start = 150, end = 800)
    df_hbtot59 %>% f_plot_hbtot(59) 
    df_hbtot59$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 60: TM1014UM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(60)
    df_hboxy60 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(60, start = 0, end = 1100)
    df_hboxy60 %>% f_plot_hboxy(60) 
    df_hboxy60$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(60) 
    df_hbtot60 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(60, start = 0, end = 1100)
    df_hbtot60 %>% f_plot_hbtot(60) 
    df_hbtot60$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 61: TM1015-01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(61)
    df_hboxy61 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(61, start = 550, end = 800)
    df_hboxy61 %>% f_plot_hboxy(61) 
    df_hboxy61$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(61) 
    df_hbtot61 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(61, start = 400, end = 800)
    df_hbtot61 %>% f_plot_hbtot(61) 
    df_hbtot61$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 62: TM1015UM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(62)
    df_hboxy62 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(62, start = 550, end = 800)
    df_hboxy62 %>% f_plot_hboxy(62) 
    df_hboxy62$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(62) 
    df_hbtot62 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(62, start = 400, end = 800)
    df_hbtot62 %>% f_plot_hbtot(62) 
    df_hbtot62$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 63: TM1016-01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(63)
    df_hboxy63 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(63, start = 0, end = 1150)
    df_hboxy63 %>% f_plot_hboxy(63) 
    df_hboxy63$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(63) 
    df_hbtot63 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(63, start = 0, end = 1150)
    df_hbtot63 %>% f_plot_hbtot(63) 
    df_hbtot63$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 64: TM1016UM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(64)
    df_hboxy64 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(64, start = 0, end = 1150)
    df_hboxy64 %>% f_plot_hboxy(64) 
    df_hboxy64$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(64) 
    df_hbtot64 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(64, start = 0, end = 1150)
    df_hbtot64 %>% f_plot_hbtot(64) 
    df_hbtot64$Hb_tot %>% var(na.rm = TRUE) 
    
# 65: TM1017-01
    
    ## signal too short/noisy for long enough segment
    
    
# 66: TM1017UM01
    
    ## signal too short/noisy for long enough segment
    
    
# 67: TM1018-01
    
    ## signal too short/noisy for long enough segment
    
    
# 68: TM1018UM01
    
    ## signal too short/noisy for long enough segment
    
    
# 69: TM1019-01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(69)
    df_hboxy69 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(69, start = 0, end = 1150)
    df_hboxy69 %>% f_plot_hboxy(69) 
    df_hboxy69$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(69) 
    df_hbtot69 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(69, start = 0, end = 1150)
    df_hbtot69 %>% f_plot_hbtot(69) 
    df_hbtot69$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 70: TM1019UM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(70)
    df_hboxy70 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(70, start = 0, end = 1150)
    df_hboxy70 %>% f_plot_hboxy(70) 
    df_hboxy70$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(70) 
    df_hbtot70 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(70, start = 0, end = 1150)
    df_hbtot70 %>% f_plot_hbtot(70) 
    df_hbtot70$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 71: TM1020-01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(71)
    df_hboxy71 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(71, start = 0, end = 1150)
    df_hboxy71 %>% f_plot_hboxy(71) 
    df_hboxy71$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(71) 
    df_hbtot71 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(71, start = 150, end = 1000)
    df_hbtot71 %>% f_plot_hbtot(71) 
    df_hbtot71$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 72: TM1020UM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(72)
    df_hboxy72 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(72, start = 0, end = 1150)
    df_hboxy72 %>% f_plot_hboxy(72) 
    df_hboxy72$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(72) 
    df_hbtot72 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(72, start = 150, end = 1000)
    df_hbtot72 %>% f_plot_hbtot(72) 
    df_hbtot72$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 73: TM1021-01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(73)
    df_hboxy73 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(73, start = 0, end = 1150)
    df_hboxy73 %>% f_plot_hboxy(73) 
    df_hboxy73$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(73) 
    df_hbtot73 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(73, start = 0, end = 1150)
    df_hbtot73 %>% f_plot_hbtot(73) 
    df_hbtot73$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 74: TM1021UM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(74)
    df_hboxy74 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(74, start = 0, end = 1150)
    df_hboxy74 %>% f_plot_hboxy(74) 
    df_hboxy74$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(74) 
    df_hbtot74 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(74, start = 0, end = 1150)
    df_hbtot74 %>% f_plot_hbtot(74) 
    df_hbtot74$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 75: TM1022UM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(75)
    df_hboxy75 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(75, start = 0, end = 800)
    df_hboxy75 %>% f_plot_hboxy(75) 
    df_hboxy75$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(75) 
    df_hbtot75 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(75, start = 0, end = 800)
    df_hbtot75 %>% f_plot_hbtot(75) 
    df_hbtot75$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 76: TM1023UM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(76)
    df_hboxy76 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(76, start = 0, end = 800)
    df_hboxy76 %>% f_plot_hboxy(76) 
    df_hboxy76$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(76) 
    df_hbtot76 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(76, start = 0, end = 800)
    df_hbtot76 %>% f_plot_hbtot(76) 
    df_hbtot76$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 77: TM1024UM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(77)
    df_hboxy77 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(77, start = 1310, end = 1580)
    df_hboxy77 %>% f_plot_hboxy(77) 
    df_hboxy77$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(77) 
    df_hbtot77 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(77, start = 1310, end = 1580)
    df_hbtot77 %>% f_plot_hbtot(77) 
    df_hbtot77$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 78: TM1025UM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(78)
    df_hboxy78 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(78, start = 0, end = 800)
    df_hboxy78 %>% f_plot_hboxy(78) 
    df_hboxy78$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(78) 
    df_hbtot78 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(78, start = 0, end = 800)
    df_hbtot78 %>% f_plot_hbtot(78) 
    df_hbtot78$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 79: TM1025UM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(79)
    df_hboxy79 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(79, start = 0, end = 800)
    df_hboxy79 %>% f_plot_hboxy(79) 
    df_hboxy79$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(79) 
    df_hbtot79 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(79, start = 0, end = 600)
    df_hbtot79 %>% f_plot_hbtot(79) 
    df_hbtot79$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 80: TM1027UM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(80)
    df_hboxy80 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(80, start = 150, end = 650)
    df_hboxy80 %>% f_plot_hboxy(80) 
    df_hboxy80$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(80) 
    df_hbtot80 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(80, start = 0, end = 700)
    df_hbtot80 %>% f_plot_hbtot(80) 
    df_hbtot80$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 81: TM1028UM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(81)
    df_hboxy81 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(81, start = 500, end = 2000)
    df_hboxy81 %>% f_plot_hboxy(81) 
    df_hboxy81$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(81) 
    df_hbtot81 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(81, start = 500, end = 2000)
    df_hbtot81 %>% f_plot_hbtot(81) 
    df_hbtot81$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 83: TM1030UM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(83)
    df_hboxy83 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(83, start = 100, end = 800)
    df_hboxy83 %>% f_plot_hboxy(83) 
    df_hboxy83$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(83) 
    df_hbtot83 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(83, start = 360, end = 800)
    df_hbtot83 %>% f_plot_hbtot(83) 
    df_hbtot83$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 84: TM1031UM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(84)
    df_hboxy84 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(84, start = 0, end = 800)
    df_hboxy84 %>% f_plot_hboxy(84) 
    df_hboxy84$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(84) 
    df_hbtot84 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(84, start = 0, end = 800)
    df_hbtot84 %>% f_plot_hbtot(84) 
    df_hbtot84$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 85: TM1032UM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(85)
    df_hboxy85 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(85, start = 0, end = 400)
    df_hboxy85 %>% f_plot_hboxy(85) 
    df_hboxy85$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(85) 
    df_hbtot85 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(85, start = 0, end = 400)
    df_hbtot85 %>% f_plot_hbtot(85) 
    df_hbtot85$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 86: TM1033UM01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(86)
    df_hboxy86 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(86, start = 0, end = 800)
    df_hboxy86 %>% f_plot_hboxy(86) 
    df_hboxy86$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(86) 
    df_hbtot86 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(86, start = 0, end = 700)
    df_hbtot86 %>% f_plot_hbtot(86) 
    df_hbtot86$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 87: TM2002HV01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(87)
    df_hboxy87 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(87, start = 0, end = 1150)
    df_hboxy87 %>% f_plot_hboxy(87) 
    df_hboxy87$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(87) 
    df_hbtot87 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(87, start = 0, end = 1150)
    df_hbtot87 %>% f_plot_hbtot(87) 
    df_hbtot87$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 88: TM2003HV01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(88)
    df_hboxy88 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(88, start = 0, end = 1150)
    df_hboxy88 %>% f_plot_hboxy(88) 
    df_hboxy88$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(88) 
    df_hbtot88 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(88, start = 0, end = 1150)
    df_hbtot88 %>% f_plot_hbtot(88) 
    df_hbtot88$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 89: TM2005HV01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(89)
    df_hboxy89 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(89, start = 0, end = 1150)
    df_hboxy89 %>% f_plot_hboxy(89) 
    df_hboxy89$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(89) 
    df_hbtot89 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(89, start = 0, end = 1150)
    df_hbtot89 %>% f_plot_hbtot(89) 
    df_hbtot89$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 90: TM2006HV01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(90)
    df_hboxy90 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(90, start = 500, end = 1150)
    df_hboxy90 %>% f_plot_hboxy(90) 
    df_hboxy90$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(90) 
    df_hbtot90 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(90, start = 0, end = 1150)
    df_hbtot90 %>% f_plot_hbtot(90) 
    df_hbtot90$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 91: TM2007HV01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(91)
    df_hboxy91 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(91, start = 250, end = 1150)
    df_hboxy91 %>% f_plot_hboxy(91) 
    df_hboxy91$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(91) 
    df_hbtot91 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(91, start = 250, end = 1000)
    df_hbtot91 %>% f_plot_hbtot(91) 
    df_hbtot91$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 92: TM2008HV01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(92)
    df_hboxy92 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(92, start = 0, end = 800)
    df_hboxy92 %>% f_plot_hboxy(92) 
    df_hboxy92$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(92) 
    df_hbtot92 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(92, start = 0, end = 800)
    df_hbtot92 %>% f_plot_hbtot(92) 
    df_hbtot92$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 93: TM2009HV01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(93)
    df_hboxy93 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(93, start = 0, end = 1150)
    df_hboxy93 %>% f_plot_hboxy(93) 
    df_hboxy93$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(93) 
    df_hbtot93 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(93, start = 0, end = 800)
    df_hbtot93 %>% f_plot_hbtot(93) 
    df_hbtot93$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 94: TM2010HV01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(94)
    df_hboxy94 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(94, start = 0, end = 1150)
    df_hboxy94 %>% f_plot_hboxy(94) 
    df_hboxy94$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(94) 
    df_hbtot94 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(94, start = 0, end = 800)
    df_hbtot94 %>% f_plot_hbtot(94) 
    df_hbtot94$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 95: TM2011HV01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(95)
    df_hboxy95 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(95, start = 0, end = 1150)
    df_hboxy95 %>% f_plot_hboxy(95) 
    df_hboxy95$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(95) 
    df_hbtot95 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(95, start = 0, end = 1150)
    df_hbtot95 %>% f_plot_hbtot(95) 
    df_hbtot95$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 96: TM2012HV01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(96)
    df_hboxy96 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(96, start = 0, end = 1150)
    df_hboxy96 %>% f_plot_hboxy(96) 
    df_hboxy96$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(96) 
    df_hbtot96 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(96, start = 0, end = 1150)
    df_hbtot96 %>% f_plot_hbtot(96) 
    df_hbtot96$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 97: TM2013HV01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(97)
    df_hboxy97 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(97, start = 400, end = 1150)
    df_hboxy97 %>% f_plot_hboxy(97) 
    df_hboxy97$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(97) 
    df_hbtot97 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(97, start = 400, end = 1150)
    df_hbtot97 %>% f_plot_hbtot(97) 
    df_hbtot97$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 98: TM2014HV01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(98)
    df_hboxy98 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(98, start = 150, end = 1000)
    df_hboxy98 %>% f_plot_hboxy(98) 
    df_hboxy98$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(98) 
    df_hbtot98 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(98, start = 150, end = 1000)
    df_hbtot98 %>% f_plot_hbtot(98) 
    df_hbtot98$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 99: TM2015HV01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(99)
    df_hboxy99 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(99, start = 150, end = 1150)
    df_hboxy99 %>% f_plot_hboxy(99) 
    df_hboxy99$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(99) 
    df_hbtot99 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(99, start = 150, end = 1150)
    df_hbtot99 %>% f_plot_hbtot(99) 
    df_hbtot99$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 100: TM2016HV01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(100)
    df_hboxy100 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(100, start = 0, end = 1150)
    df_hboxy100 %>% f_plot_hboxy(100) 
    df_hboxy100$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(100) 
    df_hbtot100 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(100, start = 0, end = 1150)
    df_hbtot100 %>% f_plot_hbtot(100) 
    df_hbtot100$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 101: TM2017HV01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(101)
    df_hboxy101 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(101, start = 0, end = 200)
    df_hboxy101 %>% f_plot_hboxy(101) 
    df_hboxy101$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(101) 
    df_hbtot101 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(101, start = 0, end = 200)
    df_hbtot101 %>% f_plot_hbtot(101) 
    df_hbtot101$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 102: TM2018HV01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(102)
    df_hboxy102 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(102, start = 0, end = 1150)
    df_hboxy102 %>% f_plot_hboxy(102) 
    df_hboxy102$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(102) 
    df_hbtot102 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(102, start = 0, end = 1150)
    df_hbtot102 %>% f_plot_hbtot(102) 
    df_hbtot102$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 103: TM2019HV01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(103)
    df_hboxy103 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(103, start = 0, end = 1150)
    df_hboxy103 %>% f_plot_hboxy(103) 
    df_hboxy103$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(103) 
    df_hbtot103 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(103, start = 0, end = 1150)
    df_hbtot103 %>% f_plot_hbtot(103) 
    df_hbtot103$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 104: TM2020HV01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(104)
    df_hboxy104 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(104, start = 0, end = 1150)
    df_hboxy104 %>% f_plot_hboxy(104) 
    df_hboxy104$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(104) 
    df_hbtot104 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(104, start = 0, end = 1150)
    df_hbtot104 %>% f_plot_hbtot(104) 
    df_hbtot104$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 105: TM2021HV01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(105)
    df_hboxy105 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(105, start = 100, end = 1150)
    df_hboxy105 %>% f_plot_hboxy(105) 
    df_hboxy105$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(105) 
    df_hbtot105 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(105, start = 0, end = 600)
    df_hbtot105 %>% f_plot_hbtot(105) 
    df_hbtot105$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 106: TM2022HV01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(106)
    df_hboxy106 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(106, start = 200, end = 1150)
    df_hboxy106 %>% f_plot_hboxy(106) 
    df_hboxy106$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(106) 
    df_hbtot106 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(106, start = 0, end = 1000)
    df_hbtot106 %>% f_plot_hbtot(106) 
    df_hbtot106$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 107: TM2023HV01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(107)
    df_hboxy107 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(107, start = 0, end = 1150)
    df_hboxy107 %>% f_plot_hboxy(107) 
    df_hboxy107$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(107) 
    df_hbtot107 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(107, start = 0, end = 1150)
    df_hbtot107 %>% f_plot_hbtot(107) 
    df_hbtot107$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 108: TM2024-01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(108)
    df_hboxy108 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(108, start = 250, end = 1000)
    df_hboxy108 %>% f_plot_hboxy(108) 
    df_hboxy108$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(108) 
    df_hbtot108 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(108, start = 0, end = 1000)
    df_hbtot108 %>% f_plot_hbtot(108) 
    df_hbtot108$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 109: TM2024HV01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(109)
    df_hboxy109 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(109, start = 250, end = 1000)
    df_hboxy109 %>% f_plot_hboxy(109) 
    df_hboxy109$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(109) 
    df_hbtot109 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(109, start = 0, end = 1000)
    df_hbtot109 %>% f_plot_hbtot(109) 
    df_hbtot109$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 110: TM2025-01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(110)
    df_hboxy110 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(110, start = 0, end = 1000)
    df_hboxy110 %>% f_plot_hboxy(110) 
    df_hboxy110$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(110) 
    df_hbtot110 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(110, start = 0, end = 1000)
    df_hbtot110 %>% f_plot_hbtot(110) 
    df_hbtot110$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 111: TM2025HV01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(111)
    df_hboxy111 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(111, start = 0, end = 1000)
    df_hboxy111 %>% f_plot_hboxy(111) 
    df_hboxy111$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(111) 
    df_hbtot111 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(111, start = 0, end = 1000)
    df_hbtot111 %>% f_plot_hbtot(111) 
    df_hbtot111$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 112: TM2026-01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(112)
    df_hboxy112 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(112, start = 0, end = 1000)
    df_hboxy112 %>% f_plot_hboxy(112) 
    df_hboxy112$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(112) 
    df_hbtot112 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(112, start = 0, end = 1000)
    df_hbtot112 %>% f_plot_hbtot(112) 
    df_hbtot112$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 113: TM2026HV01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(113)
    df_hboxy113 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(113, start = 0, end = 1000)
    df_hboxy113 %>% f_plot_hboxy(113) 
    df_hboxy113$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(113) 
    df_hbtot113 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(113, start = 0, end = 1000)
    df_hbtot113 %>% f_plot_hbtot(113) 
    df_hbtot113$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 114: TM2027-01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(114)
    df_hboxy114 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(114, start = 0, end = 1150)
    df_hboxy114 %>% f_plot_hboxy(114) 
    df_hboxy114$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(114) 
    df_hbtot114 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(114, start = 0, end = 1150)
    df_hbtot114 %>% f_plot_hbtot(114) 
    df_hbtot114$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 115: TM2027HV01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(115)
    df_hboxy115 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(115, start = 200, end = 1150)
    df_hboxy115 %>% f_plot_hboxy(115) 
    df_hboxy115$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(115) 
    df_hbtot115 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(115, start = 0, end = 1150)
    df_hbtot115 %>% f_plot_hbtot(115) 
    df_hbtot115$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 116: TM2028-01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(116)
    df_hboxy116 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(116, start = 0, end = 350)
    df_hboxy116 %>% f_plot_hboxy(116) 
    df_hboxy116$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(116) 
    df_hbtot116 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(116, start = 0, end = 300)
    df_hbtot116 %>% f_plot_hbtot(116) 
    df_hbtot116$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 117: TM2028HV01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(117)
    df_hboxy117 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(117, start = 0, end = 350)
    df_hboxy117 %>% f_plot_hboxy(117) 
    df_hboxy117$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(117) 
    df_hbtot117 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(117, start = 0, end = 300)
    df_hbtot117 %>% f_plot_hbtot(117) 
    df_hbtot117$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 118: TM2029-01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(118)
    df_hboxy118 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(118, start = 500, end = 1000)
    df_hboxy118 %>% f_plot_hboxy(118) 
    df_hboxy118$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(118) 
    df_hbtot118 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(118, start = 500, end = 1000)
    df_hbtot118 %>% f_plot_hbtot(118) 
    df_hbtot118$Hb_tot %>% var(na.rm = TRUE) 
    
    
# 119: TM2029HV01
    
    # Hb_oxy 
    df_raw_data %>% f_plot_hboxy(119)
    df_hboxy119 <- df_raw_data %>% dplyr::select(-Hb_tot) %>% f_segment_signal(119, start = 500, end = 1000)
    df_hboxy119 %>% f_plot_hboxy(119) 
    df_hboxy119$Hb_oxy %>% var(na.rm = TRUE)
    
    # Hb_tot
    df_raw_data %>% f_plot_hbtot(119) 
    df_hbtot119 <- df_raw_data %>% dplyr::select(-Hb_oxy) %>% f_segment_signal(119, start = 500, end = 1000)
    df_hbtot119 %>% f_plot_hbtot(119) 
    df_hbtot119$Hb_tot %>% var(na.rm = TRUE) 
    
    
# Once you have repeated this for all the patients, combine the filtered Hb_tot signals into one tibble and the filtered Hb_oxy signals into another
    
    df_hbtot_filt <- bind_rows(df_hbtot1, df_hbtot2, df_hbtot3, df_hbtot4, df_hbtot5, df_hbtot6, df_hbtot7, df_hbtot8, df_hbtot9, df_hbtot10,
                               df_hbtot11, df_hbtot12, df_hbtot13, df_hbtot14, df_hbtot15, df_hbtot16, df_hbtot17, df_hbtot18, df_hbtot20,
                               df_hbtot21, df_hbtot22, df_hbtot26, df_hbtot27, df_hbtot28, df_hbtot29, df_hbtot30,
                               df_hbtot31, df_hbtot32, df_hbtot33, df_hbtot34, df_hbtot35, df_hbtot36, df_hbtot37, df_hbtot38, df_hbtot39, df_hbtot40,
                               df_hbtot41, df_hbtot42, df_hbtot43, df_hbtot44, df_hbtot45, df_hbtot46, df_hbtot47, df_hbtot48, df_hbtot49, df_hbtot50,
                               df_hbtot51, df_hbtot52, df_hbtot53, df_hbtot54, df_hbtot55, df_hbtot58, df_hbtot59, df_hbtot60,
                               df_hbtot61, df_hbtot62, df_hbtot63, df_hbtot64, df_hbtot69, df_hbtot70,
                               df_hbtot71, df_hbtot72, df_hbtot73, df_hbtot74, df_hbtot75, df_hbtot76, df_hbtot77, df_hbtot78, df_hbtot79, df_hbtot80,
                               df_hbtot81, df_hbtot83, df_hbtot84, df_hbtot85, df_hbtot86, df_hbtot87, df_hbtot88, df_hbtot89, df_hbtot90,
                               df_hbtot91, df_hbtot92, df_hbtot93, df_hbtot94, df_hbtot95, df_hbtot96, df_hbtot97, df_hbtot98, df_hbtot99, df_hbtot100,
                               df_hbtot101, df_hbtot102, df_hbtot103, df_hbtot104, df_hbtot105, df_hbtot106, df_hbtot107, df_hbtot108, df_hbtot109, df_hbtot110,
                               df_hbtot111, df_hbtot112, df_hbtot113, df_hbtot114, df_hbtot115, df_hbtot116, df_hbtot117, df_hbtot118, df_hbtot119)
    
    df_hboxy_filt <- bind_rows(df_hboxy1, df_hboxy2, df_hboxy3, df_hboxy4, df_hboxy5, df_hboxy6, df_hboxy7, df_hboxy8, df_hboxy9, df_hboxy10,
                               df_hboxy11, df_hboxy12, df_hboxy13, df_hboxy14, df_hboxy15, df_hboxy16, df_hboxy17, df_hboxy18, df_hboxy20,
                               df_hboxy21, df_hboxy22, df_hboxy26, df_hboxy27, df_hboxy28, df_hboxy29, df_hboxy30,
                               df_hboxy31, df_hboxy32, df_hboxy33, df_hboxy34, df_hboxy35, df_hboxy36, df_hboxy37, df_hboxy38, df_hboxy39, df_hboxy40,
                               df_hboxy41, df_hboxy42, df_hboxy43, df_hboxy44, df_hboxy45, df_hboxy46, df_hboxy47, df_hboxy48, df_hboxy49, df_hboxy50,
                               df_hboxy51, df_hboxy52, df_hboxy53, df_hboxy54, df_hboxy55, df_hboxy58, df_hboxy59, df_hboxy60,
                               df_hboxy61, df_hboxy62, df_hboxy63, df_hboxy64, df_hboxy69, df_hboxy70,
                               df_hboxy71, df_hboxy72, df_hboxy73, df_hboxy74, df_hboxy75, df_hboxy76, df_hboxy77, df_hboxy78, df_hboxy79, df_hboxy80,
                               df_hboxy81, df_hboxy83, df_hboxy84, df_hboxy85, df_hboxy86, df_hboxy87, df_hboxy88, df_hboxy89, df_hboxy90,
                               df_hboxy91, df_hboxy92, df_hboxy93, df_hboxy94, df_hboxy95, df_hboxy96, df_hboxy97, df_hboxy98, df_hboxy99, df_hboxy100,
                               df_hboxy101, df_hboxy102, df_hboxy103, df_hboxy104, df_hboxy105, df_hboxy106, df_hboxy107, df_hboxy108, df_hboxy109, df_hboxy110,
                               df_hboxy111, df_hboxy112, df_hboxy113, df_hboxy114, df_hboxy115, df_hboxy116, df_hboxy117, df_hboxy118, df_hboxy119)
    
    

    

# Pre-process segmented signals -------------------------------------------

    
    
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
      ggtitle("Muscle hemoglobin concentration") +
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
      ggtitle("Muscle hemoglobin oxygen saturation") +
      theme(legend.position = "none")
    
    # Explore any nonsensical values, especially any negatives in the Hb_oxy. you may want to remove those patients from the analysis due to unreliable readings
    
    
    
### SAVITZY-GOLAY SIGNAL SMOOTHING
    
    
    # For some reason, the Oxiplex removes milliseconds after a few seconds. Add milliseconds back in to Time
    
    # Hb_tot
    df_hbtot_filt <- df_hbtot_filt %>% 
      group_by(number, Time) %>% 
      left_join(count(.)) %>% 
      mutate(Time = ifelse(n == 49, Time + 0.02*row_number(), Time + 0.02*(row_number() - 1))) %>% 
      dplyr::select(-n)
    
    # Hb_oxy
    df_hboxy_filt <- df_hboxy_filt %>% 
      group_by(number, Time) %>% 
      left_join(count(.)) %>% 
      mutate(Time = ifelse(n == 49, Time + 0.02*row_number(), Time + 0.02*(row_number() - 1))) %>% 
      dplyr::select(-n)
    
    
    # 4th order savitzy-golay filter
    
    library(signal)
    
    # Hb_tot
    df_hbtot_filt <- df_hbtot_filt %>% 
      group_by(number) %>% 
      mutate(sg_filt = sgolayfilt(Hb_tot, p = 4)) %>% 
      ungroup()
    
    # plot to visualize smoothing
    df_hbtot_filt %>% dplyr::filter(number == 10) %>% # can pick any arbitrary time segment to look at
      ggplot(aes(x = Time)) +
      geom_line(aes(y = Hb_tot)) +
      geom_line(aes(y = sg_filt), color = "red")
    
    # Hb_oxy 
    df_hboxy_filt <- df_hboxy_filt %>% 
      group_by(number) %>% 
      mutate(sg_filt = sgolayfilt(Hb_oxy, p = 3)) %>% 
      ungroup()
    
    # plot
    df_hboxy_filt %>% dplyr::filter(number == 10 & Time < 1500 & Time > 1490) %>% 
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
    l_hbtot_filt <- df_hbtot_filt %>% split(df_hbtot_filt$number) # split filtered signals into list of tibbles
    sampling_rate <- 50 # again, may need to adjust with different sampling rates
    butter_filter <- butter(4, W = 0.01/(sampling_rate/2), type = "high") # create butter filter
    l_hbtot_filt_pass <- 1:length(l_hbtot_filt) %>% purrr::map(~mutate(l_hbtot_filt[[.x]], high_pass = filtfilt(butter_filter, sg_filt))) # apply to SG filtered signal
    
    df_hbtot_filt_pass <- rbindlist(l_hbtot_filt_pass) %>% as_tibble() # combing again into tibble
    
    # repeat for Hb_oxy
    l_hboxy_filt <- df_hboxy_filt %>% split(df_hboxy_filt$number)
    l_hboxy_filt_pass <- 1:length(l_hboxy_filt) %>% purrr::map(~mutate(l_hboxy_filt[[.x]], high_pass = filtfilt(butter_filter, sg_filt)))
    
    df_hboxy_filt_pass <- rbindlist(l_hboxy_filt_pass) %>% as_tibble()
    
    
    
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
    
    
    f_plot_fft(df_hboxy_filt_pass, 10, zoom = FALSE)  
    f_plot_fft(df_hbtot_filt_pass, 2, zoom = TRUE)  
    
    
    
### SAVE THE FILTERED SIGNALS
    
    # This is important so you don't have to repeat all the signal pre-processing before running the DFA
    
    df_hbtot_filt_muscle <- df_hbtot_filt
    df_hbtot_filt_pass_muscle <- df_hbtot_filt_pass
    df_hboxy_filt_muscle <- df_hboxy_filt
    df_hboxy_filt_pass_muscle <- df_hboxy_filt_pass
    
    # e.g.
    save(df_hbtot_filt_muscle, df_hbtot_filt_pass_muscle, 
         df_hboxy_filt_muscle, df_hboxy_filt_pass_muscle,
         file = "filtered_muscle_signals.Rdata")

    

    
    
        
# Detrended fluctuation analysis ------------------------------------------
    
    
### LOAD IN DATA
    
    # If not already in your environment, load the processed signals
    
    load("filtered_muscle_signals.Rdata")
    
    
    
### COMPUTE CUMULATIVE SUM OF EACH SIGNAL
    
    
    # DFA is performed on the cumulative sum of the signal, not the signal itself
    
  ## Hb_tot
    
    # create list, where each object is a patient
    l_hbtot_filt_pass_muscle <- df_hbtot_filt_pass_muscle %>% split(f = df_hbtot_filt_pass_muscle$number)
    
    # take the mean of each filtered signal
    l_means_muscle <- 1:length(l_hbtot_filt_pass_muscle) %>% purrr::map(~mean(l_hbtot_filt_pass_muscle[[.x]]$sg_filt, na.rm = TRUE))
    
    # subtract the mean and take the cumulative sum
    l_hbtot_cumsum_muscle <- 1:length(l_means_muscle) %>% purrr::map(~mutate(l_hbtot_filt_pass_muscle[[.x]], 
                                                                             mean = l_means_muscle[[.x]]))
    df_hbtot_cumsum_muscle <- rbindlist(l_hbtot_cumsum_muscle) %>% as_tibble()
    
    df_hbtot_cumsum_muscle <- df_hbtot_cumsum_muscle %>% 
      mutate(subtract_mean = sg_filt - mean) %>% 
      group_by(number) %>% 
      mutate(cumsum = cumsum(subtract_mean)) %>% 
      ungroup()
    
    # plot
    df_hbtot_cumsum_muscle %>% dplyr::filter(number == 1) %>% 
      ggplot(aes(x = Time)) +
      geom_line(aes(y = cumsum)) +
      geom_line(aes(y = sg_filt), color = "red") +
      
      theme_bw() +
      ggtitle("Cumulative sum of Savitzy-Golar smoothed signal")
    
    
  ## Hb_oxy
    
    # create list, where each object is a patient
    l_hboxy_filt_pass_muscle <- df_hboxy_filt_pass_muscle %>% split(f = df_hboxy_filt_pass_muscle$number)
    
    # take the mean of each filtered signal
    l_means_muscle <- 1:length(l_hboxy_filt_pass_muscle) %>% purrr::map(~mean(l_hboxy_filt_pass_muscle[[.x]]$sg_filt, na.rm = TRUE))
    
    # subtract the mean and take the cumulative sum
    l_hboxy_cumsum_muscle <- 1:length(l_means_muscle) %>% purrr::map(~mutate(l_hboxy_filt_pass_muscle[[.x]], 
                                                                             mean = l_means_muscle[[.x]]))
    df_hboxy_cumsum_muscle <- rbindlist(l_hboxy_cumsum_muscle) %>% as_tibble()
    
    df_hboxy_cumsum_muscle <- df_hboxy_cumsum_muscle %>% 
      mutate(subtract_mean = sg_filt - mean) %>% 
      group_by(number) %>% 
      mutate(cumsum = cumsum(subtract_mean)) %>% 
      ungroup()
    
    # plot
    df_hboxy_cumsum_muscle %>% dplyr::filter(number == 1) %>% 
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
    l_hbtot_cumsum_muscle <- df_hbtot_cumsum_muscle %>% split(f = df_hbtot_cumsum_muscle$number)
    l_hboxy_cumsum_muscle <- df_hboxy_cumsum_muscle %>% split(f = df_hboxy_cumsum_muscle$number)
    
    
    # repeat for all patients!!
    
    library(doParallel) # runs faster in parallel
    registerDoParallel()
    
    
    l_dfa_hbtot_muscle <- 1:length(l_hbtot_cumsum_muscle) %>% purrr::map(~f_dfa(l_hbtot_cumsum_muscle[[.x]], num_of_windows = 10))
    l_dfa_hboxy_muscle <- 1:length(l_hboxy_cumsum_muscle) %>% purrr::map(~f_dfa(l_hboxy_cumsum_muscle[[.x]], num_of_windows = 10))
    
    
    # would be wise to save DFA results here because it takes a while to run
    save(l_dfa_hbtot_muscle, l_dfa_hboxy_muscle, file = "muscle_dfa_results.Rdata")   
    
    
    
    
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
    hbtot_alphas_muscle <- 1:length(l_dfa_hbtot_muscle) %>% purrr::map(~l_dfa_hbtot_muscle[[.x]][[2]]) %>% unlist()
    
    l_hbtot_second_alpha_muscle <- 1:length(l_dfa_hbtot_muscle) %>% purrr::map(~f_segment_slope(l_dfa_hbtot_muscle[[.x]]))
    hbtot_second_alpha_muscle <- 1:length(l_dfa_hbtot_muscle) %>% purrr::map(~f_segment_slope(l_dfa_hbtot_muscle[[.x]])) %>% unlist()
    
    # combine dfa results into df
    df_hbtot_dfa_muscle <- df_hbtot_cumsum_muscle %>% 
      count(number, Status) %>% 
      dplyr::select(-n) %>% 
      mutate(alpha = hbtot_alphas_muscle, 
             second_alpha = hbtot_second_alpha_muscle,
             Status = as.factor(Status))
    
    # Hb_oxy
    hboxy_alphas_muscle <- 1:length(l_dfa_hboxy_muscle) %>% purrr::map(~l_dfa_hboxy_muscle[[.x]][[2]]) %>% unlist()
    
    l_hboxy_second_alpha_muscle <- 1:length(l_dfa_hboxy_muscle) %>% purrr::map(~f_segment_slope(l_dfa_hboxy_muscle[[.x]]))
    hboxy_second_alpha_muscle <- 1:length(l_dfa_hboxy_muscle) %>% purrr::map(~f_segment_slope(l_dfa_hboxy_muscle[[.x]])) %>% unlist()
    
    # combine dfa results into df
    df_hboxy_dfa_muscle <- df_hboxy_cumsum_muscle %>% 
      count(number, Status) %>% 
      dplyr::select(-n) %>% 
      mutate(alpha = hboxy_alphas_muscle, 
             second_alpha = hboxy_second_alpha_muscle,
             Status = as.factor(Status))

    
    
    
    
    
# plots and statistical analyses ------------------------------------------
    
    
    library(DescTools)
    
    # plot alphas for Hb_tot and Hb_oxy
    
    
    # Hb_tot
    
    # plot
    df_hbtot_dfa_muscle %>% 
      ggplot(aes(x = Status, y = alpha)) +
      geom_violin(aes(color = Status)) +
      geom_jitter(position = position_jitter(0.2), shape = 1) +
      stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
      
      labs(y = "alpha") +
      ggtitle("Muscle hemoglobin concentration alpha") +
      
      theme_bw() +
      theme(legend.position = "none")
    
    # test for differences
    df_hbtot_dfa_muscle %>% group_by(Status) %>% summarise(median = median(alpha))
    kruskal.test(alpha ~ Status, data = df_hbtot_dfa_muscle)
    DunnTest(alpha ~ Status, data = df_hbtot_dfa_muscle)
    
    
    
    # Hb_oxy
    
    # plot
    df_hboxy_dfa_muscle %>% 
      ggplot(aes(x = Status, y = alpha)) +
      geom_violin(aes(color = Status)) +
      geom_jitter(position = position_jitter(0.2), shape = 1) +
      stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
      
      labs(y = "alpha") +
      ggtitle("Muscle hemoglobin oxygen saturation alpha") +
      
      theme_bw() +
      theme(legend.position = "none")
    
    
    # test for differences
    df_hboxy_dfa_muscle %>% group_by(Status) %>% summarise(median = median(alpha))
    kruskal.test(alpha ~ Status, data = df_hboxy_dfa_muscle)
    DunnTest(alpha ~ Status, data = df_hboxy_dfa_muscle)
    
    
    
    

# combine results for export ----------------------------------------------


    # pair number with subject_id, remove duplicates (ending in 2)
    df_num_id_muscle <- df_hbtot_cumsum_muscle %>% 
      count(number, subject_id) %>% 
      dplyr::select(-n)
    
    df_num_id_muscle <- df_num_id_muscle %>% dplyr::filter(!grepl("01 2", subject_id))
    
    # Take means
    df_hbtot_muscle_mean <- df_hbtot_filt_muscle %>% 
      dplyr::filter(!is.na(sg_filt)) %>% 
      group_by(number, Status) %>% 
      summarise(avg_Hb_conc_muscle = mean(sg_filt))
    
    df_hboxy_muscle_mean <- df_hboxy_filt_muscle %>% 
      dplyr::filter(!is.na(sg_filt)) %>% 
      group_by(number, Status) %>% 
      summarise(avg_Hb_o2sat_muscle = mean(sg_filt))
    
    # combine into df
    
    df_muscle_results <- df_num_id_muscle %>% 
      left_join(df_hbtot_dfa_muscle, by = "number") %>% 
      dplyr::rename("muscle_Hb_conc_alpha" = "alpha", "muscle_Hb_conc_second_alpha" = "second_alpha") %>% 
      left_join(df_hboxy_dfa_muscle, by = c("number", "Status")) %>% 
      rename("muscle_Hb_o2sat_alpha" = "alpha", "muscle_Hb_o2sat_second_alpha" = "second_alpha") %>% 
      
      # add hbtot and hboxy means
      left_join(df_hbtot_muscle_mean, by = "number") %>% 
      left_join(df_hboxy_muscle_mean, by = "number") %>% 
      
      # reorder
      dplyr::select(number, subject_id, Status, 
                    avg_Hb_conc_muscle, muscle_Hb_conc_alpha, muscle_Hb_conc_second_alpha,
                    avg_Hb_o2sat_muscle, muscle_Hb_o2sat_alpha, muscle_Hb_o2sat_second_alpha)
    
    # save!!!
    save(df_muscle_results, file = "muscle_dfa_results.Rdata")
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    






