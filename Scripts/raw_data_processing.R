


# read in data -------------------------------------------------------------------------

  library(tidyverse)
  library(purrr)
  

  # get names
  my_files <- list.files(path = "Data/Raw data/", pattern = "*.txt", full.names = TRUE, recursive = FALSE)
  
  # remove duplicate patient files (CM 14, 15, 16, 24part2)
  my_files <- my_files[!grepl("part2.txt|TM0014CM01.txt|TM0015CM01.txt|TM0016CM01.txt", my_files)]
  
  # read in
  l_raw_data <- 1:length(my_files) %>% 
    purrr::map(~read.delim2(file = my_files[.x], header = FALSE, sep = "\t") %>% as_tibble()) 
  
  # name with Status
  l_names <- 1:length(my_files) %>% 
    purrr::map(~substr(my_files[.x], start = 16, stop = 25)) 
  
  names(l_raw_data) <- l_names


  
# format data -------------------------------------------------------------

  library(data.table)
  library(zoo)
  
  # format to just THC and o2 sat in workable form
  l_raw_data <- 1:length(l_raw_data) %>% 
    
    purrr::map(~mutate(l_raw_data[[.x]], subject_id = l_names[[.x]])) %>% 
    purrr::map(~dplyr::filter(.x, !grepl("Patient|Data|[R|r]aw", V1))) %>% 
    purrr::map(~dplyr::select(.x, c(ncol(.x), 2:4))) %>% 
    purrr::map(~dplyr::rename(.x, "Time" = "V2", "O2_sat" = "V3", "THC" = "V4")) %>% 
    purrr::map(~mutate(.x, Time = as.numeric(Time))) %>% 
    purrr::map(~mutate(.x, THC = as.numeric(THC))) %>% 
    purrr::map(~mutate(.x, O2_sat = as.numeric(O2_sat)))
  
  # renumber patients
  l_raw_data <- 1:length(l_raw_data) %>% 
    purrr::map(~mutate(l_raw_data[[.x]], number = c(1:length(l_raw_data))[.x])) 
  
  # condense list to df
  df_raw_data <- rbindlist(l_raw_data) %>% as_tibble() 
  
  # add status
  df_raw_data <- df_raw_data %>% 
    mutate(Status = substr(subject_id, start = 7, stop = 8)) %>% 
    dplyr::select(c(number, subject_id, Status, everything()))
  
  # make sure status is HV, UM, or CM
  df_raw_data$Status[df_raw_data$Status == "HC"] <- "HV"
  df_raw_data <- df_raw_data %>% dplyr::filter(Status != "-0") # unclear what group this patient was in
  
  # # smooth data with 1s moving average
  # df_thc <- df_thc %>% 
  #   group_by(number) %>% 
  #   mutate(roll_mean = zoo::rollmean(THC, k = 50, fill = NA)) %>% 
  #   ungroup()
  

  save(df_raw_data, file = "raw_data.Rdata")


  
# explore --------------------------------------------------------------------


  # distribution of THC data by status only
  df_raw_data %>% 
    
    ggplot(aes(THC)) +
    geom_histogram() +
    facet_wrap(~Status, scales = "free")
  
  # distribution by patient number
  df_raw_data %>% 
    
    ggplot(aes(THC)) +
    geom_histogram() +
    facet_grid(~number, scales = "free") 

  # plot all signals together
  df_raw_data %>% 
    
    ggplot(aes(x = Time, y = THC)) +
    geom_line() +
    facet_wrap(~Status, scales = "free_y") +
    
    theme_bw()
  
  # facet by number
  df_raw_data %>% 
    
    ggplot(aes(x = Time, y = THC)) +
    geom_line() +
    facet_wrap(~number, scales = "free_y") +
    
    theme_bw()
  
  


# remove blips (manual) ------------------------------------------------------------


  load("raw_data.Rdata")
  
  
# remove signals that aren't at least 200s  
  sampling_rate <- 1/.02
  df_raw_data %>% group_by(number) %>% count() %>% dplyr::filter(n < 200*sampling_rate) # 21 and 26 not long enough
  df_raw_data <- df_raw_data %>% filter(!(number %in% c(21, 26)))
  
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

  
# Plot all together
  
  df_raw_data_long <- df_raw_data %>% 
    pivot_longer(5:6, values_to = "value", names_to = "variable")
  
  f_plot_all <- function(df, num) {
    
    subject_id <- df %>% 
      count(number, subject_id) %>% 
      dplyr::filter(number == num) %>% 
      .$subject_id 
    
    df %>% 
      dplyr::filter(number == num) %>% 
      ggplot(aes(x = Time, y = value)) +
      geom_line() +
      facet_wrap(~variable, scales = "free_y") +
      ggtitle(subject_id) +
      theme_bw()
    
  }
  
  f_plot_all(df_raw_data_long, 2)

  # save as pdf
  l_plot <- list()
  for (i in unique(df_raw_data_long$number)) {
    
    p <- f_plot_all(df_raw_data_long, i)
    l_plot[[i]] <- p
    
  }
  
  pdf("Outputs/brain_tracings.pdf")
  for (i in unique(df_raw_data_long$number)) {
    print(l_plot[[i]])    
  }
  dev.off()
  
# 1
  # THC
  df_raw_data %>% f_plot_thc(1) + theme_bw()# eyeball where after blips is
  df_thc1 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(1, 1200, 2500) # shorten signal to that length, must be >200s
  df_thc1$THC %>% var(na.rm = TRUE) # check var less than 1000
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(1)
  df_o2sat1 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(1, 1200, 2500)
  df_o2sat1$O2_sat %>% var(na.rm = TRUE)
  

# 2
  # THC
  df_raw_data %>% f_plot_thc(2) 
  df_thc2 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(2, 900, 2000)
  df_thc2$THC %>% var(na.rm = TRUE)

  # O2 sat
  df_raw_data %>% f_plot_o2sat(2)
  df_o2sat2 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(2, 900, 2000)
  df_o2sat2$O2_sat %>% var(na.rm = TRUE)
  
  
# 3
  # THC
  df_raw_data %>% f_plot_thc(3) 
  df_thc3 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(3, 500, 2000)
  df_thc3$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(3)
  df_o2sat3 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(3, 800, 2000)
  df_o2sat3$O2_sat %>% var(na.rm = TRUE)
  

# 4
  # THC
  df_raw_data %>% f_plot_thc(4) 
  df_thc4 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(4, 0, 2000)
  df_thc4$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(4)
  df_o2sat4 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(4, 0, 2000)
  df_o2sat4$O2_sat %>% var(na.rm = TRUE)
  
  
# 5
  # THC
  df_raw_data %>% f_plot_thc(5) # ***same as 4???
  df_thc5 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(5, 0, 2000)
  df_thc5$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(5)
  df_o2sat5 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(5, 0, 2000)
  df_o2sat5$O2_sat %>% var(na.rm = TRUE)
  
  
# 6
  # THC
  df_raw_data %>% f_plot_thc(6) # same as 4 & 5
  df_thc6 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(6, 0, 2000) 
  df_thc6$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(6)
  df_o2sat6 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(6, 0, 2000)
  df_o2sat6$O2_sat %>% var(na.rm = TRUE)
  
  
# 7
  # THC
  df_raw_data %>% f_plot_thc(7)
  df_thc7 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(7, 0, 2000) 
  df_thc7$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(7)
  df_o2sat7 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(7, 0, 2000)
  df_o2sat7$O2_sat %>% var(na.rm = TRUE)
  
  
# 8
  # THC
  df_raw_data %>% f_plot_thc(8)
  df_thc8 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(8, 0, 2000)
  df_thc8$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(8)
  df_o2sat8 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(8, 0, 2000)
  df_o2sat8$O2_sat %>% var(na.rm = TRUE)
  
  
# 9
  # THC
  df_raw_data %>% f_plot_thc(9)
  df_thc9 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(9, 0, 2000) 
  df_thc9$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(9)
  df_o2sat9 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(9, 0, 2000)
  df_o2sat9$O2_sat %>% var(na.rm = TRUE)
  
  
# 10
  # THC
  df_raw_data %>% f_plot_thc(10)
  df_thc10 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(10, 0, 2000) 
  df_thc10$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(10)
  df_o2sat10 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(10, 0, 2000)
  df_o2sat10$O2_sat %>% var(na.rm = TRUE)
  
  
# 11
  # THC
  df_raw_data %>% f_plot_thc(11)
  df_thc11 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(11, 0, 2000)
  df_thc11$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(11)
  df_o2sat11 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(11, 0, 2000)
  df_o2sat11$O2_sat %>% var(na.rm = TRUE)
  
  
# 12
  # THC
  df_raw_data %>% f_plot_thc(12)
  df_thc12 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(12, 0, 2000) 
  df_thc12$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(12)
  df_o2sat12 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(12, 0, 1800)
  df_o2sat12$O2_sat %>% var(na.rm = TRUE)
  

# 13
  # THC
  df_thc13 %>% f_plot_thc(13)
  df_thc13 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(13, 300, 1500)
  df_thc13$THC %>% var(na.rm = TRUE) 
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(13)
  df_o2sat13 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(13, 875, 1190)
  df_o2sat13$O2_sat %>% var(na.rm = TRUE)
  
  
# 14
  # THC
  df_raw_data %>% f_plot_thc(14)
  df_thc14 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(14, 400, 2000)
  df_thc14$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(14) ### signal too noisy

  
# 15
  # THC
  df_raw_data %>% f_plot_thc(15)
  df_thc15 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(15, 0, 2000)
  df_thc15$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(15)
  df_o2sat15 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(15, 0, 2000)
  df_o2sat15$O2_sat %>% var(na.rm = TRUE)
  
  
# 16
  # THC
  df_raw_data %>% f_plot_thc(16)
  df_thc16 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(16, 0, 2000)
  df_thc16$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(16)
  df_o2sat16 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(16, 0, 2000)
  df_o2sat16$O2_sat %>% var(na.rm = TRUE)
  

# 17
  # THC
  df_raw_data %>% f_plot_thc(17)
  df_thc17 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(17, 0, 2000)
  df_thc17$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(17)
  df_o2sat17 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(17, 0, 2000)
  df_o2sat17$O2_sat %>% var(na.rm = TRUE)
  
  
# 18
  # THC
  df_raw_data %>% f_plot_thc(18)
  df_thc18 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(18, 0, 2000)
  df_thc18$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(18)
  df_o2sat18 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(18, 1250, 2000)
  df_o2sat18$O2_sat %>% var(na.rm = TRUE)
  
  
# 19
  # THC
  df_raw_data %>% f_plot_thc(19)
  df_thc19 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(19, 0, 2000)
  df_thc19$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(19)
  df_o2sat19 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(19, 0, 2000)
  df_o2sat19$O2_sat %>% var(na.rm = TRUE)
  
  
# 20
  # THC
  df_raw_data %>% f_plot_thc(20)
  df_thc20 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(20, 1350, 2000)
  df_thc20$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(20)
  df_o2sat20 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(20, 1350, 2000)
  df_o2sat20$O2_sat %>% var(na.rm = TRUE)
  
  
# 22
  # THC
  df_raw_data %>% f_plot_thc(22)
  df_thc22 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(22, 0, 2000) 
  df_thc22$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(22)
  df_o2sat22 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(22, 0, 1000)
  df_o2sat22$O2_sat %>% var(na.rm = TRUE)
  
  
# 23
  # THC
  df_raw_data %>% f_plot_thc(23)
  df_thc23 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(23, 0, 2000) 
  df_thc23$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(23)
  df_o2sat23 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(23, 0, 2000)
  df_o2sat23$O2_sat %>% var(na.rm = TRUE)
  
  
# 24
  # THC
  df_raw_data %>% f_plot_thc(24)
  df_thc24 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(24, 0, 2000)
  df_thc24$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(24)
  df_o2sat24 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(24, 0, 2000)
  df_o2sat24$O2_sat %>% var(na.rm = TRUE)
  
  
# 27
  # THC
  df_raw_data %>% f_plot_thc(27)
  df_thc27 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(27, 0, 2000) 
  df_thc27$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(27)
  df_o2sat27 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(27, 0, 2000)
  df_o2sat27$O2_sat %>% var(na.rm = TRUE)
  
  
# 28
  # THC
  df_raw_data %>% f_plot_thc(28)
  df_thc28 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(28, 0, 2000) 
  df_thc28$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(28)
  df_o2sat28 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(28, 0, 2000)
  df_o2sat28$O2_sat %>% var(na.rm = TRUE)
  
  
# 29
  # THC
  df_raw_data %>% f_plot_thc(29)
  df_thc29 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(29, 0, 2000)
  df_thc29$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(29)
  df_o2sat29 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(29, 0, 2000)
  df_o2sat29$O2_sat %>% var(na.rm = TRUE)
  
  
# 30
  # THC
  df_raw_data %>% f_plot_thc(30)
  df_thc30 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(30, 0, 400) # very noisy signal
  df_thc30$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(30)
  df_o2sat30 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(30, 0, 400)
  df_o2sat30$O2_sat %>% var(na.rm = TRUE)
  
  
# 31
  # THC
  df_raw_data %>% f_plot_thc(31)
  df_thc31 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(31, 0, 500) 
  df_thc31$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(31)
  df_o2sat31 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(31, 0, 500)
  df_o2sat31$O2_sat %>% var(na.rm = TRUE)
  

# 32
  # THC
  df_raw_data %>% f_plot_thc(32)
  df_thc32 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(32, 0, 2000)
  df_thc32$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(32)
  df_o2sat32 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(32, 0, 2000)
  df_o2sat32$O2_sat %>% var(na.rm = TRUE)
  
  
# 33
  # THC
  df_raw_data %>% f_plot_thc(33)
  df_thc33 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(33, 1200, 2000)
  df_thc33$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(33)
  df_o2sat33 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(33, 1200, 2000)
  df_o2sat33$O2_sat %>% var(na.rm = TRUE)
  
  
# 34
  # THC
  df_raw_data %>% f_plot_thc(34)
  df_thc34 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(34, 0, 2000) 
  df_thc34$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(34)
  df_o2sat34 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(34, 0, 2000)
  df_o2sat34$O2_sat %>% var(na.rm = TRUE)
  
  
# 35
  # THC
  df_raw_data %>% f_plot_thc(35)
  df_thc35 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(35, 0, 1000)
  df_thc35$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(35)
  df_o2sat35 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(35, 0, 1000)
  df_o2sat35$O2_sat %>% var(na.rm = TRUE)
  
  
# 36
  # THC
  df_raw_data %>% f_plot_thc(36)
  df_thc36 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(36, 0, 2000)
  df_thc36$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(36)
  df_o2sat36 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(36, 0, 2000)
  df_o2sat36$O2_sat %>% var(na.rm = TRUE)
  
  
# 37
  # THC
  df_raw_data %>% f_plot_thc(37)
  df_thc37 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(37, 0, 2000)
  df_thc37$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(37)
  df_o2sat37 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(37, 0, 2000)
  df_o2sat37$O2_sat %>% var(na.rm = TRUE)
  
  
# 38
  # THC
  df_raw_data %>% f_plot_thc(38)
  df_thc38 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(38, 400, 2000)
  df_thc38$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(38)
  df_o2sat38 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(38, 400, 2000)
  df_o2sat38$O2_sat %>% var(na.rm = TRUE)
  
  
# 39
  # THC
  df_raw_data %>% f_plot_thc(39)
  df_thc39 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(39, 0, 2000)
  df_thc39$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(39)
  df_o2sat39 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(39, 0, 2000)
  df_o2sat39$O2_sat %>% var(na.rm = TRUE)
  
  
# 40
  # THC
  df_raw_data %>% f_plot_thc(40)
  df_thc40 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(40, 0, 2000)
  df_thc40$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(40)
  df_o2sat40 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(40, 0, 2000)
  df_o2sat40$O2_sat %>% var(na.rm = TRUE)
  
  
# 41
  # THC
  df_raw_data %>% f_plot_thc(41)
  df_thc41 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(41, 0, 2000)
  df_thc41$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(41)
  df_o2sat41 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(41, 0, 1000)
  df_o2sat41$O2_sat %>% var(na.rm = TRUE)
  

# 42
  # THC
  df_raw_data %>% f_plot_thc(42)
  df_thc42 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(42, 500, 2000)
  df_thc42$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(42)
  df_o2sat42 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(42, 500, 2000)
  df_o2sat42$O2_sat %>% var(na.rm = TRUE)
  

# 43
  # THC
  df_raw_data %>% f_plot_thc(43)
  df_thc43 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(43, 0, 2000)
  df_thc43$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(43) 
  df_o2sat43 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(43, 200, 650)
  df_o2sat43$O2_sat %>% mean(na.rm = TRUE)

    
# 44
  # THC
  df_raw_data %>% f_plot_thc(44)
  df_thc44 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(44, 300, 2000)
  df_thc44$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_o2sat44 %>% f_plot_o2sat(44) ### signal too noisy

  
# 45
  # THC
  df_raw_data %>% f_plot_thc(45) ### signal too noisy?
  df_thc45 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(45, 400, 1000)
  df_thc45$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(45) ### also very noisy
  df_o2sat45 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(45, 0, 300)
  df_o2sat45$O2_sat %>% var(na.rm = TRUE)
  
  
# 46
  # THC
  df_raw_data %>% f_plot_thc(46)
  df_thc46 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(46, 500, 2000)
  df_thc46$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(46) ### signal very noisy
  df_o2sat46 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(46, 500, 2000)
  df_o2sat46$O2_sat %>% var(na.rm = TRUE)

    
# 47
  # THC
  df_raw_data %>% f_plot_thc(47)
  df_thc47 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(47, 1000, 2000) 
  df_thc47$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(47) ### signal too noisy
  
  
# 48
  # THC
  df_raw_data %>% f_plot_thc(48)
  df_thc48 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(48, 0, 2000)
  df_thc48$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(48)
  df_o2sat48 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(48, 150, 1150)
  df_o2sat48$O2_sat %>% var(na.rm = TRUE)
  
  
# 49
  # THC
  df_thc49 %>% f_plot_thc(49)
  df_thc49 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(49, 150, 2000)
  df_thc49$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(49)
  df_o2sat49 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(49, 150, 2000)
  df_o2sat49$O2_sat %>% var(na.rm = TRUE)
  
  
# 50
  # THC
  df_raw_data %>% f_plot_thc(50)
  df_thc50 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(50, 0, 2000)
  df_thc50$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(50)
  df_o2sat50 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(50, 0, 1000) ### negative mean??
  df_o2sat50$O2_sat %>% mean(na.rm = TRUE)
  
  
# 51
  # THC
  df_raw_data %>% f_plot_thc(51)
  df_thc51 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(51, 0, 2000)
  df_thc51$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(51)
  df_o2sat51 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(51, 0, 2000)
  df_o2sat51$O2_sat %>% var(na.rm = TRUE)
  
  
# 52
  # THC
  df_raw_data %>% f_plot_thc(52)
  df_thc52 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(52, 0, 2000)
  df_thc52$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(52)
  df_o2sat52 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(52, 0, 2000)
  df_o2sat52$O2_sat %>% var(na.rm = TRUE)
  
  
# 53
  # THC
  df_raw_data %>% f_plot_thc(53)
  df_thc53 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(53, 300, 1600)
  df_thc53$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(53)
  df_o2sat53 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(53, 300, 1600)
  df_o2sat53$O2_sat %>% var(na.rm = TRUE)
  

# 54
  # THC
  df_raw_data %>% f_plot_thc(54)
  df_thc54 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(54, 700, 2000)
  df_thc54$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(54)
  df_o2sat54 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(54, 700, 2000)
  df_o2sat54$O2_sat %>% var(na.rm = TRUE)
  
  
# 55
  # THC
  df_raw_data %>% f_plot_thc(55)
  df_thc55 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(55, 0, 2000) %>% dplyr::filter(THC > 0)
  df_thc55$THC %>% var(na.rm = TRUE)
   
  # O2 sat
  df_raw_data %>% f_plot_o2sat(55)
  df_o2sat55 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(55, 0, 2000)
  df_o2sat55$O2_sat %>% var(na.rm = TRUE)
  
  
# 56
  # THC
  df_raw_data %>% f_plot_thc(56)
  df_thc56 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(56, 0, 2000) 
  df_thc56$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(56)
  df_o2sat56 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(56, 0, 2000)
  df_o2sat56$O2_sat %>% var(na.rm = TRUE)
  
  
# 57
  # THC
  df_raw_data %>% f_plot_thc(57)
  df_thc57 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(57, 0, 1000)
  df_thc57$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(57)
  df_o2sat57 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(57, 0, 2000)
  df_o2sat57$O2_sat %>% var(na.rm = TRUE)
  
  
# 58
  # THC
  df_raw_data %>% f_plot_thc(58)
  df_thc58 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(58, 0, 1000)
  df_thc58$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(58)
  df_o2sat58 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(58, 0, 1000)
  df_o2sat58$O2_sat %>% var(na.rm = TRUE)
  
  
# 59
  # THC
  df_raw_data %>% f_plot_thc(59)
  df_thc59 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(59, 1000, 2000)
  df_thc59$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(59)
  df_o2sat59 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(59, 1500, 2000)
  df_o2sat59$O2_sat %>% var(na.rm = TRUE)
  
  
# 60
  # THC
  df_raw_data %>% f_plot_thc(60)
  df_thc60 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(60, 0, 2000)
  df_thc60$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(60)
  df_o2sat60 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(60, 0, 2000)
  df_o2sat60$O2_sat %>% var(na.rm = TRUE)
  
  
# 61
  # THC
  df_raw_data %>% f_plot_thc(61)
  df_thc61 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(61, 0, 1200)
  df_thc61$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(61)
  df_o2sat61 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(61, 0, 1700)
  df_o2sat61$O2_sat %>% var(na.rm = TRUE)
  
  
# 62
  # THC
  df_raw_data %>% f_plot_thc(62)
  df_thc62 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(62, 0, 2000)
  df_thc62$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(62)
  df_o2sat62 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(62, 0, 2000)
  df_o2sat62$O2_sat %>% var(na.rm = TRUE)
  
  
# 63
  # THC
  df_raw_data %>% f_plot_thc(63)
  df_thc63 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(63, 0, 525)
  df_thc63$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(63)
  df_o2sat63 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(63, 0, 525)
  df_o2sat63$O2_sat %>% var(na.rm = TRUE)

    
# 64
  # THC
  df_thc64 %>% f_plot_thc(64) ### signal too noisy

  # O2 sat
  df_raw_data %>% f_plot_o2sat(64) ### signal too noisy

# 65
  # THC
  df_raw_data %>% f_plot_thc(65) ### all negative??
  # O2 sat
  df_raw_data %>% f_plot_o2sat(65) ### signal too noisy
  
# 66
  # THC
  df_raw_data %>% f_plot_thc(66)
  df_thc66 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(66, 0, 1800) 
  df_thc66$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(66)
  df_o2sat66 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(66, 0, 1800)
  df_o2sat66$O2_sat %>% var(na.rm = TRUE)
  
  
# 67
  # THC
  df_raw_data %>% f_plot_thc(67)
  df_thc67 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(67, 0, 2000) 
  df_thc67$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(67)
  df_o2sat67 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(67, 0, 2000)
  df_o2sat67$O2_sat %>% var(na.rm = TRUE)
  

# 68
  # THC
  df_raw_data %>% f_plot_thc(68)
  df_thc68 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(68, 0, 1400) 
  df_thc68$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(68)
  df_o2sat68 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(68, 0, 1400)
  df_o2sat68$O2_sat %>% var(na.rm = TRUE)
  
  
# 69
  # THC
  df_raw_data %>% f_plot_thc(69)
  df_thc69 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(69, 0, 2000) 
  df_thc69$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(69)
  df_o2sat69 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(69, 0, 2000)
  df_o2sat69$O2_sat %>% var(na.rm = TRUE)
  
  
# 70
  # THC
  df_raw_data %>% f_plot_thc(70)
  df_thc70 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(70, 0, 2000)
  df_thc70$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(70)
  df_o2sat70 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(70, 0, 2000)
  df_o2sat70$O2_sat %>% var(na.rm = TRUE)
  

# 71
  # THC
  df_raw_data %>% f_plot_thc(71)
  df_thc71 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(71, 750, 1550)
  df_thc71$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(71)
  df_o2sat71 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(71, 790, 1500)
  df_o2sat71$O2_sat %>% var(na.rm = TRUE)
  
  
# 72
  # THC
  df_raw_data %>% f_plot_thc(72)
  df_thc72 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(72, 0, 2000) 
  df_thc72$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(72)
  df_o2sat72 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(72, 0, 2000) %>% dplyr::filter(O2_sat < 500)
  df_o2sat72$O2_sat %>% var(na.rm = TRUE)
   
  
# 73
  # THC
  df_raw_data %>% f_plot_thc(73)
  df_thc73 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(73, 0, 2000)
  df_thc73$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(73)
  df_o2sat73 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(73, 0, 2000)
  df_o2sat73$O2_sat %>% var(na.rm = TRUE)
  
  
# 74
  # THC
  df_raw_data %>% f_plot_thc(74)
  df_thc74 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(74, 0, 2000)
  df_thc74$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(74)
  df_o2sat74 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(74, 0, 2000)
  df_o2sat74$O2_sat %>% var(na.rm = TRUE)
  
  
# 75
  # THC
  df_raw_data %>% f_plot_thc(75)
  df_thc75 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(75, 0, 2000)
  df_thc75$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(75)
  df_o2sat75 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(75, 0, 2000)
  df_o2sat75$O2_sat %>% var(na.rm = TRUE)
  
  
# 76
  # THC
  df_raw_data %>% f_plot_thc(76)
  df_thc76 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(76, 0, 2000)
  df_thc76$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(76)
  df_o2sat76 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(76, 0, 2000)
  df_o2sat76$O2_sat %>% var(na.rm = TRUE)
  
  
# 77
  # THC
  df_raw_data %>% f_plot_thc(77)
  df_thc77 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(77, 0, 2000)
  df_thc77$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(77)
  df_o2sat77 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(77, 0, 2000)
  df_o2sat77$O2_sat %>% var(na.rm = TRUE)
  
  
# 78
  # THC
  df_raw_data %>% f_plot_thc(78)
  df_thc78 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(78, 0, 2000)
  df_thc78$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(78)
  df_o2sat78 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(78, 0, 2000)
  df_o2sat78$O2_sat %>% var(na.rm = TRUE)
  
  
# 79
  # THC
  df_raw_data %>% f_plot_thc(79)
  df_thc79 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(79, 0, 2000)
  df_thc79$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(79)
  df_o2sat79 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(79, 0, 2000)
  df_o2sat79$O2_sat %>% var(na.rm = TRUE)
  
  
# 80
  # THC
  df_raw_data %>% f_plot_thc(80)
  df_thc80 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(80, 0, 2000)
  df_thc80$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(80)
  df_o2sat80 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(80, 0, 2000)
  df_o2sat80$O2_sat %>% var(na.rm = TRUE)
  
  
# 81
  # THC
  df_raw_data %>% f_plot_thc(81)
  df_thc81 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(81, 0, 2000) %>% dplyr::filter(THC > 0)
  df_thc81$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(81)
  df_o2sat81 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(81, 0, 2000)
  df_o2sat81$O2_sat %>% var(na.rm = TRUE)
  
  
# 82
  # THC
  df_raw_data %>% f_plot_thc(82)
  df_thc82 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(82, 0, 2000) 
  df_thc82$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(82)
  df_o2sat82 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(82, 0, 2000)
  df_o2sat82$O2_sat %>% var(na.rm = TRUE)
  
  
# 83
  # THC
  df_raw_data %>% f_plot_thc(83)
  df_thc83 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(83, 0, 2000)
  df_thc83$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(83)
  df_o2sat83 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(83, 0, 2000)
  df_o2sat83$O2_sat %>% var(na.rm = TRUE)
  
  
# 84
  # THC
  df_raw_data %>% f_plot_thc(84)
  df_thc84 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(84, 0, 2000)
  df_thc84$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(84)
  df_o2sat84 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(84, 0, 2000)
  df_o2sat84$O2_sat %>% var(na.rm = TRUE)
  
  
# 85
  # THC
  df_raw_data %>% f_plot_thc(85)
  df_thc85 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(85, 0, 2000)
  df_thc85$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(85)
  df_o2sat85 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(85, 0, 2000)
  df_o2sat85$O2_sat %>% var(na.rm = TRUE)
  
  
# 86
  # THC
  df_raw_data %>% f_plot_thc(86)
  df_thc86 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(86, 0, 2000)
  df_thc86$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(86)
  df_o2sat86 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(86, 0, 2000)
  df_o2sat86$O2_sat %>% var(na.rm = TRUE)
  
  
# 87
  # THC
  df_raw_data %>% f_plot_thc(87)
  df_thc87 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(87, 0, 2000)
  df_thc87$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(87)
  df_o2sat87 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(87, 0, 2000)
  df_o2sat87$O2_sat %>% var(na.rm = TRUE)
  
  
# 88
  # THC
  df_raw_data %>% f_plot_thc(88)
  df_thc88 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(88, 0, 2000)
  df_thc88$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(88)
  df_o2sat88 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(88, 0, 2000)
  df_o2sat88$O2_sat %>% var(na.rm = TRUE)
  

# 89
  # THC
  df_raw_data %>% f_plot_thc(89)
  df_thc89 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(89, 0, 2000)
  df_thc89$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(89)
  df_o2sat89 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(89, 0, 2000)
  df_o2sat89$O2_sat %>% var(na.rm = TRUE)
  
  
# 90
  # THC
  df_raw_data %>% f_plot_thc(90)
  df_thc90 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(90, 0, 2000) %>% dplyr::filter(THC > 0 & THC < 400)
  df_thc90$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(90)
  df_o2sat90 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(90, 0, 2000) %>% dplyr::filter(O2_sat > -600)
  df_o2sat90$O2_sat %>% var(na.rm = TRUE)
  
  
# 91
  # THC
  df_raw_data %>% f_plot_thc(91)
  df_thc91 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(91, 0, 800)
  df_thc91$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(91)
  df_o2sat91 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(91, 0, 800)
  df_o2sat91$O2_sat %>% var(na.rm = TRUE)
  
  
# 92
  # THC
  df_raw_data %>% f_plot_thc(92)
  df_thc92 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(92, 0, 2000)
  df_thc92$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(92)
  df_o2sat92 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(92, 0, 2000)
  df_o2sat92$O2_sat %>% var(na.rm = TRUE)
  

# 93
  # THC
  df_raw_data %>% f_plot_thc(93)
  df_thc93 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(93, 0, 2000)
  df_thc93$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(93)
  df_o2sat93 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(93, 0, 2000)
  df_o2sat93$O2_sat %>% var(na.rm = TRUE)
  
  
# 94
  # THC
  df_raw_data %>% f_plot_thc(94)
  df_thc94 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(94, 0, 2000)
  df_thc94$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(94)
  df_o2sat94 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(94, 0, 2000)
  df_o2sat94$O2_sat %>% var(na.rm = TRUE)
  
  
# 95
  # THC
  df_raw_data %>% f_plot_thc(95)
  df_thc95 <- df_raw_data %>% dplyr::select(-O2_sat) %>% f_segment_signal(95, 0, 2000)
  df_thc95$THC %>% var(na.rm = TRUE)
  
  # O2 sat
  df_raw_data %>% f_plot_o2sat(95)
  df_o2sat95 <- df_raw_data %>% dplyr::select(-THC) %>% f_segment_signal(95, 0, 2000)
  df_o2sat95$O2_sat %>% var(na.rm = TRUE)
  


# combine all  
  df_thc_filt <- bind_rows(df_thc1, df_thc2, df_thc3, df_thc6, df_thc7, df_thc8, df_thc9, df_thc10, #4, 5,  and 6 are all the exact same
                           df_thc11, df_thc12, df_thc13, df_thc14, df_thc15, df_thc16, df_thc17, df_thc18, df_thc19, df_thc20,
                           df_thc22, df_thc23, df_thc24, df_thc27, df_thc28, df_thc29, df_thc30,
                           df_thc31, df_thc32, df_thc33, df_thc34, df_thc35, df_thc36, df_thc37, df_thc38, df_thc39, df_thc40,
                           df_thc41, df_thc42, df_thc43, df_thc44, df_thc45, df_thc46, df_thc47, df_thc48, df_thc49, df_thc50,
                           df_thc51, df_thc52, df_thc53, df_thc54, df_thc55, df_thc56, df_thc57, df_thc58, df_thc59, df_thc60,
                           df_thc61, df_thc62, df_thc63, df_thc66, df_thc67, df_thc68, df_thc69, df_thc70,
                           df_thc71, df_thc72, df_thc73, df_thc74, df_thc75, df_thc76, df_thc77, df_thc78, df_thc79, df_thc80,
                           df_thc81, df_thc82, df_thc83, df_thc84, df_thc85, df_thc86, df_thc87, df_thc88, df_thc89, df_thc90,
                           df_thc91, df_thc92, df_thc93, df_thc94, df_thc95)
  

  df_o2sat_filt <- bind_rows(df_o2sat1, df_o2sat2, df_o2sat3, df_o2sat6, df_o2sat7, df_o2sat8, df_o2sat9, df_o2sat10, #4, 5,  and 6 are all the exact same
                             df_o2sat11, df_o2sat12, df_o2sat13, df_o2sat15, df_o2sat16, df_o2sat17, df_o2sat18, df_o2sat19, df_o2sat20,
                             df_o2sat22, df_o2sat23, df_o2sat24, df_o2sat27, df_o2sat28, df_o2sat29, df_o2sat30,
                             df_o2sat31, df_o2sat32, df_o2sat33, df_o2sat34, df_o2sat35, df_o2sat36, df_o2sat37, df_o2sat38, df_o2sat39, df_o2sat40,
                             df_o2sat41, df_o2sat42, df_o2sat43, df_o2sat45, df_o2sat46, df_o2sat48, df_o2sat49, df_o2sat50,
                             df_o2sat51, df_o2sat52, df_o2sat53, df_o2sat54, df_o2sat55, df_o2sat56, df_o2sat57, df_o2sat58, df_o2sat59, df_o2sat60,
                             df_o2sat61, df_o2sat62, df_o2sat63, df_o2sat66, df_o2sat67, df_o2sat68, df_o2sat69, df_o2sat70,
                             df_o2sat71, df_o2sat72, df_o2sat73, df_o2sat74, df_o2sat75, df_o2sat76, df_o2sat77, df_o2sat78, df_o2sat79, df_o2sat80,
                             df_o2sat81, df_o2sat82, df_o2sat83, df_o2sat84, df_o2sat85, df_o2sat86, df_o2sat87, df_o2sat88, df_o2sat89, df_o2sat90,
                             df_o2sat91, df_o2sat92, df_o2sat93, df_o2sat94, df_o2sat95)
  
  
# calculate average Hb_tot and Hb_oxy sat ------------------------------------------------

  # total THC
    df_THCmean <- df_thc_filt %>% 
      dplyr::filter(!is.na(THC)) %>% 
      group_by(number, Status) %>% 
      summarise(THC_mean = mean(THC))
    
    # plot
    df_THCmean  %>%  
      ggplot(aes(x = Status, y = THC_mean)) +
      geom_violin(aes(color = Status)) +
      geom_jitter(position = position_jitter(0.2), shape = 1) +
      stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
      
      theme_bw() +
      labs(y = "uM") +
      ggtitle("Cerebral tissue hemoglobin concentration") +
      theme(legend.position = "none") 
    
  
  # O2 sat
    df_o2sat_mean <- df_o2sat_filt %>% 
      dplyr::filter(!is.na(O2_sat)) %>% 
      group_by(number, Status) %>% 
      summarise(o2sat_mean = mean(O2_sat))
    
    # plot
    df_o2sat_mean %>%  
      ggplot(aes(x = Status, y = o2sat_mean)) +
      geom_violin(aes(color = Status)) +
      geom_jitter(position = position_jitter(0.2), shape = 1) +
      stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
      
      theme_bw() +
      labs(y = "percent") +
      ggtitle("Cerebral tissue hemoglobin oxygen saturation") +
      theme(legend.position = "none") ## negative values??
    
  
    # explore negative O2 sat
      df_o2sat_mean %>% dplyr::filter(o2sat_mean < 0)
      df_raw_data %>% dplyr::filter(number %in% c(3, 50)) %>% count(number, subject_id)
      # remove number 3, that is a third visit (#2 is a follow-up!!)
      # look at number 50 (TM1003UM01)
      # read in
      df_TM1003UM01 <- read.delim2(file = my_files[50], header = FALSE, sep = "\t") %>% as_tibble()
      # format
      df_TM1003UM01 <- df_TM1003UM01 %>% 
        dplyr::filter(!grepl("Patient|Data|[R|r]aw", V1)) %>% 
        dplyr::select(c(2:6)) %>% 
        dplyr::rename("Time" = "V2", "O2_sat" = "V3", "THC" = "V4", "OxyHb" = "V5", "DeoxyHb" = "V6") %>% 
        mutate_if(is.character, as.numeric)
      df_TM1003UM01 %>% summarise(avg_oxy = mean(OxyHb, na.rm = TRUE), avg_deoxy = mean(DeoxyHb, na.rm = TRUE)) # negative OxyHb value
      # remove negative values
      df_o2sat_mean %>% dplyr::filter(!(number %in% c(2, 3, 50))) %>% 
        ggplot(aes(x = Status, y = o2sat_mean)) +
        geom_violin(aes(color = Status)) +
        geom_jitter(position = position_jitter(0.2), shape = 1) +
        stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
        
        theme_bw() +
        labs(y = "percent") +
        ggtitle("Cerebral tissue hemoglobin oxygen saturation") +
        theme(legend.position = "none") ## negative values??
      

      

# export to Matlab for signal processing ----------------------------------

      
    # add milliseconds in to Time
      
      # THC
      df_thc_filt <- df_thc_filt %>% 
        group_by(number, Time) %>% 
        left_join(count(.)) %>% 
        mutate(Time = ifelse(n == 49, Time + 0.02*row_number(), Time + 0.02*(row_number() - 1))) %>% 
        dplyr::select(-n)
      
      # O2 sat
      df_o2sat_filt <- df_o2sat_filt %>% 
        group_by(number, Time) %>% 
        left_join(count(.)) %>% 
        mutate(Time = ifelse(n == 49, Time + 0.02*row_number(), Time + 0.02*(row_number() - 1))) %>% 
        dplyr::select(-n)
      
    # export to .txt to filter in Matlab
      
      # Hb_tot
      l_thc_filt <- df_thc_filt %>% split(df_thc_filt$number)
      thc_names <- df_thc_filt$subject_id %>% unique()
      names(l_thc_filt) <- thc_names
      path <- "Data/signal_segments/Hb_tot/"
      1:length(l_thc_filt) %>% map(~write.table(l_thc_filt[[.x]], file = paste0(path, names(l_thc_filt[.x]), ".txt")))
      
      # Hb_oxy
      l_o2sat_filt <- df_o2sat_filt %>% split(df_o2sat_filt$number)
      o2sat_names <- df_o2sat_filt$subject_id %>% unique()
      names(l_o2sat_filt) <- o2sat_names
      path <- "Data/signal_segments/Hb_oxy/"
      1:length(l_o2sat_filt) %>% map(~write.table(l_o2sat_filt[[.x]], file = paste0(path, names(l_o2sat_filt[.x]), ".txt")))
      
      
      
   ### moving average & signal filtering are done in MatLab, read back into R for DFA   

      
      
      
      
      
