


# read in data -------------------------------------------------------------------------

  library(tidyverse)
  library(purrr)
  

  # read in
  my_files <- list.files(path = "Data/Raw data/", pattern = "*.txt", full.names = TRUE, recursive = FALSE)
  
  l_raw_data <- 1:length(my_files) %>% 
    purrr::map(~read.delim2(file = my_files[.x], header = FALSE, sep = "\t") %>% as_tibble()) 
  
  # name with Status
  l_names <- 1:length(my_files) %>% 
    purrr::map(~substr(my_files[.x], start = 22, stop = 23)) 
  
  names(l_raw_data) <- l_names


  
# format data -------------------------------------------------------------


  library(data.table)
  library(zoo)
  
  # format to just THC in workable form
  l_thc <- 1:length(l_raw_data) %>% 
    
    purrr::map(~mutate(l_raw_data[[.x]], Status = l_names[[.x]])) %>% 
    purrr::map(~filter(.x, !grepl("Patient|Data|[R|r]aw", V1))) %>% 
    purrr::map(~dplyr::select(.x, c(ncol(.x), 2,4))) %>% 
    purrr::map(~dplyr::rename(.x, "Time" = "V2", "THC" = "V4")) %>% 
    purrr::map(~mutate(.x, Time = as.numeric(Time))) %>% 
    purrr::map(~mutate(.x, THC = as.numeric(THC)))
  
  # renumber patients
  l_thc <- 1:length(l_thc) %>% 
    purrr::map(~mutate(l_thc[[.x]], number = c(1:length(l_thc))[.x])) %>% 
    purrr::map(~dplyr::select(.x, c(number, everything())))
  
  # condense list to df
  df_thc <- rbindlist(l_thc) %>% 
    as_tibble() 
  
  # make sure status is HV, UM, or CM
  df_thc$Status[df_thc$Status == "HC"] <- "HV"
  df_thc <- df_thc %>% dplyr::filter(Status != "-0") # unclear what group this patient was in
  
  # smooth data with 1s moving average
  df_thc <- df_thc %>% 
    group_by(number) %>% 
    mutate(roll_mean = zoo::rollmean(THC, k = 50, fill = NA)) %>% 
    ungroup()
  




# explore --------------------------------------------------------------------


  # distribution of data by status only
  df_thc %>% 
    
    ggplot(aes(roll_mean)) +
    geom_histogram() +
    facet_wrap(~Status, scales = "free")
  
  # distribution by patient number
  df_thc %>% 
    
    ggplot(aes(roll_mean)) +
    geom_histogram() +
    facet_grid(~number, scales = "free") 

  # plot all signals together
  df_thc %>% 
    
    ggplot(aes(x = Time, y = roll_mean)) +
    geom_line() +
    facet_wrap(~Status, scales = "free_y") +
    
    theme_bw()
  
  # facet by number
  df_thc %>% 
    
    ggplot(aes(x = Time, y = roll_mean)) +
    geom_line() +
    facet_wrap(~number, scales = "free_y") +
    
    theme_bw()
  
  


# remove blips ------------------------------------------------------------


# remove signals that aren't at least 200s  
  
  df_thc %>% group_by(number) %>% count() %>% filter(n < 20000) # 24 and 30 not long enough
  df_thc <- df_thc %>% filter(!(number %in% c(24, 30)))

# plotting function
  
  f_plot_pt <- function(df, num) {
    
    df %>% 
      filter(number == num) %>% 
      ggplot(aes(x = Time, y = roll_mean)) +
      geom_line()
    
  }
  
# segment signal function
  
  f_segment_signal <- function(df, num, start, end) {
    
    df %>% 
      filter(number == num) %>% 
      filter(Time > start & Time < end)
  }

# 1
  df_thc %>% f_plot_pt(1) # eyeball where after blips is
  df_thc1 <- df_thc %>% f_segment_signal(1, 1200, 2500) # shorten signal to that length, must be >200s
  df_thc1$roll_mean %>% var(na.rm = TRUE) # check var less than 1000
  
# 2
  df_thc %>% f_plot_pt(2) 
  df_thc2 <- df_thc %>% f_segment_signal(2, 900, 2000)
  df_thc2$roll_mean %>% var(na.rm = TRUE)
  
# 3
  df_thc %>% f_plot_pt(3) 
  df_thc3 <- df_thc %>% f_segment_signal(3, 500, 2000)
  df_thc3$roll_mean %>% var(na.rm = TRUE)

# 4
  df_thc %>% f_plot_pt(4) 
  df_thc4 <- df_thc %>% f_segment_signal(4, 0, 2000) # no blips!
  df_thc4$roll_mean %>% var(na.rm = TRUE)
  
# 5
  df_thc %>% f_plot_pt(5) # ***same as 4???
  df_thc5 <- df_thc %>% f_segment_signal(5, 0, 2000) # no blips!
  df_thc5$roll_mean %>% var(na.rm = TRUE)
  
# 6
  df_thc %>% f_plot_pt(6) # same as 4 & 5
  df_thc6 <- df_thc %>% f_segment_signal(6, 0, 2000) # no blips!
  df_thc6$roll_mean %>% var(na.rm = TRUE)
  
# 7
  df_thc %>% f_plot_pt(7)
  df_thc7 <- df_thc %>% f_segment_signal(7, 0, 2000) # no blips!
  df_thc7$roll_mean %>% var(na.rm = TRUE)
  
# 8
  df_thc %>% f_plot_pt(8)
  df_thc8 <- df_thc %>% f_segment_signal(8, 0, 2000) # no blips!
  df_thc8$roll_mean %>% var(na.rm = TRUE)
  
# 9
  df_thc %>% f_plot_pt(9)
  df_thc9 <- df_thc %>% f_segment_signal(9, 0, 2000) # no blips!
  df_thc9$roll_mean %>% var(na.rm = TRUE)
  
# 10
  df_thc %>% f_plot_pt(10)
  df_thc10 <- df_thc %>% f_segment_signal(10, 0, 2000) # no blips!
  df_thc10$roll_mean %>% var(na.rm = TRUE)
  
# 11
  df_thc %>% f_plot_pt(11)
  df_thc11 <- df_thc %>% f_segment_signal(11, 0, 2000) # no blips!
  df_thc11$roll_mean %>% var(na.rm = TRUE)
  
# 12
  df_thc %>% f_plot_pt(12)
  df_thc12 <- df_thc %>% f_segment_signal(12, 0, 2000) # no blips!
  df_thc12$roll_mean %>% var(na.rm = TRUE)

# 13
  df_thc %>% f_plot_pt(13)
  df_thc13 <- df_thc %>% f_segment_signal(13, 0, 2000) # no blips?
  df_thc13$roll_mean %>% var(na.rm = TRUE) # technically within range
  
# 14
  df_thc %>% f_plot_pt(14)
  df_thc14 <- df_thc %>% f_segment_signal(14, 400, 2000)
  df_thc14$roll_mean %>% var(na.rm = TRUE)
  
# 15
  df_thc %>% f_plot_pt(15)
  df_thc15 <- df_thc %>% f_segment_signal(15, 0, 2000) # no blips!
  df_thc15$roll_mean %>% var(na.rm = TRUE)

# 16
  df_thc %>% f_plot_pt(16)
  df_thc16 <- df_thc %>% f_segment_signal(16, 0, 2000) # no blips!
  df_thc16$roll_mean %>% var(na.rm = TRUE)
  
# 17
  df_thc %>% f_plot_pt(17)
  df_thc17 <- df_thc %>% f_segment_signal(17, 0, 2000) # no blips!
  df_thc17$roll_mean %>% var(na.rm = TRUE)

# 18
  df_thc %>% f_plot_pt(18)
  df_thc18 <- df_thc %>% f_segment_signal(18, 0, 2000) # no blips!
  df_thc18$roll_mean %>% var(na.rm = TRUE)

# 19
  df_thc %>% f_plot_pt(19)
  df_thc19 <- df_thc %>% f_segment_signal(19, 0, 2000) # no blips!
  df_thc19$roll_mean %>% var(na.rm = TRUE)

# 20
  df_thc %>% f_plot_pt(20)
  df_thc20 <- df_thc %>% f_segment_signal(20, 0, 2000) # no blips!
  df_thc20$roll_mean %>% var(na.rm = TRUE)
  
# 21
  df_thc %>% f_plot_pt(21)
  df_thc21 <- df_thc %>% f_segment_signal(21, 0, 2000) # no blips!
  df_thc21$roll_mean %>% var(na.rm = TRUE)
  
# 22
  df_thc %>% f_plot_pt(22)
  df_thc22 <- df_thc %>% f_segment_signal(22, 0, 2000) # no blips!
  df_thc22$roll_mean %>% var(na.rm = TRUE)
  
# 23
  df_thc23 %>% f_plot_pt(23)
  df_thc23 <- df_thc %>% f_segment_signal(23, 1350, 2000)
  df_thc23$roll_mean %>% var(na.rm = TRUE)
  
# 25
  df_thc %>% f_plot_pt(25)
  df_thc25 <- df_thc %>% f_segment_signal(25, 0, 2000) # no blips!
  df_thc25$roll_mean %>% var(na.rm = TRUE)
  
# 26
  df_thc %>% f_plot_pt(26)
  df_thc26 <- df_thc %>% f_segment_signal(26, 0, 2000) # no blips!
  df_thc26$roll_mean %>% var(na.rm = TRUE)

# 27
  df_thc %>% f_plot_pt(27)
  df_thc27 <- df_thc %>% f_segment_signal(27, 0, 2000) # no blips!
  df_thc27$roll_mean %>% var(na.rm = TRUE)
  
# 28
  df_thc %>% f_plot_pt(28)
  df_thc28 <- df_thc %>% f_segment_signal(28, 0, 2000) # no blips!
  df_thc28$roll_mean %>% var(na.rm = TRUE)
  
# 31
  df_thc %>% f_plot_pt(31)
  df_thc31 <- df_thc %>% f_segment_signal(31, 0, 2000) # no blips!
  df_thc31$roll_mean %>% var(na.rm = TRUE)
  
# 32
  df_thc %>% f_plot_pt(32)
  df_thc32 <- df_thc %>% f_segment_signal(32, 0, 2000) # no blips!
  df_thc32$roll_mean %>% var(na.rm = TRUE)
  
# 33
  df_thc %>% f_plot_pt(33)
  df_thc33 <- df_thc %>% f_segment_signal(33, 0, 2000) # no blips!
  df_thc33$roll_mean %>% var(na.rm = TRUE)
  
# 34
  df_thc %>% f_plot_pt(34)
  df_thc34 <- df_thc %>% f_segment_signal(34, 0, 430)
  df_thc34$roll_mean %>% var(na.rm = TRUE)
  
# 35
  df_thc %>% f_plot_pt(35)
  df_thc35 <- df_thc %>% f_segment_signal(35, 0, 510)
  df_thc35$roll_mean %>% var(na.rm = TRUE)
  
# 36
  df_thc %>% f_plot_pt(36)
  df_thc36 <- df_thc %>% f_segment_signal(36, 0, 2000) # no blips!
  df_thc36$roll_mean %>% var(na.rm = TRUE)

# 37
  df_thc %>% f_plot_pt(37)
  df_thc37 <- df_thc %>% f_segment_signal(37, 700, 2000)
  df_thc37$roll_mean %>% var(na.rm = TRUE)
  
# 38
  df_thc %>% f_plot_pt(38)
  df_thc38 <- df_thc %>% f_segment_signal(38, 0, 2000) # no blips!
  df_thc38$roll_mean %>% var(na.rm = TRUE)
  
# 39
  df_thc %>% f_plot_pt(39)
  df_thc39 <- df_thc %>% f_segment_signal(39, 0, 1100) 
  df_thc39$roll_mean %>% var(na.rm = TRUE)
  
# 40
  df_thc %>% f_plot_pt(40)
  df_thc40 <- df_thc %>% f_segment_signal(40, 0, 2000) # no blips!
  df_thc40$roll_mean %>% var(na.rm = TRUE)
  
# 41
  df_thc %>% f_plot_pt(41)
  df_thc41 <- df_thc %>% f_segment_signal(41, 0, 2000) # no blips!
  df_thc41$roll_mean %>% var(na.rm = TRUE)
  
# 42
  df_thc %>% f_plot_pt(42)
  df_thc42 <- df_thc %>% f_segment_signal(42, 400, 2000)
  df_thc42$roll_mean %>% var(na.rm = TRUE)
  
# 43
  df_thc %>% f_plot_pt(43)
  df_thc43 <- df_thc %>% f_segment_signal(43, 0, 2000) # no blips!
  df_thc43$roll_mean %>% var(na.rm = TRUE)
  
# 44
  df_thc %>% f_plot_pt(44)
  df_thc44 <- df_thc %>% f_segment_signal(44, 0, 2000) # no blips!
  df_thc44$roll_mean %>% var(na.rm = TRUE)
  
# 45
  df_thc %>% f_plot_pt(45)
  df_thc45 <- df_thc %>% f_segment_signal(45, 0, 2000) # no blips!
  df_thc45$roll_mean %>% var(na.rm = TRUE)
  
# 46
  df_thc %>% f_plot_pt(46)
  df_thc46 <- df_thc %>% f_segment_signal(46, 0, 2000) # no blips!
  df_thc46$roll_mean %>% var(na.rm = TRUE)

# 47
  df_thc %>% f_plot_pt(47)
  df_thc47 <- df_thc %>% f_segment_signal(47, 0, 2000) # no blips!
  df_thc47$roll_mean %>% var(na.rm = TRUE)

# 48
  df_thc %>% f_plot_pt(48)
  df_thc48 <- df_thc %>% f_segment_signal(48, 200, 2000)
  df_thc48$roll_mean %>% var(na.rm = TRUE)
  
# 49
  df_thc %>% f_plot_pt(49)
  df_thc49 <- df_thc %>% f_segment_signal(49, 350, 2000)
  df_thc49$roll_mean %>% var(na.rm = TRUE)
  
# 50
  df_thc %>% f_plot_pt(50)
  df_thc50 <- df_thc %>% f_segment_signal(50, 490, 2000)
  df_thc50$roll_mean %>% var(na.rm = TRUE)
  
# 51
  df_thc %>% f_plot_pt(51)
  df_thc51 <- df_thc %>% f_segment_signal(51, 0, 2000) # no blips!
  df_thc51$roll_mean %>% var(na.rm = TRUE)
  
# 52
  df_thc %>% f_plot_pt(52)
  df_thc52 <- df_thc %>% f_segment_signal(52, 0, 2000) # no blips!
  df_thc52$roll_mean %>% var(na.rm = TRUE)
  
# 53
  df_thc %>% f_plot_pt(53)
  df_thc53 <- df_thc %>% f_segment_signal(53, 150, 2000)
  df_thc53$roll_mean %>% var(na.rm = TRUE)
  
# 54
  df_thc %>% f_plot_pt(54)
  df_thc54 <- df_thc %>% f_segment_signal(54, 0, 2000) # no blips!
  df_thc54$roll_mean %>% var(na.rm = TRUE)
  
# 55
  df_thc %>% f_plot_pt(55)
  df_thc55 <- df_thc %>% f_segment_signal(55, 0, 2000) # no blips!
  df_thc55$roll_mean %>% var(na.rm = TRUE)
  
# 56
  df_thc %>% f_plot_pt(56)
  df_thc56 <- df_thc %>% f_segment_signal(56, 0, 2100) # no blips!
  df_thc56$roll_mean %>% var(na.rm = TRUE)
  
# 57
  df_thc %>% f_plot_pt(57)
  df_thc57 <- df_thc %>% f_segment_signal(57, 0, 2000) # no blips!
  df_thc57$roll_mean %>% var(na.rm = TRUE)
  
# 58
  df_thc %>% f_plot_pt(58)
  df_thc58 <- df_thc %>% f_segment_signal(58, 700, 2000)
  df_thc58$roll_mean %>% var(na.rm = TRUE)

# 59
  df_thc %>% f_plot_pt(59)
  df_thc59 <- df_thc %>% f_segment_signal(59, 0, 2000) # no blips!
  df_thc59$roll_mean %>% var(na.rm = TRUE)
  
# 60
  df_thc %>% f_plot_pt(60)
  df_thc60 <- df_thc %>% f_segment_signal(60, 0, 2000) # no blips!
  df_thc60$roll_mean %>% var(na.rm = TRUE)
  
# 61
  df_thc %>% f_plot_pt(61)
  df_thc61 <- df_thc %>% f_segment_signal(61, 0, 2000) # no blips!
  df_thc61$roll_mean %>% var(na.rm = TRUE)
  
# 62
  df_thc %>% f_plot_pt(62)
  df_thc62 <- df_thc %>% f_segment_signal(62, 0, 1150)
  df_thc62$roll_mean %>% var(na.rm = TRUE)
  
# 63
  df_thc %>% f_plot_pt(63)
  df_thc63 <- df_thc %>% f_segment_signal(63, 1350, 2000)
  df_thc63$roll_mean %>% var(na.rm = TRUE)
  
# 64
  df_thc %>% f_plot_pt(64)
  df_thc64 <- df_thc %>% f_segment_signal(64, 0, 2000) # no blips!
  df_thc64$roll_mean %>% var(na.rm = TRUE)
  
# 65
  df_thc %>% f_plot_pt(65)
  df_thc65 <- df_thc %>% f_segment_signal(65, 0, 1700)
  df_thc65$roll_mean %>% var(na.rm = TRUE)
  
# 66
  df_thc %>% f_plot_pt(66)
  df_thc66 <- df_thc %>% f_segment_signal(66, 0, 2000) # no blips!
  df_thc66$roll_mean %>% var(na.rm = TRUE)
  
# 67
  df_thc %>% f_plot_pt(67)
  df_thc67 <- df_thc %>% f_segment_signal(67, 0, 525)
  df_thc67$roll_mean %>% var(na.rm = TRUE)
  
# 68
  df_thc %>% f_plot_pt(68)
  df_thc68 <- df_thc %>% f_segment_signal(68, 1350, 1600)
  df_thc68$roll_mean %>% var(na.rm = TRUE)
  
# 69
  df_thc %>% f_plot_pt(69) # all negative??

# 70
  df_thc %>% f_plot_pt(70)
  df_thc70 <- df_thc %>% f_segment_signal(70, 0, 1900)
  df_thc70$roll_mean %>% var(na.rm = TRUE)
  
# 71
  df_thc %>% f_plot_pt(71)
  df_thc71 <- df_thc %>% f_segment_signal(71, 0, 2000) # no blips!
  df_thc71$roll_mean %>% var(na.rm = TRUE)
  
# 71
  df_thc %>% f_plot_pt(71)
  df_thc71 <- df_thc %>% f_segment_signal(71, 0, 2000) # no blips!
  df_thc71$roll_mean %>% var(na.rm = TRUE)
  
# 72
  df_thc %>% f_plot_pt(72)
  df_thc72 <- df_thc %>% f_segment_signal(72, 0, 1400) 
  df_thc72$roll_mean %>% var(na.rm = TRUE)

# 73
  df_thc %>% f_plot_pt(73)
  df_thc73 <- df_thc %>% f_segment_signal(73, 0, 2000) # no blips!
  df_thc73$roll_mean %>% var(na.rm = TRUE)
  
# 74
  df_thc %>% f_plot_pt(74)
  df_thc74 <- df_thc %>% f_segment_signal(74, 0, 2000) # no blips!
  df_thc74$roll_mean %>% var(na.rm = TRUE)
  
# 75
  df_thc %>% f_plot_pt(75)
  df_thc75 <- df_thc %>% f_segment_signal(75, 0, 1800)
  df_thc75$roll_mean %>% var(na.rm = TRUE)

# 76
  df_thc %>% f_plot_pt(76)
  df_thc76 <- df_thc %>% f_segment_signal(76, 0, 2000) # no blips!
  df_thc76$roll_mean %>% var(na.rm = TRUE)
  
# 77
  df_thc %>% f_plot_pt(77)
  df_thc77 <- df_thc %>% f_segment_signal(77, 0, 2000) # no blips!
  df_thc77$roll_mean %>% var(na.rm = TRUE)
  
# 78
  df_thc %>% f_plot_pt(78)
  df_thc78 <- df_thc %>% f_segment_signal(78, 0, 2000) # no blips!
  df_thc78$roll_mean %>% var(na.rm = TRUE)
  
# 79
  df_thc %>% f_plot_pt(79)
  df_thc79 <- df_thc %>% f_segment_signal(79, 0, 2000) # no blips!
  df_thc79$roll_mean %>% var(na.rm = TRUE)
  
# 80
  df_thc %>% f_plot_pt(80)
  df_thc80 <- df_thc %>% f_segment_signal(80, 0, 2000) # no blips!
  df_thc80$roll_mean %>% var(na.rm = TRUE)
  
# 81
  df_thc %>% f_plot_pt(81)
  df_thc81 <- df_thc %>% f_segment_signal(81, 0, 2000) # no blips!
  df_thc81$roll_mean %>% var(na.rm = TRUE)
  
# 82
  df_thc %>% f_plot_pt(82)
  df_thc82 <- df_thc %>% f_segment_signal(82, 0, 2000) # no blips!
  df_thc82$roll_mean %>% var(na.rm = TRUE)
  
# 83
  df_thc %>% f_plot_pt(83)
  df_thc83 <- df_thc %>% f_segment_signal(83, 0, 2000) # no blips!
  df_thc83$roll_mean %>% var(na.rm = TRUE)
  
# 84
  df_thc %>% f_plot_pt(84)
  df_thc84 <- df_thc %>% f_segment_signal(84, 0, 2000) # no blips!
  df_thc84$roll_mean %>% var(na.rm = TRUE)
  
# 85
  df_thc %>% f_plot_pt(85)
  df_thc85 <- df_thc %>% f_segment_signal(85, 0, 2000) # no blips!
  df_thc85$roll_mean %>% var(na.rm = TRUE)
  
# 86
  df_thc %>% f_plot_pt(86)
  df_thc86 <- df_thc %>% f_segment_signal(86, 0, 2000) # no blips!
  df_thc86$roll_mean %>% var(na.rm = TRUE)
  
# 87
  df_thc %>% f_plot_pt(87)
  df_thc87 <- df_thc %>% f_segment_signal(87, 0, 2000) # no blips!
  df_thc87$roll_mean %>% var(na.rm = TRUE)
  
# 88
  df_thc %>% f_plot_pt(88)
  df_thc88 <- df_thc %>% f_segment_signal(88, 0, 2000) # no blips!
  df_thc88$roll_mean %>% var(na.rm = TRUE)
  
# 89
  df_thc %>% f_plot_pt(89)
  df_thc89 <- df_thc %>% f_segment_signal(89, 0, 2000) # no blips!
  df_thc89$roll_mean %>% var(na.rm = TRUE)
  
# 90
  df_thc %>% f_plot_pt(90)
  df_thc90 <- df_thc %>% f_segment_signal(90, 0, 2000) # no blips!
  df_thc90$roll_mean %>% var(na.rm = TRUE)
  
# 91
  df_thc %>% f_plot_pt(91)
  df_thc91 <- df_thc %>% f_segment_signal(91, 0, 2000) # no blips!
  df_thc91$roll_mean %>% var(na.rm = TRUE)
  
# 92
  df_thc %>% f_plot_pt(92)
  df_thc92 <- df_thc %>% f_segment_signal(92, 0, 2000) # no blips!
  df_thc92$roll_mean %>% var(na.rm = TRUE)

# 93
  df_thc %>% f_plot_pt(93)
  df_thc93 <- df_thc %>% f_segment_signal(93, 0, 2000) # no blips!
  df_thc93$roll_mean %>% var(na.rm = TRUE)
  
# 94
  df_thc %>% f_plot_pt(94)
  df_thc94 <- df_thc %>% f_segment_signal(94, 0, 2000) # no blips!
  df_thc94$roll_mean %>% var(na.rm = TRUE)
  
# 95
  df_thc %>% f_plot_pt(95)
  df_thc95 <- df_thc %>% f_segment_signal(95, 0, 2000) # no blips!
  df_thc95$roll_mean %>% var(na.rm = TRUE)
  
# 96
  df_thc %>% f_plot_pt(96)
  df_thc96 <- df_thc %>% f_segment_signal(96, 0, 2000) # no blips!
  df_thc96$roll_mean %>% var(na.rm = TRUE)

# 97
  df_thc %>% f_plot_pt(97)
  df_thc97 <- df_thc %>% f_segment_signal(97, 0, 2000) # no blips!
  df_thc97$roll_mean %>% var(na.rm = TRUE)
  
# 98
  df_thc %>% f_plot_pt(98)
  df_thc98 <- df_thc %>% f_segment_signal(98, 0, 2000) # no blips!
  df_thc98$roll_mean %>% var(na.rm = TRUE)
  
# 99
  df_thc %>% f_plot_pt(99)
  df_thc99 <- df_thc %>% f_segment_signal(99, 0, 2000) # no blips!
  df_thc99$roll_mean %>% var(na.rm = TRUE)


# combine all  
  df_thc_filt <- rbind(df_thc1, df_thc2, df_thc3, df_thc6, df_thc7, df_thc8, df_thc9, df_thc10, #4, 5,  and 6 are all the exact same
                       df_thc11, df_thc12, df_thc13, df_thc14, df_thc15, df_thc16, df_thc17, df_thc18, df_thc19, df_thc20,
                       df_thc21, df_thc22, df_thc23, df_thc25, df_thc26, df_thc27, df_thc28,
                       df_thc31, df_thc32, df_thc33, df_thc34, df_thc35, df_thc36, df_thc37, df_thc38, df_thc39, df_thc40,
                       df_thc41, df_thc42, df_thc43, df_thc44, df_thc45, df_thc46, df_thc47, df_thc48, df_thc49, df_thc50,
                       df_thc51, df_thc52, df_thc53, df_thc54, df_thc55, df_thc56, df_thc57, df_thc58, df_thc59, df_thc60,
                       df_thc61, df_thc62, df_thc63, df_thc64, df_thc65, df_thc66, df_thc67, df_thc68, df_thc70,
                       df_thc71, df_thc72, df_thc73, df_thc74, df_thc75, df_thc76, df_thc77, df_thc78, df_thc79, df_thc80,
                       df_thc81, df_thc82, df_thc83, df_thc84, df_thc85, df_thc86, df_thc87, df_thc88, df_thc89, df_thc90,
                       df_thc91, df_thc92, df_thc93, df_thc94, df_thc95, df_thc96, df_thc97, df_thc98, df_thc99)
  

  
  
     
# calculate THC and THC sd ------------------------------------------------

  # total THC
  df_THCmean <- df_thc_filt %>% 
    filter(!is.na(roll_mean)) %>% 
    group_by(number, Status) %>% 
    summarise(THC_mean = mean(roll_mean))
  
  df_THCmean  %>%  
    ggplot(aes(x = Status, y = THC_mean)) +
    geom_violin(aes(color = Status)) +
    geom_jitter(position = position_jitter(0.2), shape = 1) +
    stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
    theme_bw() +
    theme(legend.position = "none") 
    
  
  # THC sd --> not included in paper
  df_THCsd <- df_thc_filt %>% 
    filter(!is.na(roll_mean)) %>% 
    group_by(number, Status) %>% 
    summarise(THC_sd = sd(roll_mean))
  
  df_THCsd  %>%  
    ggplot(aes(x = Status, y = THC_sd)) +
    geom_violin(aes(color = Status)) +
    geom_jitter(position = position_jitter(0.2), shape = 1) +
    stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
    theme_bw() +
    theme(legend.position = "none") 
  
  

# savitzy-golay signal smoothing  ---------------------------------------------------------

  
  # 4th order savitzy-golay filter
  library(signal)
  
  df_thc_filt <- df_thc_filt %>% 
    dplyr::select(-roll_mean) %>% 
    group_by(number) %>% 
    mutate(sg_filt = sgolayfilt(THC, p = 3)) %>% 
    ungroup()
  
  # plot
  df_thc_filt %>% filter(number == 10 & Time < 1500 & Time > 1490) %>% 
    ggplot(aes(x = Time)) +
    geom_line(aes(y = THC)) +
    geom_line(aes(y = sg_filt), color = "red")
  
  


# FFT ---------------------------------------------------------------------


  f_plot_fft <- function(df, num) {
    
    # compute fft
    
    fft <- df %>% filter(number == num) %>% .$sg_filt %>% fft() 
    
    # plot power spectra
    
    freq <- 50  #sample frequency in Hz 
    duration <- df %>% filter(number == num) %>% nrow()/freq # length of signal in seconds
    amo <- Mod(fft)
    freqvec <- 1:length(amo)
    
    freqvec <- freqvec/duration 
    df <- tibble(freq = freqvec, power = amo)
    df <- df[(1:as.integer(0.5*freq*duration)),]
    
    df %>% filter(freq > 0.01 & freq < 1.5) %>% 
      
      ggplot(aes(x = freq, y = power)) + 
      geom_line(stat = "identity") +
      geom_vline(xintercept = 0.01, lty = 2, color = "red", alpha = 0.5) +
      geom_vline(xintercept = 0.1, lty = 2, color = "red", alpha = 0.5) +
      
      #scale_x_continuous(breaks = seq(0, 1.5, 0.1)) +
      theme_bw() +
      
      theme(axis.text.x = element_text(angle = 45, hjust = 0.9))

  }
  
 f_plot_fft(df_thc_filt_pass, 71)   
    
 

 

# band-pass filter --------------------------------------------------------

 # filter each signal to below 0.1 Hz (physiological oscillations, including hemodynamics/cerebral autoregulation)
 # & above 0.01 (anything below is head displacements and motion noise) **** to include or not to include????
 
 library(seewave)
 library(dplR)
 library(data.table)
 
 # df_thc_filt <- df_thc_filt_pass[,1:5] # uncomment to recreate df_thc_filt
 
 l_thc_filt_pass <- unique(df_thc_filt$number) %>% 
   purrr::map(~dplyr::filter(df_thc_filt, number == .x, !is.na(sg_filt))) %>% 
   purrr::map(~mutate(.x, band_pass = pass.filt(.x$sg_filt, W = c(0.01, .1), type = "pass")))
 
 df_thc_filt_pass <- rbindlist(l_thc_filt_pass) %>% as_tibble()
 
 # plot
 df_thc_filt_pass %>% filter(number == 10 & Time < 1500 & Time > 1490) %>% 
   ggplot(aes(x = Time)) +
   geom_line(aes(y = THC)) +
   geom_line(aes(y = sg_filt), color = "red") +
   geom_line(aes(y = band_pass), color = "blue")

 

# Hilbert transformation --------------------------------------------------

 
  
# final signal processing -----------------------------------------------------------------------

 # add milliseconds in to Time
 df_thc_filt_pass <- df_thc_filt_pass %>% 
   group_by(number, Time) %>% 
   left_join(count(.)) %>% 
   mutate(Time = ifelse(n == 49, Time + 0.02*row_number(), Time + 0.02*(row_number() - 1))) %>% 
   dplyr::select(-n)
 
 
 # save filtered signal
 save(df_thc_filt_pass, file = "filtered_signals.Rdata")
 
 load("filtered_signals.Rdata")
 
 

 
