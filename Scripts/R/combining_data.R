

### Generate master datatable by combining data across spreadsheets


# libraries ---------------------------------------------------------------


  library(tidyverse)
  library(readxl)



# load data ---------------------------------------------------------------


  df_hbox <- read_xlsx("Data/Hans_datatable_exports/malawi key data v13Dec2021.xlsx") %>% as_tibble() %>% 
    mutate(Status = ifelse(Status == "HV", "HC", Status)) # clinical data
  
  df_bp_alphas_muscle <- read.csv("Data/final_dfa_results/bandpass_filtered_alphas_muscle.csv") %>% as_tibble() # muscle DFA res
  df_bp_alphas_brain <- read.csv("Data/final_dfa_results/bandpass_filtered_alphas.csv") %>% as_tibble() # brain DFA res
  df_bp_alphas_brain_missing_um <- read.csv("Data/final_dfa_results/cerebral_missing_UM_alphas.csv") %>% as_tibble()
  
  # muscle averages
  
      # Hb_tot
      hbtot_muscle_path <- "Data/signal_segments/muscle/Hb_tot"
      hbtot_muscle_files <- list.files(hbtot_muscle_path)
      l_hbtot_filt_muscle <- 1:length(hbtot_muscle_files) %>% 
        purrr::map(~read.table(paste0(hbtot_muscle_path, sep = "/", hbtot_muscle_files[.x])))
      
      df_med_hbtot_muscle <- 1:length(l_hbtot_filt_muscle) %>% 
        purrr::map(~mutate(l_hbtot_filt_muscle[[.x]], muscle_hb_tot = median(Hb_tot)) %>% 
                     dplyr::select(c(subject_id, muscle_hb_tot)) %>% 
                     unique()) %>% 
        rbindlist() %>% as_tibble()
      
      # Hb_oxy
      hboxy_muscle_path <- "Data/signal_segments/muscle/Hb_oxy"
      hboxy_muscle_files <- list.files(hboxy_muscle_path)
      l_hboxy_filt_muscle <- 1:length(hboxy_muscle_files) %>% 
        purrr::map(~read.table(paste0(hboxy_muscle_path, sep = "/", hboxy_muscle_files[.x])))
      
      df_med_hboxy_muscle <- 1:length(l_hboxy_filt_muscle) %>% 
        purrr::map(~mutate(l_hboxy_filt_muscle[[.x]], muscle_hb_oxy = median(Hb_oxy)) %>% 
                     dplyr::select(c(subject_id, muscle_hb_oxy)) %>% 
                     unique()) %>% 
        rbindlist() %>% as_tibble()
      
    # brain averages
      
      # Hb_tot
      hbtot_brain_path <- "Data/signal_segments/Hb_tot"
      hbtot_brain_files <- list.files(hbtot_brain_path)
      l_hbtot_filt_brain <- 1:length(hbtot_brain_files) %>% 
        purrr::map(~read.table(paste0(hbtot_brain_path, sep = "/", hbtot_brain_files[.x])))
      
      df_med_hbtot_brain <- 1:length(l_hbtot_filt_brain) %>% 
        purrr::map(~mutate(l_hbtot_filt_brain[[.x]], cerebral_hb_tot = median(THC)) %>% 
                     dplyr::select(c(subject_id, cerebral_hb_tot)) %>% 
                     unique()) %>% 
        rbindlist() %>% as_tibble()
      
      # Hb_oxy
      hboxy_brain_path <- "Data/signal_segments/Hb_oxy"
      hboxy_brain_files <- list.files(hboxy_brain_path)
      l_hboxy_filt_brain <- 1:length(hboxy_brain_files) %>% 
        purrr::map(~read.table(paste0(hboxy_brain_path, sep = "/", hboxy_brain_files[.x])))
      
      df_med_hboxy_brain <- 1:length(l_hboxy_filt_brain) %>% 
        purrr::map(~mutate(l_hboxy_filt_brain[[.x]], cerebral_hb_oxy = median(O2_sat)) %>% 
                     dplyr::select(c(subject_id, cerebral_hb_oxy)) %>% 
                     unique()) %>% 
        rbindlist() %>% as_tibble()
      
      

# format and combine data -------------------------------------------------------------

  df_hbox %>% colnames()
      # fix Status and subject_id to match other dfs
      df_hbox <- df_hbox %>% 
        rename("subject_id" = "Subject ID and Visit") %>% 
        dplyr::filter(Status %in% c("CM", "HC", "UM") & session.no == 1 | subject_id == "TM0003CM01") %>% 
        dplyr::filter(!grepl("blood", subject_id)) %>% 
        mutate(Status = factor(Status, levels = c("HC", "UM", "CM"))) %>% 
        select(c(subject_id, Status, everything()))
      
      # remove old alphas, add averages & new alphas
      subjs <- df_bp_alphas_brain_missing_um$subject_id
      
      df_master <- df_hbox %>% 
        select(-c(contains("alpha"), cerebral_hb_oxy, cerebral_hb_tot)) %>% 
        
        # add muscle alphas & medians
        full_join(df_bp_alphas_muscle %>% 
                     select(-c(X, contains("breakpoint"))) %>% 
                    full_join(df_med_hboxy_muscle, by = "subject_id") %>% 
                    full_join(df_med_hbtot_muscle, by = "subject_id"), 
                   by = c("subject_id", "Status")) %>% 
        
        # add brain alphas & medians
        full_join(df_bp_alphas_brain %>% 
                    select(-c(X, contains("breakpoint"))) %>% 
                    rename_with(~paste0("cerebral_", .x), 3:8) %>% 
                    full_join(df_med_hboxy_brain, by = "subject_id") %>% 
                    full_join(df_med_hbtot_brain, by = "subject_id"), 
                  by = c("subject_id", "Status")) %>% 
        
        # replace missing UM values
        full_join(df_bp_alphas_brain_missing_um, by = c("subject_id", "Status")) %>%
        mutate(cerebral_hbtot_overall_a = ifelse(subject_id %in% subjs, cerebral_hbtot_overall_a.y, cerebral_hbtot_overall_a.x),
               cerebral_hbtot_short_a = ifelse(subject_id %in% subjs, cerebral_hbtot_short_a.y, cerebral_hbtot_short_a.x),
               cerebral_hbtot_long_a = ifelse(subject_id %in% subjs, cerebral_hbtot_long_a.y, cerebral_hbtot_long_a.x),
               cerebral_hboxy_overall_a = ifelse(subject_id %in% subjs, cerebral_hboxy_overall_a.y, cerebral_hboxy_overall_a.x),
               cerebral_hboxy_short_a = ifelse(subject_id %in% subjs, cerebral_hboxy_short_a.y, cerebral_hboxy_short_a.x),
               cerebral_hboxy_long_a = ifelse(subject_id %in% subjs, cerebral_hboxy_long_a.y, cerebral_hboxy_long_a.x),
               cerebral_hb_tot = ifelse(subject_id %in% subjs, cerebral_hb_tot.y, cerebral_hb_tot.x),
               cerebral_hb_oxy = ifelse(subject_id %in% subjs, cerebral_hb_oxy.y, cerebral_hb_oxy.x)) %>% 
        select(-contains(".y"), -contains(".x")) %>%
        
        # convert lactate to numeric
        mutate(Lactate = as.numeric(Lactate)) %>% 
        
        # remove duplicate subject_ids
        filter(!is.na(Status)) %>% 
        
        # reorganize
        select(c(subject_id, Status, session.no, order(colnames(.)))) %>% 
        select(-contains("breakpoint"))




# filter  ---------------------------------------------------------



  # Remove patients labeled as "UM" that instead meet the WHO criteria for severe malaria
  
    # in final analyses, remove UM patient with lactate > 5 

    # All UM patients have hematocrit > 15%
    # df_master %>% dplyr::filter(Status == "UM"  & Hematocrit <= 15)
  
  # remove duplicate patient ID
  df_master <- df_master %>% unique()
  


# double check ------------------------------------------------------------


  # any duplicated patients?
  df_master %>% count(subject_id) %>% dplyr::filter(n > 1)


  # in final analysis, remove any Hb_oxy < 20
  # df_master <- df_master %>% dplyr::filter(cerebral_hb_oxy > 20)
  
  
  
# export as csv -----------------------------------------------------------

  write.csv(df_master, file = "Data/master_datatable.csv")


  
