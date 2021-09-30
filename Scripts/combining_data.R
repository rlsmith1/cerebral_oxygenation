



# libraries ---------------------------------------------------------------


  library(tidyverse)



# load data ---------------------------------------------------------------


  df_hbox <- read.csv("Data/98 first visit records with cerebral oxygenation stats_trimmed.csv") %>% as_tibble()
  load("results.Rdata")
  
  

# format data -------------------------------------------------------------

  
  # get rid of "Other"
  df_hbox <- df_hbox %>% 
    dplyr::filter(Status != "Other") %>% 
    mutate(Status = factor(Status, levels = c("HV", "UM", "CM")))
  
  # identify duplicated rows (TM0003, TM2001)
  # df_hbox <- df_hbox %>% dplyr::filter(duplicated(Subject.ID..NIAID.) == FALSE)
  
  # Convert age to numeric
  df_hbox <- df_hbox %>%
    mutate(Age = gsub(" Years, ", "_", Age)) %>%
    mutate(Age = gsub(" Months ", "", Age)) %>%
    separate(Age, into = c("Years", "Months"), sep = "_") %>%
    mutate(Years = as.numeric(Years)) %>%
    mutate(Months = as.numeric(Months)/12) %>%
    mutate(Age = Years + Months) %>%
    dplyr::select(-c(Years, Months))
  
  # Convert sex to factor
  df_hbox <- df_hbox %>% mutate(Sex = as.factor(Sex))
  
  # Trim for variables of interest (VOI)
  df_hbox_trim <- df_hbox %>% 
    dplyr::select(c(Subject.ID..NIAID., Status, Subject.ID.and.Visit, Admission.Date, DOB, Sex, Age, Glucose, Hematocrit, Lactate, Temperature,
                    Arginine.umol.L, Haptoglobin..mg.dl., Hemoglobin..uM., 
                    BP.Diastolic, BP.Systolic, O2.Sat, RR, HR)) %>% 
    dplyr::filter(!grepl("blood", Subject.ID.and.Visit)) %>% 
    dplyr::rename("subject_id" = "Subject.ID..NIAID.", "visit" = "Subject.ID.and.Visit") %>% 
    dplyr::mutate(visit = substr(visit, start = 9, stop = 10))
  

  
  
# combine DFA results from brain and muscle, initial visit & follow-up --------

  load("results.Rdata")
  load("follow_up_results.Rdata")
  load("muscle_dfa_results.Rdata")
  load("muscle_follow_up_results.Rdata")
  

  df_results <- df_results %>% mutate(subject_id = substr(subject_id, start = 1, stop = 6)) %>% select(-number)
  df_fu_results <- df_fu_results %>% mutate(subject_id = substr(subject_id, start = 1, stop = 6)) %>% select(-number)
  df_muscle_results <- df_muscle_results %>% mutate(subject_id = substr(subject_id, start = 1, stop = 6)) %>% select(-number)
  df_muscle_fu_results <- df_muscle_fu_results %>% mutate(subject_id = substr(subject_id, start = 1, stop = 6)) %>% select(-number)
  
  # combine brain data
  df_brain_dfa_results <- df_results %>% 
    left_join(df_fu_results, by = c("subject_id", "Status")) %>% 
    select(-visit)
  
  # combine muscle data
  df_muscle_dfa_results <- df_muscle_results %>% 
    left_join(df_muscle_fu_results, by = c("subject_id", "Status"))
  
  
  # combine all!
  df_all_dfa_results <- df_brain_dfa_results %>% 
    right_join(df_muscle_dfa_results, by = c("subject_id", "Status"))
  
  
# join clinical data with dfa results -------------------------------------


  df_master <- df_hbox_trim %>% 
    left_join(df_all_dfa_results, by = c("subject_id", "Status"))


# filter  ---------------------------------------------------------



  # Remove patients labeled as "UM" that instead meet the WHO criteria for severe malaria
  
    # remove UM patient with lactate > 5 
    df_master <- df_master %>% dplyr::filter(!(Status == "UM" & Lactate > 5))
    
    # All UM patients have hematocrit > 15%
    # df_master %>% dplyr::filter(Status == "UM"  & Hematocrit <= 15)
  
  # remove duplicate patient ID
  df_master <- df_master %>% unique()
  


# double check ------------------------------------------------------------


  # any duplicated patients?
  df_master %>% count(subject_id) %>% dplyr::filter(n > 1) # TM2001


  # any negative Hb o2 sat?
  df_master %>% dplyr::filter(avg_Hb_o2sat < 20)
  
  
  
  
# export as csv -----------------------------------------------------------

  write.csv(df_master, file = "Data/master_datatable.csv")

  
  

  
  
  
  
  
  
  
  
