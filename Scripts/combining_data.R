



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
    dplyr::select(c(Subject.ID.and.Visit, Status, Sex, Age, Glucose, Hematocrit, Lactate, Temperature,
                    Arginine.umol.L, Haptoglobin..mg.dl., Hemoglobin..uM., 
                    BP.Diastolic, BP.Systolic, O2.Sat, RR, HR)) %>% 
    dplyr::rename("subject_id" = "Subject.ID.and.Visit") 
  
  

# join clinical data with dfa results -------------------------------------


  
  df_master <- df_hbox_trim %>% 
    left_join(df_results, by = c("subject_id", "Status")) %>% 
    dplyr::select(-number)


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
  df_master %>% count(subject_id) %>% count(n)
  
  # any negative Hb o2 sat?
  df_master %>% dplyr::filter(avg_Hb_o2sat < 20)
  
  
# export as csv -----------------------------------------------------------

  write.csv(df_master, file = "Data/master_datatable.csv")
    
  

  
  
  
  
  
  
  
  
