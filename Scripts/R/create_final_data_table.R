

library(tidyverse)

# LOAD DATA
df_hbox <- read_csv("Data/Hans_datatable_exports/malawi key data v04Jan2022.csv") %>% 
  as_tibble() %>% 
  janitor::clean_names() %>% 
  mutate(status = replace(status, status == "HV", "HC")) %>% 
  dplyr::rename("subject_id" = "subject_id_session") %>% 
  
  dplyr::filter(status %in% c("CM", "HC", "UM") & session_no == 1 | subject_id == "TM0003CM01") %>% 
  dplyr::filter(!grepl("blood", subject_id)) %>% 
  mutate(status = factor(status, levels = c("HC", "UM", "CM"))) %>% 
  mutate_at(c("bp_dias", "bp_sys", "glucose", "lactate"), as.numeric) %>% # converts remaining character variables to numerics
  select(subject_id, status, everything()) %>% 
  group_by(status)

# SUBSET FOR VARIABLES USED IN ANALYSIS
voi <- c(
  "subject_id", "status", "sex",
  "age_calc", "temperature", "hr", "bp_sys", "bp_dias", "resp_rate", "pulse_ox", "hct", "lactate", "parasites", "glucose", "aa_arg",
  "brain_swell",
  "cerebral_hb_tot", "cerebral_hb_tot_alpha2", "cerebral_hb_oxy", "cerebral_hb_oxy_alpha2",
  "muscle_hb_tot", "muscle_hb_tot_alpha2", "muscle_hb_oxy", "muscle_hb_oxy_alpha2"
)

df_final <- df_hbox %>% 
  dplyr::select_at(voi)

# EXPORT
write.csv(df_final, file = "data_for_upload/12June2023_key_data.csv", row.names = FALSE)

# QUANTIFY MISSING DATA

df_missing_vals <- tibble(
  
  status = c("HC", "UM", "CM")
  
) %>% 
  
  mutate(NAs = 
           
           map(
             
             .x = status,
             .f = ~ df_final %>% 
               ungroup %>% 
               filter(status == .x) %>% 
               dplyr::select(-subject_id, -status) %>% 
               is.na %>% 
               colSums() %>% 
               enframe
             
           )
         
  ) %>% 
  unnest(cols = c(NAs)) %>% 
  pivot_wider(id_cols = status, names_from = name, values_from = value)


write.csv(df_missing_vals, file = "data_for_upload/frequency_of_missing_vals.csv", row.names = FALSE)

