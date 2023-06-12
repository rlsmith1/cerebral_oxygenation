
### Code to generate Table 1 Clinical characteristics of the study population from raw data


# libraries ---------------------------------------------------------------

library(tidyverse)

# data load & format --------------------------------------------------------------------

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

# data format -------------------------------------------------------------

### TABLE FOR NUMBER OF SUBJECTS IN EACH PATIENT GROUP
df_status <- df_hbox %>% 
  count() %>% 
  pivot_wider(names_from = "status", values_from = "n") %>% 
  mutate(variable = "Number of subjects") %>% 
  dplyr::select(variable, everything()) %>% 
  mutate_if(is.numeric, as.character)

### TABLE FOR % FEMALE IN EACH PATIENT GROUP
df_sex <- df_hbox %>% 
  count(sex) %>% 
  pivot_wider(id_cols = "status", names_from = "sex", values_from = "n") %>% 
  mutate(perc_fem = round(Female/(Female + Male)*100, 2)) %>% 
  dplyr::select(-c(Female, Male)) %>% 
  
  pivot_wider(names_from = "status", values_from = "perc_fem") %>% 
  mutate(variable = "Sex (% female)") %>% 
  dplyr::select(variable, everything()) %>% 
  mutate_if(is.numeric, as.character)

### NUMERIC VARIABLES

# select variables for table
clinical_voi <- c("age_calc", "temperature", "hr", "bp_sys", "bp_dias", "resp_rate", "pulse_ox", "hct", "lactate", "glucose", "aa_arg")

# calculate median, first and third quartiles
df_clinical_res <- df_hbox %>% 
  select_at(clinical_voi) %>% 
  summarise_all(funs(median = median, q1 = quantile(., probs = 0.25), q3 = quantile(., probs = 0.75)), na.rm = TRUE) %>% 
  mutate_if(is.numeric, round, 2)

# write function to paste together values
f_unite_cols <- function(col) {
  df_clinical_res %>% 
    unite(!! col, starts_with(col), sep = "_") %>% 
    dplyr::select(!! col) 
}


# paste together and transpose for clinical data table --------------------


df_clinical_res <- 1:length(clinical_voi) %>% 
  purrr::map_dfc(~f_unite_cols(clinical_voi[.x])) %>% 
  mutate_all(funs(str_replace(., "_", " ("))) %>% 
  mutate_all(funs(str_replace(., "_", ", "))) %>% 
  mutate_all(funs(paste0(., ")"))) %>% 
  
  t() %>% 
  as_tibble(rownames = "variable") %>% 
  dplyr::rename("HC" = "V1", "UM" = "V2", "CM" = "V3") %>% 
  add_row(df_status, .before = 1) %>% 
  add_row(df_sex, .before = 2) 

# export to CSV
write.csv(df_clinical_res, file = "Outputs/table_1.csv")



# CM specific variables ---------------------------------------------------

### IDENTIFY AND SELECT VARIABLES OF INTEREST
cm_voi <- c("bcs_session", "brain_swell", "hrp2")

df_cm_voi <- df_hbox %>% 
  ungroup %>% 
  filter(status == "CM") %>% 
  select_at(cm_voi) %>% 
  mutate_at(c("bcs_session", "brain_swell"), as.factor) 

### BCS COUNTS
df_cm_voi %>% count(bcs_session)

# BRAIN SWELL SCORE COUNTS
df_cm_voi %>% count(brain_swell)

# HRP
df_cm_hrp2 <- df_cm_voi %>% 
  select(cm_voi[3]) %>% 
  summarise_all(funs(median = median, 
                     q1 = quantile(., probs = 0.25), 
                     q3 = quantile(., probs = 0.75)), 
                na.rm = TRUE) %>% 
  unite(hrp2, 1:2, sep = " (") %>% 
  unite(hrp2, 1:2, sep = ", ") %>% 
  mutate(hrp2 = paste0(hrp2, ")")) 

