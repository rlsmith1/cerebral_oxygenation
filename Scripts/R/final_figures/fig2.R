
### Code to generate Figure 2: correlation between brain swell score and Hb_tot


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


# plot --------------------------------------------------------------------

df_hbox %>%
  
  # log10 transform data that are not normal
  mutate(log10_cerebral_hb_tot = log10(cerebral_hb_tot)) %>% 
  
  # plot
  ggplot(aes(x = log10_cerebral_hb_tot, y = brain_swell)) +
  geom_point(aes(size = hct), shape = 1) +
  
  geom_smooth(method = "lm", lty = 2, color = "#696969", alpha = 0.2, se = FALSE) +
  theme_bw() +
  ggtitle("1. Relationship between brain swell score and log10-transformed cerebral Hb_tot") +
  theme(strip.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

