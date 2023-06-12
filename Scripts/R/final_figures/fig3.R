
### Code to generate Figure 1: cerebral and muscle Hb_tot alpha and Hb_oxy alpha in each patient group


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

# EXTRACT HB_TOT ALPHA AND HB_OXY ALPHA FOR THE BRAIN AND MUSCLE
df_fig2 <- df_hbox %>% 
  select(contains("alpha2")) %>%
  ungroup() %>% 
  pivot_longer(2:ncol(.), names_to = "variable", values_to = "value") %>% 
  separate(variable, into = c("tissue", "variable"), sep = "_", extra = "merge") %>% 
  mutate(variable = str_remove(variable, "2"))

# CREATE COLOR VECTOR
my_col <- c("CM" = "#00BFC4", "HC" = "#7CAE00", "UM" = "#F8766D")

# PLOT
df_fig2 %>% 
  
  ggplot(aes(x = status, y = value)) +
  geom_violin(aes(color = status)) +
  geom_jitter(position = position_jitter(0.2), shape = 1, size = 2) +
  stat_summary(fun = "median", geom = "crossbar", aes(color = status), size = 0.2, width = 0.5) +
  facet_grid(variable ~ tissue, scales = "free_y") +
  
  scale_color_manual(values = my_col) +
  labs(x = "") +
  theme_bw() +
  theme(legend.position = "none")  +
  ylim(c(0, 1.5)) +
  
  theme(# text = element_text(family = "Arial"),
    strip.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 15))

