
### Code to generate Figure 1: cerebral and muscle Hb_tot and Hb_oxy in each patient group


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

# EXTRACT HB_TOT AND HB_OXY FOR THE BRAIN AND MUSCLE
df_fig1 <- df_hbox %>% 
  select(contains(c("hb_tot","hb_oxy")), -contains(c("alpha", "bsl"))) %>%
  ungroup() %>% 
  pivot_longer(2:ncol(.), names_to = "variable", values_to = "value") %>% 
  separate(variable, into = c("tissue", "variable"), sep = "_", extra = "merge") %>% 
  
  mutate(
    tissue = ifelse(tissue == "cerebral", "Brain", "Muscle"),
    variable = ifelse(variable == "hb_oxy", "Hb[oxy]", "Hb[tot]"
    )
  )


# CREATE COLOR VECTOR
my_col <- c("CM" = "royalblue2", "UM" = "goldenrod2", "HC" = "tomato2")

# SET THEME
theme_set(
  theme_bw() +
    theme(
      legend.position = "none",
      text = element_text(family = "Tahoma"),
      axis.title = element_text(size = 15),
      axis.title.y = element_blank(),
      axis.text = element_text(size = 15),
      strip.text.x = element_text(size = 20),
      strip.text.y = element_text(size = 20, vjust = 0),
      strip.background = element_blank(),
      strip.placement = "outside"
    )
)

# PLOT
df_fig1 %>% 
  
  ggplot(aes(x = status, y = value, color = status)) +
  geom_violin(linewidth = 1, trim = FALSE) +
  #geom_boxplot(width = 0.3) +
  geom_point(position = position_jitter(width = 0.15), 
             color = "black", shape = 1, size = 2) +
  stat_summary(fun = "median", geom = "crossbar", 
               linewidth = 0.3, width = 0.5) +
  facet_grid(variable ~ tissue, scales = "free_y", switch = "y",
             labeller = label_parsed
  ) +
  scale_color_manual(values = my_col) +
  labs(x = "")

# SAVE AS PDF
my_path <- "manuscript drafts/revision 3 21Aug2023/final_figures/"
ggsave(paste0(my_path, "Fig1.pdf"),
       height = 6, width = 10)
