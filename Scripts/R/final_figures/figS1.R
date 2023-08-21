

### Code to generate fig S1: DFA proof-of-concept on simulated signals


# libraries ---------------------------------------------------------------

library(tidyverse)

# data load & format --------------------------------------------------------------------

load("objects/sims.Rdata")

wn_alpha <- 1:length(l_wn_dfa) %>% purrr::map_df(~l_wn_dfa[[.x]][[2]]) %>% dplyr::rename("white" = "log10(window_size)")
pn_alpha <- 1:length(l_pn_dfa) %>% purrr::map_df(~l_pn_dfa[[.x]][[2]]) %>% dplyr::rename("pink" = "log10(window_size)")
bn_alpha <- 1:length(l_bn_dfa) %>% purrr::map_df(~l_bn_dfa[[.x]][[2]]) %>% dplyr::rename("brown" = "log10(window_size)")


# plot --------------------------------------------------------------------

df_figs1 <- bind_cols(wn_alpha, pn_alpha, bn_alpha) %>%
  pivot_longer(1:3, names_to = "noise", values_to = "alpha")

# CREATE COLOR VECTOR
my_col <- c("white" = "white", "brown" = "sienna2", "pink" = "pink")

# SET THEME
theme_set(
  theme_dark() +
    theme(
      legend.position = "none",
      text = element_text(family = "Tahoma"),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 15)
    )
)

# PLOT
df_figs1 %>% 
  
  ggplot(aes(x = noise, y = alpha, color = noise)) +
  geom_violin(fill = "transparent") +
  geom_jitter(position = position_jitter(width = 0.1), shape = 1, size = 2, alpha = 0.25) +
  stat_summary(fun = "median", geom = "crossbar", aes(color = noise), size = 0.2, width = 0.5) +
  
  scale_color_manual(values = my_col) +
  labs(x = "noise type", y = "alpha (\u03b1)",
       title = "DFA pipeline simulation results") +
  scale_y_continuous(limits = c(0, 1.7))

# SAVE AS PDF
my_path <- "manuscript drafts/revision 3 21Aug2023/final_figures/"
ggsave(paste0(my_path, "FigS1.pdf"),
       height = 4, width = 5)
