
### Code to generate Tables 2 & 3 Logistic regression model results


# libraries ---------------------------------------------------------------

library(tidyverse)
library(tidymodels)

# data --------------------------------------------------------------------


load("objects/20221125_final_glm_models.Rdata") # models compared and results generated in in scripts/markdown_files/20221115_revise_model.R

# table 2 -----------------------------------------------------------------

df_CI <- confint(hbtot_final2$fit) %>% exp() %>% as_tibble(rownames = "term")

hbtot_final2 %>%
  tidy(exponentiate = TRUE) %>%
  filter(p.value < 0.05) %>% 
  left_join(df_CI, by = "term") %>%
  mutate(estimate = ifelse(estimate < 10^-2, scientific(estimate, digits = 3), round(estimate, 2)),
         conf.low = ifelse(`2.5 %` < 10^-2, scientific(`2.5 %`, digits = 3), round(`2.5 %`, 2)),
         conf.high = ifelse(`97.5 %` < 10^-2, scientific(`97.5 %`, digits = 3), round(`97.5 %`, 2))) %>%
  mutate(p.value = scientific(p.value, digits = 3),
         statistic = round(statistic, 2)) %>%
  unite("95% CI", c(conf.low, conf.high), sep = ", ") %>%
  dplyr::select(c(term, estimate, `95% CI`, statistic, p.value)) %>% 
  dplyr::rename("Term" = "term", "Odds ratio" = "estimate",
                "t-statistic" = "statistic", "p-value" = "p.value") %>%
  filter(Term != "(Intercept)")


# table 3 -----------------------------------------------------------------

df_CI <- confint(hbtot_alpha_final2$fit) %>% exp() %>% as_tibble(rownames = "term")

hbtot_alpha_final2 %>%
  tidy(exponentiate = TRUE) %>%
  filter(p.value < 0.05) %>% 
  left_join(df_CI, by = "term") %>%
  mutate(estimate = ifelse(estimate < 10^-2, scientific(estimate, digits = 3), round(estimate, 2)),
         conf.low = ifelse(`2.5 %` < 10^-2, scientific(`2.5 %`, digits = 3), round(`2.5 %`, 2)),
         conf.high = ifelse(`97.5 %` < 10^-2, scientific(`97.5 %`, digits = 3), round(`97.5 %`, 2))) %>%
  mutate(p.value = scientific(p.value, digits = 3),
         statistic = round(statistic, 2)) %>%
  unite("95% CI", c(conf.low, conf.high), sep = ", ") %>%
  dplyr::select(c(term, estimate, `95% CI`, statistic, p.value)) %>% 
  dplyr::rename("Term" = "term", "Odds ratio" = "estimate",
                "t-statistic" = "statistic", "p-value" = "p.value") %>%
  filter(Term != "(Intercept)")

