

### Code to generate fig S2: walk-through detrended fluctuation analysis



# libraries ---------------------------------------------------------------

library(tidyverse)
library(segmented)
library(patchwork)


# data load (brain DFA results) -------------------------------------------

load("objects/brain_bandpass_filtered_dfa_results.Rdata")

# set plot theme ----------------------------------------------------------

theme_set(
  theme_bw() +
    theme(
      legend.position = "none",
      text = element_text(family = "Tahoma"),
      plot.title = element_text(size = 20),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 15)
    )
)

# SET PATH FOR SAVING FIGURE
my_path <- "manuscript drafts/revision 3 21Aug2023/final_figures/"

# fig S2a -----------------------------------------------------------------

df_figs2a <- read.table("Data/signal_segments/Hb_tot/TM0001CM01.txt")

df_figs2a %>% 
  
  ggplot(aes(x = Time)) +
  geom_line(aes(y = THC)) +
  labs(y = "Hemoglobin concentration (μM)", x = "Time (s)") +
  ggtitle("Raw signal")

# SAVE
ggsave(paste0(my_path, "FigS2a.pdf"),
       height = 6, width = 8)

# fig S2b -----------------------------------------------------------------

df_figs2b <- read.table("Data/filtered_signals/0.2s_mean_0.001_1_filt/Hb_tot/TM0001CM01.txt") %>% 
  mutate(time = row_number())

df_figs2b %>% 
  ggplot(aes(x = time)) +
  geom_line(aes(y = V1)) +
  labs(y = "Bandpass filtered signal", x = "Time") +
  ggtitle("Pre-processed signal")

# SAVE
ggsave(paste0(my_path, "FigS2b.pdf"),
       height = 6, width = 8)


# fig S2c -----------------------------------------------------------------

df_figs2c <- df_figs2b %>% 
  mutate(mean = mean(V1),
         subtract_mean = V1 - mean,
         cumsum = cumsum(subtract_mean))

df_figs2c %>% 
  
  ggplot(aes(x = time, y = cumsum)) +
  geom_line() +
  geom_vline(xintercept = 1, lty = 2, color = "black") +
  geom_vline(xintercept = 1000, lty = 2, color = "black") +
  
  labs(y = "Cumulative sum", x = "Time") +
  ggtitle("Cumulative signal")

# SAVE
ggsave(paste0(my_path, "FigS2c.pdf"),
       height = 6, width = 8)

# fig S2d -----------------------------------------------------------------

fit <- lm(cumsum ~ time, data = df_figs2c %>% dplyr::filter(time > 1 & time < 1000))

# Window 
p_figs2c_window <- df_figs2c %>% 
  dplyr::filter(time > 1 & time < 1000) %>% 
  
  ggplot(aes(x = time, y = cumsum)) +
  geom_line(size = 2) +
  geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2], color = "black", lty = 2, size = 2) +
  labs(y = "Cumulative signal", x = "") +
  theme(
    legend.position = "none",
    text = element_text(family = "Tahoma"),
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_blank()
  )

# Detrended
p_figs2c_detrended <- df_figs2c %>% 
  dplyr::filter(time > 1 & time < 1000) %>% 
  
  ggplot(aes(x = time, y = fit$residuals)) +
  geom_line(size = 2) +
  geom_hline(yintercept = 0, color = "black", lty = 2, size = 2) +
  
  labs(y = "Detrended signal", x = "Time") +
  theme(
    legend.position = "none",
    text = element_text(family = "Tahoma"),
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_blank()
  )

# Combine
p_figs2c_window / p_figs2c_detrended +
  plot_annotation(title = "Detrending the cumulative signal")

# SAVE
ggsave(paste0(my_path, "FigS2d.pdf"),
       height = 10, width = 8)


# fig S2e -----------------------------------------------------------------

df_figs2e <- l_dfa_hbtot[[1]][[1]]

# find slope of segmented line
df_figs2e <- df_figs2e %>% mutate(log10_avg_fluct = log10(avg_fluctuation), log10_window_size = log10(window_size))
fit_figs2e <- lm(log10_avg_fluct ~ log10_window_size, data = df_figs2e) 
segfit <- segmented(fit_figs2e)
slopes <- coef(segfit) # the coefficients are the differences in slope in comparison to the previous slope

# first line: 
#y = b0 + b1*x; y = intercept1 + slope1 * x

# second line:
# y = c0 + c1*x; y = intercept2 + slope2 * x

# At the breakpoint (break1), the segments b and c intersect: b0 + b1*x = c0 + c1*x
b0 <- slopes[[1]]
b1 <- slopes[[2]]

c1 <- slopes[[2]] + slopes[[3]]
break1 <- segfit$psi[[2]]

# Solve for c0 (intercept of second segment):
c0 <- b0 + b1 * break1 - c1 * break1

# PLOT
df_figs2e %>% 
  
  ggplot(aes(x = log10_window_size, y = log10_avg_fluct)) +
  geom_point(size = 3, shape = 1) +
  geom_abline(intercept = c0, slope = c1, color = "black", lty = 2) +
  annotate("text", x = -1, y = -2, label = paste0("alpha = ", round(c1, 2)), size = 8, vjust = 1, hjust = 0) +
  
  labs(x = "log10(window length)", y = "log10(average fluctuation)") +
  ggtitle("DFA results")

# SAVE
ggsave(paste0(my_path, "FigS2e.pdf"),
       height = 6, width = 8)





