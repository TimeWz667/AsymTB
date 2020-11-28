rm(list = ls())

library(ggplot2)
library(tidyverse)
library(tidybayes)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))

ext <- ".png"


### Country list -----
source("data/country_list.R")


### Meta
load(file = "out/PN_All.rdata")

g_PN <- rbind(Dur_sim_all %>% mutate(Year = 2019), 
              Dur_data_all %>% filter(Year == 2019)) %>%
  mutate(Def = factor(Def, levels = c("D1", "D2"))) %>%
  ggplot() +
  geom_pointrange(aes(x = ISO, y = m, ymin = l, ymax = u, colour = Type), position = position_dodge2(0.3), size = rel(0.4)) +
  scale_x_discrete("Country") + 
  scale_y_continuous("Delay to case-detection, months", labels = function(x) x * 12, breaks = c(0, 1, 2, 4, 6, 8, 10, 12)) +
  facet_grid(Def~., scales = "free_y", 
             labeller = labeller(Def = c(D1 = "From TB-detectable", D2 = "From symptom onset"))) +
  scale_colour_hue("", h = c(40, 200), labels = c("Duration by P:N ratio", "Duration by model prediction")) +
  expand_limits(y = 0) +
  theme(legend.position = "bottom")
g_PN


ggsave(plot = g_PN, paste0("docs/figs/duration/PN", ext), width = 7.5, height = 6.5)
