rm(list = ls())

library(ggplot2)
library(tidyverse)
library(tidybayes)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))

ext <- ".png"


### Country list -----
countries <- c(
  KHM = "Cambodia",
  KEN = "Kenya",
  LAO = "Lao People's Democratic Republic", 
  MWI = "Malawi", 
  PAK = "Pakistan", 
  PHL = "Philippines", 
  TZA = "United Republic of Tanzania", 
  UGA = "Uganda", #
  VNM = "Viet Nam", 
  ZMB = "Zambia"
)


countries_cs <- c(
  MWI = "Malawi", 
  PHL = "Philippines", 
  TZA = "United Republic of Tanzania", 
  UGA = "Uganda", #
  ZMB = "Zambia"
)


countries_lab <- c(
  "Cambodia", "Kenya", "Lao PDR", 
  "Malawi", "Pakistan", "Philippines", 
  "Tanzania", "Uganda" , "Viet Nam", "Zambia"
)
names(countries_lab) <- countries


### Meta
load(file = "out/PN_All.rdata")

g_PN <- Dur_sim_all %>%
  mutate(Def = factor(Def, levels = c("D1", "D2", "D3"))) %>%
  ggplot() +
  geom_pointrange(aes(x = Def, y = m, ymin = l, ymax = u, colour = Def, shape = Type)) +
  geom_point(data = Dur_data_all, aes(x = Def, y = m, shape = Type)) +
  scale_x_discrete("Delay definition") + 
  scale_y_continuous("Delay to case-detection, months", labels = function(x) x * 12, breaks = c(0, 1, 2, 4, 6, 8)) +
  scale_colour_discrete("Definition") +
  facet_grid(.~ISO) +
  expand_limits(y = 0) +
  theme(legend.position = "bottom")



ggsave(plot = g_PN, paste0("docs/figs/duration/PN", ext), width = 7.5, height = 4.5)
