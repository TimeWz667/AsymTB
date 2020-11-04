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
  KEN = "Kenya",
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

g_PN <- rbind(Dur_sim_all, 
              Dur_data_all %>% mutate(l = m, u = m)
              ) %>%
  mutate(Def = factor(Def, levels = c("D1", "D2", "D3")),
         Cat = paste0(Def, ", ", Type)) %>%
  filter(!is.na(m) & Type == "Fitted") %>%
  ggplot() +
  geom_pointrange(aes(x = Def, y = m, ymin = l, ymax = u), position = position_dodge2(0.9)) +
  scale_x_discrete("Delay definition") + 
  scale_y_log10("Delay to case-detection, months", labels = function(x) x * 12, breaks = c(0, 1, 2, 4, 6, 8)) +
#  scale_colour_hue("Definition", l = 30, c = 100) +
  facet_grid(.~ISO) +
  #coord_flip() + 
  expand_limits(y = 1) +
  theme(legend.position = "bottom")
g_PN


ggsave(plot = g_PN, paste0("docs/figs/duration/PN", ext), width = 7.5, height = 4.5)
