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



### Main ----
for (i in 1:length(countries)) {
  iso <- names(countries)[i]
  country <- countries[i]
  
  ### Load data ----
  load(paste0("out/Duration_", iso, ".rdata"))
  
  gs <- list()

  
  gs$g_Duration <- Durations %>%
    ggplot() + 
    stat_halfeye(aes(x = value * 12, y = Sex, fill = Sex), .width = .95, size = 2/3, alpha = 0.8) +
    scale_y_discrete("") +
    scale_x_continuous("Duration, month", breaks = c(0, 4, 8, 12, 16, 24, 36)) + 
    expand_limits(x = 0) +
    facet_grid(Stage ~., labeller = labeller(Stage = c(sym = "Asymptomatic TB", care = "Symptomatic TB"))) + 
    theme(legend.position = "none", panel.grid.minor.x = element_blank())
  
  gs$g_TTE <- TTE %>%
    ggplot() + 
    stat_halfeye(aes(x = value * 12, y = Sex, fill = Sex), .width = .95, size = 2/3, alpha = 0.8) +
    scale_y_discrete("") +
    scale_x_continuous("Time since TB-detectable, month", breaks = c(0, 4, 8, 12, 16, 24, 36)) + 
    expand_limits(x = 0) +
    facet_grid(Stage ~., labeller = labeller(Stage = c(sym = "Symptom onset", care = "Notification"))) + 
    theme(legend.position = "none", panel.grid.minor.x = element_blank())
  
  save(gs, file = paste0("out/g_Duration_", iso, ".rdata"))
  
  ggsave(plot = gs$g_Duration, paste0("docs/figs/duration/Duration_", iso, ext), width = 7.5, height = 4.5)
  ggsave(plot = gs$g_TTE, paste0("docs/figs/duration/TTE_", iso, ext), width = 7.5, height = 4.5)

}


### Meta
load(file = "out/Duration_All.rdata")


countries_lab <- c(
  "Cambodia", "Kenya", "Lao PDR", 
  "Malawi", "Pakistan", "Philippines", 
  "Tanzania", "Uganda" , "Viet Nam", "Zambia"
)
names(countries_lab) <- countries

gs <- list()


gs$g_Dist <- sim_dur %>%
  filter(Sex == "Total") %>%
  ggplot(aes(x = Duration * 12, y = reorder(Country, DA))) +
  stat_halfeye(aes(fill = Phase), .width = .95, size = 2/3, alpha = 0.8) +
  scale_y_discrete("Country", labels = countries_lab) +
  scale_x_continuous("Time since TB-detectable, month", breaks = c(0, 4, 8, 12, 16, 24, 36)) + 
  scale_fill_discrete("End point", labels = c(A = "Symptom onset", S = "Care-seeking intention", C = "Notification")) +
  theme(legend.position = "bottom", panel.grid.minor.x = element_blank()) +
  expand_limits(y = 0)


gs$g_Dist_Sex <- sim_dur %>%
  filter(Sex != "Total") %>%
  ggplot(aes(x = Duration * 12, y = reorder(Country, DA))) +
  stat_halfeye(aes(fill = Phase), .width = .95, size = 2/3, alpha = 0.8) +
  scale_y_discrete("Country", labels = countries_lab) +
  scale_x_continuous("Time since TB-detectable, month", breaks = c(0, 4, 8, 12, 16, 24, 36)) + 
  scale_fill_discrete("End point", labels = c(A = "Symptom onset", S = "Care-seeking intention", C = "Notification")) +
  facet_grid(.~Sex) +
  theme(legend.position = "bottom", panel.grid.minor.x = element_blank()) +
  expand_limits(y = 0)


save(gs, file = "out/g_Duration_All.rdata")

ggsave(plot = gs$g_Dist, paste0("docs/figs/duration/TTE_All", ext), width = 5.5, height = 4.5)
ggsave(plot = gs$g_Dist_Sex, paste0("docs/figs/duration/TTE_All_Sex", ext), width = 7.5, height = 4.5)

