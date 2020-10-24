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



### Main ----
for (i in 1:length(countries)) {
  iso <- names(countries)[i]
  country <- countries[i]
  
  ### Load data ----
  load(paste0("out/TTE_", iso, ".rdata"))
  
  gs <- list()

  
  gs$g_Dist <- TTE %>%
    ggplot(aes(x = Duration * 12)) +
    stat_halfeye(aes(fill = Phase), .width = .95, size = 2/3, alpha = 0.8) +
    scale_y_discrete("Density") +
    scale_x_continuous("Time since TB-detectable, month", breaks = c(0, 4, 8, 12, 16, 24, 36)) + 
    scale_fill_discrete("End point", labels = c(A = "Symptom onset", S = "Care-seeking intention", C = "Notification")) +
    facet_grid(Sex~.) +
    theme(legend.position = "bottom", panel.grid.minor.x = element_blank()) +
    expand_limits(x = 0)
  
  gs$g_Asym <- TTE %>%
    filter(Phase == "A") %>%
    ggplot(aes(x = Duration * 12)) +
    stat_halfeye(aes(fill = Phase), .width = .95, size = 2/3, alpha = 0.8) +
    scale_y_discrete("Density") +
    scale_x_continuous("Time since TB-detectable, month", breaks = c(0, 4, 8, 12, 16, 24, 36)) + 
    scale_fill_discrete("End point", labels = c(A = "Symptom onset")) +
    facet_grid(Sex~.) +
    theme(legend.position = "bottom", panel.grid.minor.x = element_blank()) +
    expand_limits(x = 0)
  
  
  save(gs, file = paste0("out/g_TTE_", iso, ".rdata"))
  
  ggsave(plot = gs$g_Dist, paste0("docs/figs/duration/TTE_", iso, ext), width = 7.5, height = 4.5)
  ggsave(plot = gs$g_Asym, paste0("docs/figs/duration/Asym_", iso, ext), width = 7.5, height = 4.5)

}


### Meta
load(file = "out/TTE_All.rdata")


gs <- list()


gs$g_Dist <- TTE_all %>%
  filter(Sex == "Total") %>%
  ggplot(aes(x = Duration * 12, y = reorder(Country, DA))) +
  stat_halfeye(aes(fill = Phase), .width = .95, size = 2/3, alpha = 0.8) +
  scale_y_discrete("Country", labels = countries_lab) +
  scale_x_continuous("Time since TB-detectable, month", breaks = c(0, 4, 8, 12, 16, 24, 36)) + 
  scale_fill_discrete("End point", labels = c(A = "Symptom onset", S = "Care-seeking intention", C = "Notification")) +
  theme(legend.position = "bottom", panel.grid.minor.x = element_blank()) +
  expand_limits(x = 0)


gs$g_Dist_Sex <- TTE_all %>%
  filter(Sex != "Total") %>%
  ggplot(aes(x = Duration * 12, y = reorder(Country, DA))) +
  stat_halfeye(aes(fill = Phase), .width = .95, size = 2/3, alpha = 0.8) +
  scale_y_discrete("Country", labels = countries_lab) +
  scale_x_continuous("Time since TB-detectable, month", breaks = c(0, 4, 8, 12, 16, 24, 36)) + 
  scale_fill_discrete("End point", labels = c(A = "Symptom onset", S = "Care-seeking intention", C = "Notification")) +
  facet_grid(.~Sex) +
  theme(legend.position = "bottom", panel.grid.minor.x = element_blank()) +
  expand_limits(x = 0)


gs$g_Asym <- TTE_all %>%
  filter(Sex == "Total" & Phase == "A") %>%
  ggplot(aes(x = Duration * 12, y = reorder(Country, DA))) +
  stat_halfeye(aes(fill = Phase), .width = .95, size = 2/3, alpha = 0.8) +
  scale_y_discrete("Country", labels = countries_lab) +
  scale_x_continuous("Time since TB-detectable, month", breaks = c(0, 4, 8, 12, 16, 24, 36)) + 
  scale_fill_discrete("End point", labels = c(A = "Symptom onset")) +
  theme(legend.position = "bottom", panel.grid.minor.x = element_blank()) +
  expand_limits(x = 0)

gs$g_Asym_Sex <- TTE_all %>%
  filter(Sex != "Total" & Phase == "A") %>%
  ggplot(aes(x = Duration * 12, y = reorder(Country, DA))) +
  stat_halfeye(aes(fill = Phase), .width = .95, size = 2/3, alpha = 0.8) +
  scale_y_discrete("Country", labels = countries_lab) +
  scale_x_continuous("Time since TB-detectable, month", breaks = c(0, 4, 8, 12, 16, 24, 36)) + 
  scale_fill_discrete("End point", labels = c(A = "Symptom onset")) +
  facet_grid(.~Sex) +
  theme(legend.position = "bottom", panel.grid.minor.x = element_blank()) +
  expand_limits(x = 0)


save(gs, file = "out/g_TTE_All.rdata")

ggsave(plot = gs$g_Dist, paste0("docs/figs/duration/TTE_All", ext), width = 5.5, height = 4.5)
ggsave(plot = gs$g_Dist_Sex, paste0("docs/figs/duration/TTE_Sex_All", ext), width = 7.5, height = 4.5)

ggsave(plot = gs$g_Asym, paste0("docs/figs/duration/Asym_All", ext), width = 5.5, height = 4.5)
ggsave(plot = gs$g_Asym_Sex, paste0("docs/figs/duration/Asym_Sex_All", ext), width = 7.5, height = 4.5)

