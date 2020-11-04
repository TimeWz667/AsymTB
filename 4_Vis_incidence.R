rm(list = ls())

library(ggplot2)


theme_set(theme_bw() + theme(text = element_text(family = "sans")))


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

countries_iso <- names(countries)
names(countries_iso) <- countries

countries_lab <- c(
  "Cambodia", "Kenya", "Lao PDR", 
  "Malawi", "Pakistan", "Philippines", 
  "Tanzania", "Uganda" , "Viet Nam", "Zambia"
)
names(countries_lab) <- countries


ext <- ".png"


for (i in 1:length(countries)) {
  iso <- names(countries)[i]
  country <- countries[i]
  
  ### Load data ----
  load(paste0("out/Incidence_", iso, ".rdata"))

  g_Inc <- inc %>%
    mutate(Index = factor(Index, levels = c("A", "S", "WHO"))) %>%
    ggplot() +
    geom_pointrange(aes(x = Index, y = m, ymin = l, ymax = u, colour = Index), size = 1.2) +
    scale_x_discrete("Incidence definition", labels = c(A = "A", S = "S", WHO = "WHO")) + 
    scale_y_continuous("Incidence rate, per 100 000", labels = function(x) x * 1E5) +
    scale_colour_discrete("", labels = c(A = "Est. Asymptomatic TB", S = "Est. Symptomatic TB", WHO = "WHO estimates")) +
    facet_grid(.~Sex) +
    expand_limits(y = 0) +
    theme(legend.position = "bottom")
  
  save(g_Inc, file = paste0("out/g_Incidence_", iso, ".rdata"))
  
  ggsave(plot = g_Inc, paste0("docs/figs/incidence/Inc_", iso, ext), width = 5.5, height = 3.5)
}



load("out/Incidence_All.rdata")


g_Inc <- incs %>%
  filter(Sex == "Total") %>%
  mutate(Index = factor(Index, levels = c("A", "S", "WHO"))) %>%
  ggplot() +
  geom_pointrange(aes(x = Index, y = m, ymin = l, ymax = u, colour = Index)) +
  scale_x_discrete("Incidence definition", labels = c(A = "A", S = "S", WHO = "WHO")) + 
  scale_y_continuous("Incidence rate, per 100 000", labels = function(x) x * 1E5) +
  scale_colour_discrete("", labels = c(A = "Est. Asymptomatic TB", S = "Est. Symptomatic TB", WHO = "WHO estimates")) +
  facet_grid(.~Country, labeller = labeller(Country = countries_iso)) +
  expand_limits(y = 0) +
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank())

save(g_Inc, file = paste0("out/g_Incidence_All", ext))
ggsave(plot = g_Inc, paste0("docs/figs/incidence/Inc_All", ext), width = 7.5, height = 5.0)
