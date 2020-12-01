rm(list = ls())

library(ggplot2)
library(tidyverse)
library(tidybayes)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))

ext <- glue::as_glue(".png")


### Country list -----
source("data/country_list.R")
load("out/Durations_all.rdata")

### Main ----
for (i in 1:length(countries)) {
  iso <- names(countries)[i]
  country <- countries[i]
  
  ### Load data ----
  TTE <- Durations$Durations_Sex %>%
    filter(ISO == iso) %>%
    select(Key, Sex, TTS, TTC, TTN) %>%
    pivot_longer(c("TTS", "TTC", "TTN"), names_to = "Phase", values_to = "Duration") %>%
    filter(Duration > 0)
  
  gs <- list()

  
  gs$g_Dist <- TTE %>%
    ggplot(aes(x = Duration * 12)) +
    stat_halfeye(aes(fill = Phase), .width = .95, size = 2/3, alpha = 0.8) +
    scale_y_discrete("Density") +
    scale_x_continuous("Time since TB-detectable, month", breaks = c(0, 4, 8, 12, 16, 24, 36)) + 
    scale_fill_discrete("End point", labels = c(TTS = "Symptom onset", TTC = "Care-seeking", TTN = "Notification")) +
    facet_grid(Sex~.) +
    theme(legend.position = "bottom", panel.grid.minor.x = element_blank()) +
    expand_limits(x = 0)
  
  gs$g_Asym <- TTE %>%
    filter(Phase == "TTS") %>%
    ggplot(aes(x = Duration * 12)) +
    stat_halfeye(aes(fill = Phase), .width = .95, size = 2/3, alpha = 0.8) +
    scale_y_discrete("Density") +
    scale_x_continuous("Time since TB-detectable, month", breaks = c(0, 4, 8, 12, 16, 24, 36)) + 
    scale_fill_discrete("End point", labels = c(TTS = "Symptom onset")) +
    facet_grid(Sex~.) +
    theme(legend.position = "bottom", panel.grid.minor.x = element_blank()) +
    expand_limits(x = 0)
  
  
  save(gs, file = paste0("out/g_TTE_", iso, ".rdata"))
  
  ggsave(plot = gs$g_Dist, paste0("docs/figs/duration/TTE_", iso, ext), width = 7.5, height = 4.5)
  ggsave(plot = gs$g_Asym, paste0("docs/figs/duration/Asym_", iso, ext), width = 7.5, height = 4.5)

}


### Meta
gs <- list()



TD <- Durations$Durations_All %>%
  group_by(Country, ISO) %>%
  summarise(m = mean(TTN), l = quantile(TTN, 0.025), u = quantile(TTN, 0.975), DA = mean(TTS)) %>%
  mutate(TD = sprintf("%i (%i-%i)", round(m * 12), round(l * 12), round(u * 12)), Key = 1) %>%
  arrange(DA) %>%
  select(-DA) %>% select(Country, TD, u, Key)


country_i <- 1:length(countries)
names(country_i) <- TD$Country


gs$g_Dist <- Durations$Durations_All %>%
  select(Country, ISO, Key, TTS, TTC, TTN) %>%
  mutate(DA = TTS, CountryI = country_i[Country]) %>%
  pivot_longer(c("TTS", "TTC", "TTN"), names_to = "Phase", values_to = "Duration") %>%
  mutate(Phase = factor(Phase, levels = c("TTS", "TTC", "TTN"))) %>%
  filter(Duration > 0 & Key < 200) %>%
  left_join(TD) %>%
  ggplot(aes(x = Duration * 12, y = reorder(Country, DA))) +
  stat_halfeye(aes(fill = Phase), .width = .95, size = 2/3, alpha = 0.8) +
  geom_text(aes(x = 52, label = TD), vjust = -.5, hjust = 0.5, size=3) +
  scale_y_discrete("Country", labels = countries_lab) +
  scale_x_continuous("Time since TB-detectable, month", 
                     breaks = c(0, 4, 8, 12, 16, 24, 36, 52),
                     labels = c(0, 4, 8, 12, 16, 24, 36, "Total duration (95% CrI)")) + 
  scale_fill_discrete("End point", labels = c(TTS = "Symptom onset", TTC = "Care-seeking", TTN = "Notification")) +
  theme(legend.position = "bottom", panel.grid.minor.x = element_blank()) +
  expand_limits(x = c(0, 56))


gs$g_Dist_Sex <- Durations$Durations_Sex %>%
  select(Country, ISO, Key, Sex, TTS, TTC, TTN) %>%
  mutate(DA = TTS) %>%
  pivot_longer(c("TTS", "TTC", "TTN"), names_to = "Phase", values_to = "Duration") %>%
  mutate(Phase = factor(Phase, levels = c("TTS", "TTC", "TTN"))) %>%
  filter(Duration > 0) %>%
  ggplot(aes(x = Duration * 12, y = reorder(Country, DA))) +
  stat_halfeye(aes(fill = Phase), .width = .95, size = 2/3, alpha = 0.8) +
  scale_y_discrete("Country", labels = countries_lab) +
  scale_x_continuous("Time since TB-detectable, month", breaks = c(0, 4, 8, 12, 16, 24, 36)) + 
  scale_fill_discrete("End point", labels = c(TTS = "Symptom onset", TTC = "Care-seeking", TTN = "Notification")) +
  facet_grid(.~Sex) +
  theme(legend.position = "bottom", panel.grid.minor.x = element_blank()) +
  expand_limits(x = 0)


gs$g_Asym <- Durations$Durations_All %>%
  select(Country, ISO, Key, TTS) %>%
  mutate(DA = TTS) %>%
  pivot_longer(c("TTS"), names_to = "Phase", values_to = "Duration") %>%
  filter(Duration > 0) %>%
  ggplot(aes(x = Duration * 12, y = reorder(Country, DA))) +
  stat_halfeye(aes(fill = Phase), .width = .95, size = 2/3, alpha = 0.8) +
  scale_y_discrete("Country", labels = countries_lab) +
  scale_x_continuous("Time since TB-detectable, month", breaks = c(0, 4, 8, 12, 16, 24, 36)) + 
  scale_fill_discrete("End point", labels = c(TTS = "Symptom onset")) +
  theme(legend.position = "bottom", panel.grid.minor.x = element_blank()) +
  expand_limits(x = 0)

gs$g_Asym_Sex <- Durations$Durations_Sex %>%
  select(Country, Sex, ISO, Key, TTS) %>%
  mutate(DA = TTS) %>%
  pivot_longer(c("TTS"), names_to = "Phase", values_to = "Duration") %>%
  filter(Duration > 0) %>%
  ggplot(aes(x = Duration * 12, y = reorder(Country, DA))) +
  stat_halfeye(aes(fill = Phase), .width = .95, size = 2/3, alpha = 0.8) +
  scale_y_discrete("Country", labels = countries_lab) +
  scale_x_continuous("Time since TB-detectable, month", breaks = c(0, 4, 8, 12, 16, 24, 36)) + 
  scale_fill_discrete("End point", labels = c(A = "Symptom onset")) +
  facet_grid(.~Sex) +
  theme(legend.position = "bottom", panel.grid.minor.x = element_blank()) +
  expand_limits(x = 0)


gs$g_Dur_Dur <- Durations$Durations_All %>%
  select(Country, ISO, Key, DurA, DurS, DurC) %>%
  filter(Key < 200) %>%
  mutate(DA = DurA, DT = DurS + DurC,Country = countries_lab[Country]) %>%
  ggplot(aes(x = DA * 12, y = DT * 12)) +
  geom_point(aes(colour = Country), alpha = 0.2) +
  scale_y_continuous("Symptomat onset to notification, month", breaks = c(0, 4, 8, 12, 16, 24, 36)) +
  scale_x_continuous("TB-detectable to symptom onset, month", breaks = c(0, 4, 8, 12, 16, 24, 36)) +
  facet_wrap(.~Country) +
  theme(legend.position = "none")


save(gs, file = "out/g_TTE_All.rdata")

ggsave(plot = gs$g_Dist, paste0("docs/figs/duration/TTE_All", ext), width = 5.5, height = 4.5)
ggsave(plot = gs$g_Dist, paste0("docs/figs/TTE", ext), width = 6.5, height = 4.5)
ggsave(plot = gs$g_Dist_Sex, paste0("docs/figs/duration/TTE_Sex_All", ext), width = 7.5, height = 4.5)

ggsave(plot = gs$g_Asym, paste0("docs/figs/duration/Asym_All", ext), width = 5.5, height = 4.5)
ggsave(plot = gs$g_Asym_Sex, paste0("docs/figs/duration/Asym_Sex_All", ext), width = 7.5, height = 4.5)
ggsave(plot = gs$g_Dur_Dur, paste0("docs/figs/duration/DurationDuration", ext), width = 6.5, height = 6.5)

