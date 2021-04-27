library(ggplot2)
library(tidyverse)
library(tidybayes)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))

ext <- glue::as_glue(".png")


### Country list -----
source("data/country_list.R")


### Data
load("out/summary_data.rdata")
load("out/Durations_All.rdata")

## PN ratio vs. duration
g_PN <- Durations$Durations_All %>%
  group_by(Country, ISO) %>%
  summarise(across(c(Delay_All, Delay_S, PN_All, PN_S, 
                     Delay_Sp, Delay_SymSp, PN_Sp, PN_SymSp), 
                   list(M = mean, L = ~quantile(.x, 0.025), U = ~quantile(.x, 0.975)))) %>%
  pivot_longer(-c(Country, ISO)) %>%
  separate(name, c("Est", "Def", "Index"), "_") %>%
  pivot_wider(names_from = Index, values_from = value) %>%
  mutate(Smear = ifelse(endsWith(Def, "Sp"), "Sp", "BC"),
         Sym = ifelse(Def %in% c("All", "Sp"), "All", "Sym")) %>%
  ggplot() +
  geom_pointrange(aes(x = ISO, y = M, ymin = L, ymax = U, colour = Est), position = position_dodge2(0.7), size = rel(0.4)) +
  scale_x_discrete("Country") + 
  scale_y_continuous("Time to case-detection, months", 
                     labels = function(x) x * 12, breaks = c(0, 1, 2, 4, 6, 8, 10, 12)) +
  facet_grid(Sym~Smear,
             labeller = labeller(Smear = c(BC = "Regardless of smear status", Sp = "Smear-positive"),
                                 Sym = c(All = "Asymptomatic / symptomatic TB", Sym = "Symptomatic TB"))) +
  scale_colour_hue("", h = c(40, 200), labels = c(PN="Duration by P:N ratio", Delay="Duration by model prediction")) +
  expand_limits(y = 0) +
  theme(legend.position = "bottom") 

g_PN


ggsave(plot = g_PN, "docs/figs/duration/PN_All" + ext, width = 9.5, height = 7.5)




## Incidence

g_Inc <- Durations$Durations_All %>%
  group_by(Country, ISO) %>%
  summarise(across(c(IncR_All, IncR_S), 
                   list(M = mean, L = ~quantile(.x, 0.025), U = ~quantile(.x, 0.975)))) %>%
  left_join(Tab_dat %>% select(Country, starts_with("IncR"))) %>%
  rename(IncR_WHO_M = IncR_M, IncR_WHO_L = IncR_L, IncR_WHO_U = IncR_U) %>%
  pivot_longer(-c(Country, ISO)) %>%
  separate(name, c("Est", "Def", "Index"), "_") %>%
  mutate(Def = factor(Def, levels = c("All", "S", "WHO"))) %>%
  pivot_wider(names_from = Index, values_from = value) %>%
  ggplot() +
  geom_pointrange(aes(x = Def, y = M, ymin = L, ymax = U, colour = Def)) +
  scale_x_discrete("Incidence definition") + 
  scale_y_continuous("Incidence rate, per 100 000", labels = function(x) x * 1E5) +
  scale_colour_discrete("", labels = c(All = "Est. Asymptomatic TB", S = "Est. Symptomatic TB", WHO = "WHO estimates")) +
  facet_grid(.~ISO) +
  expand_limits(y = 0) +
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank()) 


ggsave(plot = g_Inc, "docs/figs/incidence/Inc_All" + ext, width = 7.5, height = 5.0)


## CDR


g_CDR <- Durations$Durations_All %>%
  group_by(Country, ISO) %>%
  summarise(across(c(CDR_A, CDR_S), 
                   list(M = mean, L = ~quantile(.x, 0.025), U = ~quantile(.x, 0.975)))) %>%
  left_join(Tab_dat %>% 
              mutate(CDR_WHO_M = CDR_M / 100, CDR_WHO_L = CDR_L / 100, CDR_WHO_U = CDR_U / 100) %>% 
              select(Country, starts_with("CDR_WHO"))) %>%
  pivot_longer(-c(Country, ISO)) %>%
  separate(name, c("Est", "Def", "Index"), "_") %>%
  mutate(Def = factor(Def, levels = c("A", "S", "WHO"))) %>%
  pivot_wider(names_from = Index, values_from = value) %>%
  ggplot() +
  geom_pointrange(aes(x = Def, y = M, ymin = L, ymax = U, colour = Def)) +
  scale_x_discrete("Incidence definition") + 
  scale_y_continuous("Case-detection ratio, %", labels = scales::percent, breaks = seq(0, 1, 0.25)) +
  scale_colour_discrete("", labels = c(A = "Est. Asymptomatic TB", S = "Est. Symptomatic TB", WHO = "WHO estimates")) +
  facet_grid(.~ISO) +
  expand_limits(y = c(0, 1)) +
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank())


ggsave(plot = g_CDR, "docs/figs/incidence/CDR_All" + ext, width = 7.5, height = 5.0)

