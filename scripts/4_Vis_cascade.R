rm(list = ls())

library(ggplot2)
library(tidyverse)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))

ext <- glue::as_glue(".png")


### Country list -----
source("data/country_list.R")


colour_as <- c(Death = "#636363",
               LTFU = "#969696",
               SelfCured = "#d9d9d9",
               Cured = "#bcbddc",
               Tr = "#4BB35D",
               S = "#fd8d3c",
               A = "#9ecae1")
  
labels_as <- c(SelfCured = "Self Cured",
               Cured = "Completed",
               Tr = "On treatment",
               S = "Sym.TB",
               A = "Asym. TB")

colour_asc <- c(Death = "#636363",
                LTFU = "#969696",
                SelfCured = "#d9d9d9",
                Cured = "#bcbddc",
                Tr = "#4BB35D",
                C = "#fd8d3c",
                S = "#fdd0a2",
                A = "#9ecae1")

labels_asc <- c(SelfCured = "Self Cured",
               Cured = "Completed",
               Tr = "On treatment",
               C = "Sym. with\n care-seeking",
               S = "Sym. without\n care-seeking",
               A = "Asym. TB")


colour_co_asc = c(complete = "#bcbddc",
           care = "#4BB35D",
           aware = "#fd8d3c",
           sym = "#fdd0a2", 
           inc = "#9ecae1")
labels_co_asc = c(complete = "Treatment complete",
           care = "Treatment start",
           aware = "Care-seeking start",
           sym = "Symptom onset",
           inc = "Incidence")

colour_co_as = c(complete = "#bcbddc",
                  care = "#4BB35D",
                  sym = "#fdd0a2", 
                 inc = "#9ecae1")
labels_co_as = c(complete = "Treatment complete",
                  care = "Treatment start",
                  sym = "Symptom onset",
                 inc = "Incidence")


for (i in 1:length(countries)) {
  iso <- glue::as_glue(names(countries)[i])
  country <- countries[i]
  
  ### Load data ----
  load("out/Cascade_" + iso + ".rdata")

  gs <- list()
  
  if (country %in% countries_cs) {
    col <- colour_asc
    lab <- labels_asc 
    
    col_co <- colour_co_asc
    lab_co <- labels_co_asc
  } else {
    col <- colour_as
    lab <- labels_as
    
    col_co <- colour_co_as
    lab_co <- labels_co_as
  }
  
  end <-  Cohort$End %>% filter(Sex != "Total") %>% mutate(t = 2.2)
  
  
  gs$g_Cohort_Sex <- Cohort$Cohort %>% filter(Sex != "Total") %>%
    ggplot() +
    geom_bar(aes(x = t, y = value, fill = name), stat = "identity", width = 0.02) +
    geom_bar(data = end, aes(x = t, y = value, fill = name), stat = "identity", width = 0.1) +
    #facet_grid(Sex~.) +
    scale_fill_manual("State", 
                      values = col,
                      labels = lab) +
    scale_x_continuous("Time since TB-detectable, year", breaks = c(0.5, 1, 1.5, 2, 2.2), labels = c(0.5, 1, 1.5, 2, "End")) +
    scale_y_continuous("State Distribution, %", labels = scales::percent) + 
    facet_grid(Sex~.) +
    labs(title = country)
  
  end <-  Cohort$End %>% filter(Sex == "Total") %>% mutate(t = 2.2)
  
  gs$g_Cohort_Total <-  Cohort$Cohort %>% filter(Sex == "Total") %>%
    ggplot() +
    geom_bar(aes(x = t, y = value, fill = name), stat = "identity", width = 0.02) +
    geom_bar(data = end, aes(x = t, y = value, fill = name), stat = "identity", width = 0.1) +
    scale_fill_manual("State", 
                      values = col,
                      labels = lab) +
    scale_x_continuous("Time since TB-detectable, year", breaks = c(0.5, 1, 1.5, 2, 2.2), labels = c(0.5, 1, 1.5, 2, "End")) +
    scale_y_continuous("State Distribution, %", labels = scales::percent) + 
    facet_grid(Sex~.) +
    labs(title = country)


  save(gs, file = "out/g_Cascade_" + iso + ".rdata")
  
  ggsave(plot = gs$g_Cohort_Sex, "docs/figs/Cascade/CohortSex_" + iso + ext, width = 7.5, height = 7.5)
  ggsave(plot = gs$g_Cohort_Total, "docs/figs/Cascade/Cohort_" + iso + ext, width = 7.5, height = 4.5)
  #ggsave(plot = gs$g_Cascade, "docs/figs/Cascade/Cascade_" + iso + ext, width = 6.5, height = 4.5)
}





load("out/Cascade_All.rdata")


dat <- Cascade$Cascade_All %>%
  group_by(Country, ISO) %>%
  summarise(across(everything(), list(M = mean, 
                                    L = ~quantile(., 0.025),
                                    U = ~quantile(., 0.975)))) %>%
  pivot_longer(starts_with(c("C_", "G_")), values_to = "pr") %>%
  separate(name, c("Index", "Stage", "Stat"), "_") %>%
  select(-Key_M, -Key_L, -Key_U) %>%
  pivot_wider(c(Country, ISO, Index, Stage), values_from = pr, names_from = Stat) %>%
  filter(M > 0) %>%
  mutate(Stage = factor(Stage, levels = c("inc", "sym", "aware", "care", "complete")),
         L = ifelse(Stage == "inc" | (Stage == "complete" & Index == "G"), NA, L),
         U = ifelse(Stage == "inc" | (Stage == "complete" & Index == "G"), NA, U))


dat_sex <- Cascade$Cascade_Sex %>%
  group_by(Country, ISO, Sex) %>%
  summarise(across(everything(), list(M = mean, 
                                      L = ~quantile(., 0.025),
                                      U = ~quantile(., 0.975)))) %>%
  pivot_longer(starts_with(c("C_", "G_")), values_to = "pr") %>%
  separate(name, c("Index", "Stage", "Stat"), "_") %>%
  select(-Key_M, -Key_L, -Key_U) %>%
  pivot_wider(c(Country, ISO, Sex, Index, Stage), values_from = pr, names_from = Stat) %>%
  filter(M > 0) %>%
  mutate(Stage = factor(Stage, levels = c("inc", "sym", "aware", "care", "complete")),
         L = ifelse(Stage == "inc" | (Stage == "complete" & Index == "G"), NA, L),
         U = ifelse(Stage == "inc" | (Stage == "complete" & Index == "G"), NA, U))


gs <- list()

gs$g_cohort <- Cohort$Cohort %>%
  filter(Sex == "Total") %>%
  ggplot() +
  geom_bar(aes(x = t, y = value, fill = name), stat = "identity", width = 0.02) +
  geom_bar(data = Cohort$End %>% filter(Sex == "Total"), aes(x = 2.2, y = value, fill = name), stat = "identity", width = 0.1) +
  scale_x_continuous("Time since TB-detectable, year", breaks = c(0, 0.5, 1, 1.5, 2, 2.2), 
                     labels = c(0, 0.5, 1, 1.5, 2, "End")) +
  scale_y_continuous("State Distribution, %", labels = scales::percent) + 
  scale_fill_manual("Stage",
                    values = colour_asc,
                    labels = labels_asc,
                    guide = guide_legend(reverse = TRUE)) +
  facet_wrap(.~Country, labeller = labeller(Country = countries_lab)) +
  theme(legend.position = "bottom")


gs$g_cascade <- ggplot(dat %>% filter(Index == "C")) +
  geom_bar(aes(x = Stage, y = M, fill = Stage), stat = "identity", width = 0.6) +
  geom_errorbar(aes(x = Stage, ymin = L, ymax = U), width = 0.3) +
  scale_y_continuous("Proportion, %", limits = c(0, 1), labels = scales::percent) +
  scale_x_discrete("", labels = NULL) + 
  scale_fill_manual("Endpoint",
                    values = colour_co_asc,
                    labels = labels_co_asc) +
  facet_wrap(.~Country, labeller = labeller(Country = countries_lab)) +
  theme(legend.position = "bottom")


gs$g_gaps <- ggplot(dat %>% filter(Index == "G")) +
  geom_bar(aes(x = Stage, y = M, fill = Stage), stat = "identity", width = 0.6) +
  geom_errorbar(aes(x = Stage, ymin = L, ymax = U), width = 0.3) +
  scale_y_continuous("Proportion, %", limits = c(0, 1), labels = scales::percent) +
  scale_x_discrete("", labels = NULL) + 
  scale_fill_manual("Endpoint",
                    values = colour_co_asc,
                    labels = labels_co_asc) +
  facet_wrap(.~Country, labeller = labeller(Country = countries_lab)) +
  theme(legend.position = "bottom")


gs$g_cohort_sex <- Cohort$Cohort %>%
  filter(Sex != "Total") %>%
  ggplot() +
  geom_bar(aes(x = t, y = value, fill = name), stat = "identity", width = 0.02) +
  geom_bar(data = Cohort$End %>% filter(Sex != "Total"), aes(x = 2.2, y = value, fill = name), stat = "identity", width = 0.1) +
  scale_x_continuous("Time since TB-detectable, year", breaks = c(0, 0.5, 1, 1.5, 2, 2.2), 
                     labels = c(0, 0.5, 1, 1.5, 2, "End")) +
  scale_y_continuous("State Distribution, %", labels = scales::percent) + 
  scale_fill_manual("Stage",
                    values = colour_asc,
                    labels = labels_asc,
                    guide = guide_legend(reverse = TRUE)) +
  facet_grid(Country~Sex, labeller = labeller(Country = countries_lab)) +
  theme(legend.position = "bottom")


gs$g_cascade_sex <- ggplot(dat_sex %>% filter(Index == "C")) +
  geom_bar(aes(x = Stage, y = M, fill = Stage), stat = "identity", width = 0.6) +
  geom_errorbar(aes(x = Stage, ymin = L, ymax = U), width = 0.3) +
  scale_y_continuous("Proportion, %", limits = c(0, 1), labels = scales::percent) +
  scale_x_discrete("", labels = NULL) + 
  scale_fill_manual("Endpoint",
                    values = colour_co_asc,
                    labels = labels_co_asc) +
  facet_grid(Country~Sex) +
  theme(legend.position = "bottom")

gs$g_gaps_sex <- ggplot(dat_sex %>% filter(Index == "G")) +
  geom_bar(aes(x = Stage, y = M, fill = Stage), stat = "identity", width = 0.6) +
  geom_errorbar(aes(x = Stage, ymin = L, ymax = U), width = 0.3) +
  scale_y_continuous("Proportion, %", limits = c(0, 1), labels = scales::percent) +
  scale_x_discrete("", labels = NULL) + 
  scale_fill_manual("Endpoint",
                    values = colour_co_asc,
                    labels = labels_co_asc) +
  facet_grid(Country~Sex) +
  theme(legend.position = "bottom")


ggsave(plot = gs$g_cascade, "docs/figs/Cascade" + ext, width = 9.5, height = 7.5)
ggsave(plot = gs$g_gaps, "docs/figs/Gaps" + ext, width = 9.5, height = 7.5)
ggsave(plot = gs$g_cohort, "docs/figs/Cohort" + ext, width = 9.5, height = 7.5)

ggsave(plot = gs$g_cascade_sex, "docs/figs/cascade/Cascade_Sex" + ext, width = 7.5, height = 15.5)
ggsave(plot = gs$g_gaps_sex, "docs/figs/cascade/Gaps_Sex" + ext, width = 7.5, height = 15.5)
ggsave(plot = gs$g_cohort_sex, "docs/figs/cascade/Cohort_Sex" + ext, width = 7.5, height = 15.5)



