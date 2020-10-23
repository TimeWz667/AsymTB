rm(list = ls())

library(ggplot2)
library(tidyverse)

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
  UGA = "Uganda", 
  ZMB = "Zambia"
)

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
           sym = "#fdd0a2")
labels_co_asc = c(complete = "Treatment complete",
           care = "Treatment start",
           aware = "Care-seeking start",
           sym = "Symptom onset")

colour_co_as = c(complete = "#bcbddc",
                  care = "#4BB35D",
                  sym = "#fdd0a2")
labels_co_as = c(complete = "Treatment complete",
                  care = "Treatment start",
                  sym = "Symptom onset")



for (i in 1:length(countries)) {
  iso <- names(countries)[i]
  country <- countries[i]
  
  ### Load data ----
  load(paste0("out/Cascade_", iso, ".rdata"))

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
  
  
  gs$g_Cohort_Sex <- Cohort$TS %>% filter(Sex != "Total") %>%
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
  
  gs$g_Cohort_Total <-  Cohort$TS %>% filter(Sex == "Total") %>%
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
  
  gs$g_Cascade <- ggplot(Cascade %>% filter(Index == "C")) +
    geom_bar(aes(x = Stage, y = M, fill = Stage), stat = "identity", width = 0.6) +
    geom_errorbar(aes(x = Stage, ymin = L, ymax = U), width = 0.3) +
    scale_y_continuous("Proportion (%)", limits = c(0, 1), labels = scales::percent) +
    scale_x_discrete("Stage", 
                     labels = lab_co) + 
    scale_fill_manual("Arrivals",
                      values = col_co,
                      labels = lab_co) +
    theme(legend.position = "none") +
    labs(title = country)

  save(gs, file = paste0("out/g_Cascade_", iso, ".rdata"))
  
  ggsave(plot = gs$g_Cohort_Sex, paste0("docs/figs/Cascade/CohortSex_", iso, ext), width = 7.5, height = 7.5)
  ggsave(plot = gs$g_Cohort_Total, paste0("docs/figs/Cascade/Cohort_", iso, ext), width = 7.5, height = 4.5)
  ggsave(plot = gs$g_Cascade, paste0("docs/figs/Cascade/Cascade_", iso, ext), width = 6.5, height = 4.5)
}




