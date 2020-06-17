library(ggplot2)
library(tidyverse)

source("R/vis_p_pacf.R")

theme_set(theme_bw() + theme(text = element_text(family = "sans")))



## Visualise _____
countries <- c("Malawi", "Kenya")

for(country in countries) {
  gs_p_4_05 <- visualise_p_pacf(p_pacf = p_pacf_05[[4]])
  gs_p_8_05 <- visualise_p_pacf(p_pacf = p_pacf_05[[8]])
  gs_p_4_20 <- visualise_p_pacf(p_pacf = p_pacf_20[[4]])
  gs_p_8_20 <- visualise_p_pacf(p_pacf = p_pacf_20[[8]])
  
  save(gs_p_4_05, gs_p_4_20, gs_p_8_05, gs_p_8_20, file = paste0("Case", country, "/Vis_p_pacf.rdata"))
}




dur_asym <- 8
country = "Malawi"
load(file = paste0("Case", country, "/PACF_population.rdata"))

p_pacf_1 <- p_pacf_20[[dur_asym]]


country = "Kenya"
load(file = paste0("Case", country, "/PACF_population.rdata"))

p_pacf_2 <- p_pacf_20[[dur_asym]] 

lab_1 <- "Malawi"
lab_2 <- "Kenya"


root <- "output"
ext <- ".pdf"
width <- 8
height <- 7


gs_bi <- visualise_p_pacf(p_pacf_1, p_pacf_2, lab_1, lab_2)


ggsave(plot = gs_bi$g_occur, filename = paste0(root, "/Forecasts", ext), width = width*1.4, height = height)

ggsave(plot = gs_bi$g_avert, filename = paste0(root, "/Averted", ext), width = width, height = height)


