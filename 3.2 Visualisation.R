library(ggplot2)
library(tidyverse)

source("R/vis_i_pacf.R")

theme_set(theme_bw() + theme(text = element_text(family = "sans")))



## Visualise _____
countries <- c("Malawi", "Kenya")

for (country in countries) {
  load(file = paste0("Case", country, "/PACF_individual.rdata"))
  
  gs_i_4_8 <- visualise_i_pacf_bi(i_pacf[[4]], i_pacf[[8]], dur1 = 4, dur2 = 8)
  gs_i_6_12 <- visualise_i_pacf_bi(i_pacf[[6]], i_pacf[[12]], dur1 = 6, dur2 = 12)
  gs_i <- visualise_i_pacf_12(i_pacf)
  
  save(gs_i_4_8, gs_i_6_12, gs_i, file = paste0("Case", country, "/Vis_i_pacf.rdata"))
}



# load(paste0("data/Input_", country, ".rdata"))
country <- "Malawi"
load(file = paste0("Case", country, "/PACF_individual.rdata"))
idx_m <- data.table::rbindlist(lapply(1:12, function(dur_asym) {
  extract_idx(i_pacf[[dur_asym]], dur_asym)
})) %>% filter(name %in% c("Death", "Dur", "Inf")) %>%
  mutate(name = factor(name, levels = c("Death", "Dur", "Inf")), country = country)

country <- "Kenya"
load(file = paste0("Case", country, "/PACF_individual.rdata"))
idx_y <- data.table::rbindlist(lapply(1:12, function(dur_asym) {
  extract_idx(i_pacf[[dur_asym]], dur_asym)
})) %>% filter(name %in% c("Death", "Dur", "Inf")) %>%
  mutate(name = factor(name, levels = c("Death", "Dur", "Inf")), country = country)

sel_dur = c(4, 8)


idx <- rbind(idx_y, idx_m) %>%
  mutate(m = ifelse(name == "Death", 100*m, m),
         l = ifelse(name == "Death", 100*l, l),
         u = ifelse(name == "Death", 100*u, u))

g_outcome <- idx %>% filter(DurAsym %in% sel_dur) %>%
  ggplot() +
  geom_pointrange(aes(x = ACF, y = m, ymin = l, ymax = u, colour = Sex), position = position_dodge(0.2)) +
  facet_grid(name~country+DurAsym, scales = "free_y", 
             labeller = labeller(DurAsym = label_both, name = c(Dur="Duration", "Inf"="Transmissibility", Death="Untreated death"))) +
  scale_x_discrete("Inclusion criteria for ACF") +
  scale_y_continuous("per contact-year / year / percentage") +
  expand_limits(y = 0) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))


## Output -----

root <- "output"
ext <- ".pdf"
width <- 8
height <- 7

ggsave(plot = g_outcome, filename = paste0(root, "/PACF_I", ext), width = width, height = height)
