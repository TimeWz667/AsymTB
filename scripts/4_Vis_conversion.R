library(tidyverse)
library(ggplot2)
library(latex2exp)
library(ggpubr)



theme_set(theme_bw() + theme(text = element_text(family = "sans")))

ext <- glue::as_glue(".png")

source("data/country_list.R")


load("out/Conversion_All.rdata")

dat_cor <- bind_rows(lapply(names(countries), function(iso) {
  reg <- Conversion$Reg[[iso]]
  reg <- coef(reg)
  pars <- Conversion$Pars %>% filter(ISO == iso)
  names(reg) <- c("intrp", "beta")
  reg <- as.list(reg)
  reg$ISO <- iso
  reg$Country <- countries_lab[countries[iso]]
  
  reg
})) %>%
  left_join(Conversion$Wts) 

line_cor <- dat_cor %>% 
  filter(!(ISO %in% c("MWI", "VNM"))) %>%
  summarise(intrp = weighted.mean(intrp, NotiN), beta = weighted.mean(beta, NotiN))





gs <- list()


gs$g_cor <- dat_cor %>%
  ggplot() +
  geom_point(aes(x = beta, y = intrp, size = NotiN), alpha = 0.6) +
  geom_text(aes(x = beta, y = intrp, label = ISO), hjust = 0.4, vjust = 1.8, size = rel(3)) +
  scale_y_continuous("Smear-positive at symtomatic onset, % \n without smear conversion", 
                     limits = c(0, 1), labels = scales::percent) +
  scale_x_continuous(TeX("Slope of (conversion rate, percentage smear-positive)"), 
                     limits = c(-1, 0)) +
  scale_size_continuous("Notified cases, 2019", labels = function(x) paste0(x * 1E-6, " million")) +
  guides(colour = "none") +
  #expand_limits(x = c(-1, 0)) +
  #theme(legend.position = c(0.05, 0.1), legend.direction = "horizontal", 
  #      legend.justification = "left", legend.background = element_rect(fill = "grey90"))
  theme(legend.position = "bottom")


xlab <- latex2exp::TeX("Smear-positive at symtomatic onset ($\\pi$), %")
ylab <- latex2exp::TeX("Conversion rate ($\\tau$), per year")


gs$g_tr_sp <- Conversion$Pars %>%
  mutate(Country = countries_lab[countries[ISO]]) %>%
  ggplot() +
  geom_point(data = Conversion$Pars %>%
               mutate(Country = countries_lab[countries[ISO]]), aes(x = r_tr, y = p_sp), alpha = 0.1) + 
  scale_y_continuous(xlab, 
                     limits = c(0, 1), labels = scales::percent) +
  scale_x_continuous(ylab) +
  facet_wrap(.~Country) +
  scale_fill_viridis_c(alpha = 0.7, na.value = "white")


gs$g_den <- Conversion$Dens %>% 
  left_join(Conversion$Wts) %>%
  filter(!(ISO %in% c("MWI", "VNM"))) %>%
  mutate(Country = countries_lab[countries[ISO]], wt = NotiN) %>%
  group_by(x, y) %>%
  summarise(z = weighted.mean(z, wt)) %>%
  filter(z > 0.5) %>%
  ggplot() +
  geom_tile(aes(x, y, fill = z)) +
  annotate("segment", x = 0, xend = - line_cor$intrp / line_cor$beta, y = line_cor$intrp, yend = 0, size = 1.5) + 
  annotate("text", x = - line_cor$intrp / line_cor$beta / 2, y = line_cor$intrp,
           label = latex2exp::TeX(sprintf("$\\pi = %.3f - %.3f \\tau$", line_cor$intrp, abs(line_cor$beta))), parse = TRUE,
           hjust = 0) +
  scale_y_continuous(xlab, 
                     limits = c(0, 1), labels = scales::percent) +
  scale_x_continuous(ylab) +
  scale_fill_viridis_c("", alpha = 0.7, na.value = "white", direction = -1, breaks = c(0.9, 2.5), labels = c("Low", "High")) +
  #theme(legend.position = c(0.95, 0.85), legend.direction = "horizontal", 
  #      legend.justification = "right", legend.background = element_rect(fill = "grey90"))
  theme(legend.position = "bottom")


gs$g_bind <- ggarrange(
  gs$g_tr_sp + labs(title = "A. Joint probability density of percentage smear-positive and conversion rate (smear- to smear+)"),
  ggarrange(
    gs$g_cor + labs(title = "B. Correlation"),
    gs$g_den + labs(title = "C. Pooled joint probability density"),
    widths = c(2, 1.3)
  ),
  nrow = 2, heights = c(1.1, 1.3)
)


ggsave(plot = gs$g_cor, "docs/figs/conversion/Conversion_cor" + ext, width = 6.5, height = 4.7)
ggsave(plot = gs$g_tr_sp, "docs/figs/conversion/Conversion_cross" + ext, width = 7.5, height = 7.5)
ggsave(plot = gs$g_den, "docs/figs/conversion/Conversion_density" + ext, width = 6.5, height = 4.7)
ggsave(plot = gs$g_bind, "docs/figs/Conversion" + ext, width = 10, height = 11)

