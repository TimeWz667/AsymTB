library(ggplot2)
library(tidyverse)
library(ggpubr)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))

ext <- glue::as_glue(".png")


n_iter <- 1E4

prior <- tibble(
  prv0 = runif(n_iter, 0, 1),
  r_det = 1 / rgamma(n_iter, 1, 1),
  r_sym = 1 / rgamma(n_iter, 1, 1),
  r_sc = runif(n_iter, 0.1, 0.3)
) %>%
  mutate(
    ra = r_sc + 0.002,
    rs = r_sc + 0.002,
    
    a0 = rs + r_det,
    s0 = r_sym,

    pra = a0 / (a0 + s0),
    prs = s0 / (a0 + s0),
    prva = prv0 * pra,
    prvs = prv0 * prs,
    nr = prvs * r_det,
    inc = (ra * pra + (rs + r_det) * prs) * prv0,
    cdr = nr / inc,
    dura = 1 / (ra + r_sym),
    durs = 1 / (rs + r_det),
    GapA = r_sym * dura,
    GapS = r_det * durs
  )



gs <- list()

gs$g_cdr <- prior %>%
  ggplot() +  
  geom_density(aes(x = cdr), fill = "grey") +
  scale_x_continuous("Case-detection ratio, %", limits = c(0, 1), labels = scales::percent) +
  scale_y_continuous("p.d.f.")

gs$g_dura <- prior %>%
  ggplot() +  
  geom_density(aes(x = dura), fill = "grey") +
  scale_x_continuous("Asymptomatic period, months", labels = function(x) x * 12, limits = c(0, 5)) +
  scale_y_continuous("p.d.f.")

gs$g_durs <- prior %>%
  ggplot() +  
  geom_density(aes(x = durs), fill = "grey") +
  scale_x_continuous("Symptomatic period, months", labels = function(x) x * 12, limits = c(0, 5)) +
  scale_y_continuous("p.d.f.")

gs$g_gapa <- prior %>%
  ggplot() +  
  geom_density(aes(x = GapA), fill = "grey") +
  scale_x_continuous("Progress to symptomatic stage, %", limits = c(0, 1), labels = scales::percent) +
  scale_y_continuous("p.d.f.")

gs$g_gaps <- prior %>%
  ggplot() +  
  geom_density(aes(x = GapS), fill = "grey") +
  scale_x_continuous("Progress to case-detection, %", limits = c(0, 1), labels = scales::percent) +
  scale_y_continuous("p.d.f.")


gs$g_bind <- ggarrange(
  gs$g_dura, gs$g_durs,
  gs$g_gapa, gs$g_gaps,
  gs$g_cdr,
  nrow = 5
)


ggsave(plot = gs$g_bind, "docs/figs/Prior" + ext, width = 6, height = 7)

