library(ggplot2)


source("R/visualise.R")

theme_set(theme_bw() + theme(text = element_text(family = "sans")))


ext <- ".jpg"
width <- 6
height <- 7


countries <- c("Malawi", "Kenya")

for (country in countries) {
  root <- paste0("Case", country)
  load(paste0("data/Input_", country, ".rdata"))
  load(paste0(root, "/Fitted.rdata"))
  
  gs_fitted <- visualise_fitted(fit_asym = fit_asym, fit_baseline = fit_baseline, inc = inc)

  save(gs_fitted, file = paste0(root, "/Vis_Fitted.rdata"))
  
  # ggsave(plot = gs_fitted$error, filename = paste0(root, "/Figures/Error", ext), width = width, height = height)
  ggsave(plot = gs_fitted$delay, filename = paste0(root, "/Figures/Delay", ext), width = width, height = height)
  ggsave(plot = gs_fitted$duration, filename = paste0(root, "/Figures/Duration", ext), width = width, height = 0.6 * height)
  ggsave(plot = gs_fitted$inc, filename = paste0(root, "/Figures/Incidence", ext), width = width, height = 0.6 * height)
} 