library(ggplot2)


source("R/vis_fitted.R")

theme_set(theme_bw() + theme(text = element_text(family = "sans")))


ext <- ".pdf"
width <- 6
height <- 7


countries <- c("Malawi", "Kenya")
gf <- list()

for (country in countries) {
  root <- paste0("Case", country)
  load(paste0("data/Input_", country, ".rdata"))
  load(paste0(root, "/Fitted.rdata"))
  
  gs_fitted <- visualise_fitted(fit_asym = fit_asym, fit_baseline = fit_baseline, inc = inc)
  gf[[country]] <- gs_fitted
  save(gs_fitted, file = paste0(root, "/Vis_Fitted.rdata"))
  
  # ggsave(plot = gs_fitted$error, filename = paste0(root, "/Figures/Error", ext), width = width, height = height)
  ggsave(plot = gs_fitted$delay, filename = paste0(root, "/Figures/Delay", ext), width = width, height = height)
  ggsave(plot = gs_fitted$duration, filename = paste0(root, "/Figures/Duration", ext), width = width, height = 0.6 * height)
  ggsave(plot = gs_fitted$inc, filename = paste0(root, "/Figures/Incidence", ext), width = width, height = 0.6 * height)
}


ext <- ".pdf"
width <- 8
height <- 8

gs_fitted_bind <- visualise_bind(gf)
ggsave(plot = gs_fitted_bind$duration, filename = paste0("output/Duration", ext), width = width, height = height)
ggsave(plot = gs_fitted_bind$inc, filename = paste0("output/Incidence", ext), width = width, height = height)
