library(ggplot2)

source("R/visualise.R")

theme_set(theme_bw() + theme(text = element_text(family = "sans")))

root <- "output/Malawi/Figures"
ext <- ".jpg"
width <- 6
height <- 7


## Fitted -----
load("output/Malawi/Fitted.rdata")
load("data/Input_Malawi.rdata")


gs_fitted <- visualise_fitted(fit_asym = fit_asym, fit_baseline = fit_baseline, inc = inc)

ggsave(plot = gs_fitted$error, filename = paste0(root, "Error", ext), width = width, height = height)
ggsave(plot = gs_fitted$delay, filename = paste0(root, "Delay", ext), width = width, height = height)
ggsave(plot = gs_fitted$duration, filename = paste0(root, "Duration", ext), width = width, height = 0.6 * height)
ggsave(plot = gs_fitted$inc, filename = paste0(root, "Incidence", ext), width = width, height = 0.6 * height)


## Periodic active case finding -----
library(gridExtra)
load(file = "output/Malawi/Intv.rdata")

gs_pacf_none <- visualise_pacf(intv_none, "None", c(4, 8, 12))
gs_pacf_sp <- visualise_pacf(intv_sp, "All smear-positive TB * 20%", c(4, 8, 12))
gs_pacf_sym <- visualise_pacf(intv_sym, "All symptomatic TB * 20%", c(4, 8, 12))
gs_pacf_all <- visualise_pacf(intv_all, "All active TB * 20%", c(4, 8, 12))

indices <- c("DelayRed", "Yield", "Smear", "InfRed", "Death")

for (i in 1:length(indices)) {
  key <- indices[i]
  gs <- grid.arrange(gs_pacf_sp[[i]], # gs_pacf_sn[[i]], 
                     gs_pacf_sym[[i]], gs_pacf_all[[i]],
                     ncol = 1)
  ggsave(plot = gs, filename = paste0(root, indices[i], ext), width = 1.5 * width, height = 2 * height)
}


for (i in 1:length(indices)) {
  key <- indices[i]

  ggsave(plot = gs_pacf_sp[[i]], filename = paste0(root, "ACF_Sp/", indices[i], ext), 
         width = width, height = 0.6 * height)
  # ggsave(plot = gs_pacf_sn[[i]], filename = paste0(root, "ACF_Sn/", indices[i], ext), 
  #       width = width, height = 0.6 * height)
  ggsave(plot = gs_pacf_sym[[i]], filename = paste0(root, "ACF_Sym/", indices[i], ext), 
         width = width, height = 0.6 * height)
  ggsave(plot = gs_pacf_all[[i]], filename = paste0(root, "ACF_All/", indices[i], ext), 
         width = width, height = 0.6 * height)
  
}


gs_bind <- vis_bind(intv_none, intv_sp, intv_sp, intv_sym, intv_all, 
                    sel_cycle = 4, sel_duration = c(4, 8, 12))

for (key in names(gs_bind)) {
  ggsave(plot = gs_bind[[key]], filename = paste0(root, "Cross", key, ext), 
         width = width, height = 0.55 * height)
}


## Individual effects -----
load("Output/Malawi/IndividualEffects.rdata")

gs_i_4_8 <- visualise_i_pacf_bi(i_pacf[[4]], i_pacf[[8]], dur1 = 4, dur2 = 8)
gs_i_6_12 <- visualise_i_pacf_bi(i_pacf[[6]], i_pacf[[12]], dur1 = 6, dur2 = 12)
gs_i <- visualise_i_pacf_12(i_pacf)
save(gs_i_4_8, gs_i_6_12, gs_i, file = "Output/Malawi/Figures/g_i.rdata")

## Population effects -----
load("Output/Malawi/PopulationEffects.rdata")

gs_p_4_05 <- visualise_p_pacf(p_pacf = p_pacf_05[[4]])
gs_p_8_05 <- visualise_p_pacf(p_pacf = p_pacf_05[[8]])
gs_p_4_20 <- visualise_p_pacf(p_pacf = p_pacf_20[[4]])
gs_p_8_20 <- visualise_p_pacf(p_pacf = p_pacf_20[[8]])

save(gs_p_4_05, gs_p_4_20, gs_p_8_05, gs_p_8_20, file = "Output/Malawi/Figures/g_p.rdata")  


