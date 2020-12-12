
# Prior to model fitting
source("1_collect_data.R")
source("1_test_simulation_model.R")
source("1_test_theoretical_model.R")


# Model fitting
source("2_Fit_ASC.R")
source("2_Fit_ASC_Sp.R")
source("2_Fit_ANP.R")
source("2_Fit_Full_BLT.R")
source("2_Fit_Full_KEN.R")

source("4_Vis_gof.R")
source("4_Tab_posterior.R")


# Inference
source("3_Inference.R")


# Visualisation / Table making
source("4_Tab_basic.R")
source("4_Tab_conversion.R")

source("4_Vis_duration.R")
source("4_Vis_cascade.R")
source("4_Vis_conversion.R")
source("4_Vis_burden.R")


# Sensitivity analysis
source("2_Fit_Sens_Mor.R")
source("2_Fit_Sens_Sym.R")
source("2_Fit_Sens_HIV.R")
source("4_Tab_Mor.R")
source("4_Tab_Sym.R")
source("4_Tab_HIV.R")


# Making profiles
source("5_Render_Profiles.R")



