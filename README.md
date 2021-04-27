# AsymTB
Bayesian empirical models for considering asymptomatic TB


## Theoretical models



## Scripts for model fitting

- [2_Fit_ANP.R](scripts/2_Fit_ANP.R): Asymptomatic->Smr-/Smr+ model
- [2_Fit_ASC.R](scripts/2_Fit_ASC.R): Asymptomatic->Symptomatic with/without care-seeking model
- [2_Fit_Full_BLT.R](scripts/2_Fit_Full_BLT.R): For Blantyre data with strata of age, sex, HIV
- [2_Fit_Full_KEN.R](scripts/2_Fit_Full_KEN.R): For Kenya data with strata of age, sex, HIV



## Model inference

- [3_Inference.R](scripts/3_Inference.R): inference with the posterior distribution


## Sensitivity analyses

- [2_Fit_Sens_Sym.R](scripts/2_Fit_Sens_Sym.R): definitions of TB-like symptoms
- [2_Fit_Sens_Mor.R](scripts/2_Fit_Sens_Mor.R): assumptions for TB-related deaths
- [2_Fit_Sens_HIV.R](scripts/2_Fit_Sens_HIV.R): assumptions for TB-related deaths of PLHIV


## Table outputs

- [4_Tab_basic.R](scripts/4_Tab_basic.R) Data summary
- [4_Tab_posterior.R](scripts/4_Tab_posterior.R) Posterior distributions by country
- [4_Tab_HIV.R](scripts/4_Tab_HIV.R) Summary of the HIV inference 
- [4_Tab_covariates.R](scripts/4_Tab_covariates.R) Regression analysis for the asymptomatic duration



## Data visualisation

- [4_Vis_gof.R](scripts/4_Vis_gof.R) Goodness of fit figures
- [4_Vis_duration.R](scripts/4_Vis_duration.R) Time to event distributions
- [4_Vis_conversion.R](scripts/4_Vis_conversion.R) Smear-type conversion
- [4_Vis_burden.R](scripts/4_Vis_burden.R) Figures for burden estimates
- [4_Vis_cascade.R](scripts/4_Vis_cascade.R) Health-care cascade



