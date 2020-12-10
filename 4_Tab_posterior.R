library(tidyverse)
library(rstan)


source("data/country_list.R")



variable_as_map <- c(
  "r_sym" = "rate to symptomatic onset",
  "r_det" = "rate to case-detection",
  "r_sc" = "self-cure rate",
  "adr" = "annual decline rate, %"
)

variable_asc_map <- c(
  "r_sym" = "rate to symptomatic onset",
  "r_aware" = "rate to initial care-seeking",
  "r_det" = "rate to case-detection",
  "r_sc" = "self-cure rate",
  "adr" = "annual decline rate, %"
)

variable_anp_map <- c(
  "r_sym" = "rate to symptomatic onset",
  "r_tr" = "smear conversion rate",
  "p_sp" = "smear+ at symptom onset",
  "r_det_sn" = "rate to case-detection, smear-",
  "r_det_sp" = "rate to case-detection, smear+",
  "r_sc" = "rate to self-cure",
  "adr" = "annual decline rate, %"
)


for (i in 1:length(countries)) {
  iso <- glue::as_glue(names(countries)[i])
  country <- countries[i]
  
  load("out/ASC/post_" + iso + ".rdata")
  
  if (country %in% countries_cs) {
    tab1 <- summary(fitted_asc_uni, pars = c("r_sym", "r_det", "r_aware", "adr"))$summary
    tab1[c("adr[1]", "adr[2]"), c(1, 4:8)] <- tab1[c("adr[1]", "adr[2]"), c(1, 4:8)] * 100
    tab1[c("adr[1]", "adr[2]"), 2:3] <- tab1[c("adr[1]", "adr[2]"), 2:3] * 100
    
    tab1 <- bind_cols(name = rownames(tab1), tab1[, c(1, 3, 4, 8:10)]) %>%
      tidyr::extract(name, c("Variable", "Sex"), "(\\w+)\\[(\\d)\\]")%>%
      mutate(Sex = ifelse(Sex == 1, "Female", "Male"), 
             Variable = variable_asc_map[Variable],
             Variable = factor(Variable, levels = variable_asc_map))
    
    tab2 <- summary(fitted_asc_uni, pars = c("r_sc"))$summary
    tab2 <- cbind(tibble(Variable = rownames(tab2)), t(tab2[, c(1, 3, 4, 8:10)])) %>%
      mutate(Sex = "Shared", 
             Variable = variable_asc_map[Variable],
             Variable = factor(Variable, levels = variable_asc_map))
    
    tab <- bind_rows(tab1, tab2) %>% 
      rename(Mean = mean, SD = sd, ESS = n_eff, Rhat = Rhat) %>%
      mutate(`95% CrI` = sprintf("(%.2f, %.2f)", `2.5%`, `97.5%`)) %>%
      mutate(Sex = factor(Sex, levels = c("Female", "Male", "Shared"))) %>%
      select(- c(`2.5%`, `97.5%`)) %>%
      relocate(Variable, Sex, Mean, SD, `95% CrI`, ESS, Rhat) %>%
      arrange(Sex, Variable)
    
    tab.tex <- knitr::kable(tab, "latex", digits = c(NA, NA, 2, 2, NA, 0, 3), 
                            booktabs = T, align = "llccccc", 
                            linesep = c(rep("", 4), "\\addlinespace"),
                            caption = "Posterior distributions of parameters, Model B")
  } else {
    tab1 <- summary(fitted_as_uni, pars = c("r_sym", "r_det", "adr"))$summary
    tab1[c("adr[1]", "adr[2]"), c(1, 4:8)] <- tab1[c("adr[1]", "adr[2]"), c(1, 4:8)] * 100
    tab1[c("adr[1]", "adr[2]"), 2:3] <- tab1[c("adr[1]", "adr[2]"), 2:3] * 100
    
    tab1 <- bind_cols(name = rownames(tab1), tab1[, c(1, 3, 4, 8:10)]) %>%
      tidyr::extract(name, c("Variable", "Sex"), "(\\w+)\\[(\\d)\\]")%>%
      mutate(Sex = ifelse(Sex == 1, "Female", "Male"), 
             Variable = variable_as_map[Variable],
             Variable = factor(Variable, levels = variable_as_map))
    
    tab2 <- summary(fitted_as_uni, pars = c("r_sc"))$summary
    tab2 <- cbind(tibble(Variable = rownames(tab2)), t(tab2[, c(1, 3, 4, 8:10)])) %>%
      mutate(Sex = "Shared", 
             Variable = variable_as_map[Variable],
             Variable = factor(Variable, levels = variable_as_map))
    
    tab <- bind_rows(tab1, tab2) %>% 
      rename(Mean = mean, SD = sd, ESS = n_eff, Rhat = Rhat) %>%
      mutate(`95% CrI` = sprintf("(%.2f, %.2f)", `2.5%`, `97.5%`)) %>%
      mutate(Sex = factor(Sex, levels = c("Female", "Male", "Shared"))) %>%
      select(- c(`2.5%`, `97.5%`)) %>%
      relocate(Variable, Sex, Mean, SD, `95% CrI`, ESS, Rhat) %>%
      arrange(Sex, Variable)
    
    tab.tex <- knitr::kable(tab, "latex", digits = c(NA, NA, 2, 2, NA, 0, 3), 
                            booktabs = T, align = "llccccc", 
                            linesep = c(rep("", 3), "\\addlinespace"),
                            caption = "Posterior distributions of parameters, Model A")
  }
  tab.tex <- gsub("Rhat", "\\\\^\\{R\\}", tab.tex)
  tab.tex <- gsub("begin\\{table\\}", "begin\\{table\\}\\[h\\]", tab.tex)
  
  write(tab.tex, "docs/tabs/posterior/ASC_" + iso + ".tex")
  
  
  load("out/ANP/post_" + iso + ".rdata")
  
  tab1 <- summary(fitted_anp_uni, pars = c("r_sym", "r_det_sn", "r_det_sp", "adr"))$summary
  tab1[c("adr[1]", "adr[2]"), c(1, 4:8)] <- tab1[c("adr[1]", "adr[2]"), c(1, 4:8)] * 100
  tab1[c("adr[1]", "adr[2]"), 2:3] <- tab1[c("adr[1]", "adr[2]"), 2:3] * 100
  
  tab1 <- bind_cols(name = rownames(tab1), tab1[, c(1, 3, 4, 8:10)]) %>%
    tidyr::extract(name, c("Variable", "Sex"), "(\\w+)\\[(\\d)\\]")%>%
    mutate(Sex = ifelse(Sex == 1, "Female", "Male"), 
           Variable = variable_anp_map[Variable],
           Variable = factor(Variable, levels = variable_anp_map))
  
  
  tab2 <- summary(fitted_anp_uni, pars = c("r_sc", "r_tr", "p_sp"))$summary
  tab2 <- bind_cols(Variable = rownames(tab2), tab2[, c(1, 3, 4, 8:10)]) %>%
    mutate(Sex = "Shared", 
           Variable = variable_anp_map[Variable],
           Variable = factor(Variable, levels = variable_anp_map))
  
  tab <- bind_rows(tab1, tab2) %>% 
    rename(Mean = mean, SD = sd, ESS = n_eff, Rhat = Rhat) %>%
    mutate(`95% CrI` = sprintf("(%.2f, %.2f)", `2.5%`, `97.5%`)) %>%
    mutate(Sex = factor(Sex, levels = c("Female", "Male", "Shared"))) %>%
    select(- c(`2.5%`, `97.5%`)) %>%
    relocate(Variable, Sex, Mean, SD, `95% CrI`, ESS, Rhat) %>%
    arrange(Sex, Variable)
  
  tab.tex <- knitr::kable(tab, "latex", digits = c(NA, NA, 2, 2, NA, 0, 3), 
                          booktabs = T, align = "llccccc", 
                          linesep = c(rep("", 3), "\\addlinespace"),
                          caption = "Posterior distributions of parameters, Model C")
  
  tab.tex <- gsub("Rhat", "\\\\^\\{R\\}", tab.tex)
  tab.tex <- gsub("begin\\{table\\}", "begin\\{table\\}\\[h\\]", tab.tex)
  
  write(tab.tex, "docs/tabs/posterior/ANP_" + iso + ".tex")
  
}




