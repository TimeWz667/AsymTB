library(ggplot2)
library(tidyverse)
library(tidybayes)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))


### Country list -----
source("data/country_list.R")



for (i in 1:length(countries)) {
  country <- countries[i]
  
  iso <- names(countries)[i]
  
  gss <- list()
  
  load(file = paste0("out/g_Cascade_", iso, ".rdata"))
  gss <- c(gs)
  load(file = paste0("out/g_TTE_", iso, ".rdata"))
  gss <- c(gss, gs)
  load(file = paste0("out/g_Fitted_", iso, ".rdata"))
  gss <- c(gss, gs)
  load(file = paste0("out/g_Incidence_", iso, ".rdata"))
  gss$g_Inc <- g_Inc
  
  rmarkdown::render("docs/country/temp_profile.Rmd",
                    output_format = rmarkdown::html_document(toc = T, toc_float = T, smooth_scroll = T, toc_depth = 2),
                    output_file = paste0("Profile_", iso, ".html"),
                    envir = list(gs = gss, country = country, iso = iso))
  
  if (iso %in% c("KEN")) {
    rmarkdown::render("docs/country/temp_profile.Rmd",
                    output_format = rmarkdown::pdf_document(toc = T, toc_depth = 2),
                    output_file = paste0("Profile_", iso, ".pdf"),
                    envir = list(gs = gss, country = country, iso = iso))
  }
  
}



