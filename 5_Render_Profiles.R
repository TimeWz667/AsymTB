
library(ggplot2)
library(tidyverse)
library(tidybayes)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))


### Country list -----
countries <- c(
  KHM = "Cambodia",
  KEN = "Kenya",
  LAO = "Lao People's Democratic Republic", 
  MWI = "Malawi", 
  PAK = "Pakistan", 
  PHL = "Philippines", 
  TZA = "United Republic of Tanzania", 
  UGA = "Uganda", #
  VNM = "Viet Nam", 
  ZMB = "Zambia"
)


for (i in 1:length(countries)) {
  country <- countries[i]
  
  iso <- names(countries)[i]
  
  gss <- list()
  
  load(file = paste0("out/g_Cascade_", iso, ".rdata"))
  gss <- c(gs)
  load(file = paste0("out/g_Duration_", iso, ".rdata"))
  gss <- c(gss, gs)
  load(file = paste0("out/g_Fitted_", iso, ".rdata"))
  gss <- c(gss, gs)
  
  rmarkdown::render("docs/country/temp_profile.Rmd",
                    output_format = rmarkdown::html_document(toc = T, toc_float = T, smooth_scroll = T, toc_depth = 2),
                    output_file = paste0("Profile_", iso, ".html"),
                    envir = list(gs = gss, country = country, iso = iso))
}



