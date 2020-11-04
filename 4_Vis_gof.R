library(ggplot2)
library(tidyverse)
library(rstan)

#### Countries
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

ext <- ".png"

n_sel <- 300

#### 

for (i in 1:length(countries)) {
  iso <- glue::as_glue(names(countries)[i])
  country <- countries[i]
  
  ### Load data ----
  load("data/Input_" + iso + ".rdata")
  load("out/ASC/Post_" + iso + ".rdata")
  
  
  ### Transform data ----
  dat_prv <- rbind(
    prevalence %>%
      mutate(Prv_A = Asym / N, Prv_S = Sym / N) %>%
      pivot_longer(starts_with("Prv_")) %>%
      separate(name, c("Index", "Stage"), "_") %>%
      select(Year, Sex, Index, Stage, value),
    prevalence %>%
      group_by(Year) %>%
      summarise(Prv_A = sum(Asym) / sum(N), Prv_S = sum(Sym) / sum(N), Sex = "Total") %>%
      pivot_longer(starts_with("Prv_")) %>%
      separate(name, c("Index", "Stage"), "_") %>%
      select(Year, Sex, Index, Stage, value)
  )
  
  
  dat_cnr <- rbind(
    notification %>%
      mutate(value = n_all / Pop, Index = "CNR") %>%
      select(Year, Sex, Index, value),
    notification %>%
      group_by(Year) %>%
      summarise(value = sum(n_all) / sum(Pop), Sex = "Total", Index = "CNR") %>%
      select(Year, Sex, Index, value)
  )
  
  
  fore <- local({
    tt <- 2010:2019 - dat_as$YearSurveyed
    
    ext <- extract(fitted_as_uni, c("prv_a", "prv_s", "adr", "r_det"))
  
    data.table::rbindlist(lapply(2010:2019, function(t) {
      prv_s <- exp(- ext$adr * (t - dat_as$YearSurveyed)) * ext$prv_s
      prv_a <- exp(- ext$adr * (t - dat_as$YearSurveyed)) * ext$prv_a
      cnr <- prv_s * ext$r_det
      
      prv_s <- as.table(prv_s[1:n_sel, ])
      prv_a <- as.table(prv_a[1:n_sel, ])
      cnr <- as.table(cnr[1:n_sel, ])
      dimnames(prv_s)[[1]] <- dimnames(prv_a)[[1]] <- dimnames(cnr)[[1]] <- 1:n_sel
      dimnames(prv_s)[[2]] <- dimnames(prv_a)[[2]] <- dimnames(cnr)[[2]] <- c("Female", "Male")
      
      
      prv_s <- data.frame(prv_s)
      colnames(prv_s) <- c("Sim", "Sex", "value")
      prv_s$Year = t
      prv_s$Stage <- "S"
      
      prv_a <- data.frame(prv_a)
      colnames(prv_a) <- c("Sim", "Sex", "value")
      prv_a$Year = t
      prv_a$Stage <- "A"
      
      
      cnr <- data.frame(cnr)
      colnames(cnr) <- c("Sim", "Sex", "value")
      cnr$Year = t
      cnr$Stage <- "N"
      
      rbind(prv_a, prv_s, cnr)
    }))
    
  })  
  
  
  gs <- list()
  
  gs$g_CNR <- ggplot(fore %>% filter(Stage == "N")) +
    geom_line(aes(x = Year, y = value, group = Sim), alpha = 0.2, colour = "pink") +
    geom_point(data = dat_cnr %>% filter(Sex != "Total"), aes(x = Year, y = value)) +
    scale_y_continuous("Case notification rate, per 100 000", labels = function(x) x * 1E5) +
    scale_x_continuous("Year") + 
    facet_grid(.~Sex) +
    expand_limits(y = 0, x = 2010)
  
  
  gs$g_Prv <- ggplot(fore %>% filter(Stage != "N")) +
    geom_line(aes(x = Year, y = value, group = Sim), alpha = 0.2, colour = "pink") +
    geom_point(data = dat_prv %>% filter(Sex != "Total"), aes(x = Year, y = value)) +
    scale_y_continuous("Prevalence, per 100 000", labels = function(x) x * 1E5) +
    scale_x_continuous("Year") + 
    facet_grid(Stage~Sex, labeller = labeller(Stage = c(A = "Asymptomatic TB", S = "Symptomatic TB"))) +
    expand_limits(y = 0, x = 2010)

  
  save(gs, file = "out/g_Fitted_" + iso + ".rdata")
  
  ggsave(plot = gs$g_CNR, "docs/figs/fit/CNR_" + iso + ext, width = 6.5, height = 3.5)
  ggsave(plot = gs$g_Prv, "docs/figs/fit/Prv_" + iso + ext, width = 6.5, height = 5.5)
}
