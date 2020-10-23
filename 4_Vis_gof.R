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

#### 

for (i in 1:length(countries)) {
  iso <- names(countries)[i]
  country <- countries[i]
  
  ### Load data ----
  load(paste0("data/Input_", iso, ".rdata"))
  load(paste0("out/ASC/Post_", iso, ".rdata"))
  
  
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
  
  
  dat_inc <- incidence %>% 
    left_join(rbind(
      notification %>%
        filter(Year == 2018) %>%
        group_by(Year) %>%
        summarise(Pop = sum(Pop), Sex = "Total"),
      notification %>%
        filter(Year == 2018) %>%
        select(Year, Pop, Sex)
    )) %>%
    mutate(value = m / Pop, vu = u / Pop, vl = l / Pop, Index = "Inc") %>%
    select(Year, Sex, Index, value, vl, vu)
  
  
  fore_inc <- local({
    inc_a <- as.table(extract(fitted_as_uni, c("inc_a"))$inc_a[1:200, , length(dat_as$Years)])
    dimnames(inc_a)[[1]] <- 1:200
    dimnames(inc_a)[[2]] <- c("Female", "Male")
    inc_a <- data.frame(inc_a)
    
    colnames(inc_a) <- c("Sim", "Sex", "value")
    inc_a$Year <- 2018
    inc_a$Index <- "asym"
    
    
    inc_s <- as.table(extract(fitted_as_uni, c("inc_s"))$inc_s[1:200, , length(dat_as$Years)])
    dimnames(inc_s)[[1]] <- 1:200
    dimnames(inc_s)[[2]] <- c("Female", "Male")
    inc_s <- data.frame(inc_s)
    
    colnames(inc_s) <- c("Sim", "Sex", "value")
    inc_s$Year <- 2018
    inc_s$Index <- "sym"
    
    rbind(inc_a, inc_s) %>%
      group_by(Sex, Year, Index) %>%
      summarise_at("value", list(m = median, vl = ~quantile(., 0.025), vu = ~quantile(., 0.975))) %>%
      rename(value = m)
  })
  
  
  fore_cnr <- local({
    noti <- as.table(extract(fitted_as_uni, c("noti"))$noti[1:200, , ])
    dimnames(noti)[[1]] <- 1:200
    dimnames(noti)[[2]] <- c("Female", "Male")
    dimnames(noti)[[3]] <- dat_as$Years
    noti <- data.frame(noti)
    
    colnames(noti) <- c("Sim", "Sex", "Year", "value")
    noti$Year <- as.numeric(as.character(noti$Year))
    
    noti
  })
  
  
  fore_prv <- local({
    ext <- extract(fitted_as_uni, c("prv", "pr_s", "pr_a"))
    
    prv_s <- ext$prv
    prv_a <- ext$prv
    
    for(i in 1:dim(ext$prv)[3]) {
      prv_s[, , i] <- prv_s[, , i] * ext$pr_s
      prv_a[, , i] <- prv_a[, , i] * ext$pr_a
    }
    
    prv_s <- as.table(prv_s[1:200, , ])
    dimnames(prv_s)[[1]] <- 1:200
    dimnames(prv_s)[[2]] <- c("Female", "Male")
    dimnames(prv_s)[[3]] <- dat_as$Years
    prv_s <- data.frame(prv_s)
    colnames(prv_s) <- c("Sim", "Sex", "Year", "value")
    prv_s$Year <- as.numeric(as.character(prv_s$Year))
    prv_s$Stage <- "S"
    
    prv_a <- as.table(prv_a[1:200, , ])
    dimnames(prv_a)[[1]] <- 1:200
    dimnames(prv_a)[[2]] <- c("Female", "Male")
    dimnames(prv_a)[[3]] <- dat_as$Years
    prv_a <- data.frame(prv_a)
    colnames(prv_a) <- c("Sim", "Sex", "Year", "value")
    prv_a$Year <- as.numeric(as.character(prv_a$Year))
    prv_a$Stage <- "A"
    
    rbind(prv_a, prv_s)
  })
  
  
  gs <- list()
  
  
  gs$g_Inc <- rbind(fore_inc, dat_inc %>% filter(Sex != "Total")) %>%
    mutate(Index = factor(Index, levels = c("asym", "sym", "Inc"))) %>%
    ggplot() +
    geom_pointrange(aes(x = Index, y = value, ymin = vl, ymax = vu, colour = Index), size = 1.5) +
    scale_x_discrete("Incidence definition", labels = c(asym = "A", sym = "S", Inc = "WHO")) + 
    scale_y_continuous("Incidence rate, per 100 000", labels = function(x) x * 1E5) +
    scale_colour_discrete("", labels = c(asym = "Est. Asymptomatic TB", sym = "Est. Symptomatic TB", Inc = "WHO estimates")) +
    facet_grid(.~Sex) +
    expand_limits(y = 0) +
    theme(legend.position = "none")
  
  
  gs$g_CNR <- ggplot(fore_cnr) +
    geom_line(aes(x = Year, y = value, group = Sim), alpha = 0.2, colour = "pink") +
    geom_point(data = dat_cnr %>% filter(Sex != "Total"), aes(x = Year, y = value)) +
    scale_y_continuous("Case notification rate, per 100 000", labels = function(x) x * 1E5) +
    scale_x_continuous("Year") + 
    facet_grid(.~Sex) +
    expand_limits(y = 0, x = 2010)
  
  
  gs$g_Prv <- ggplot(fore_prv) +
    geom_line(aes(x = Year, y = value, group = Sim), alpha = 0.2, colour = "pink") +
    geom_point(data = dat_prv %>% filter(Sex != "Total"), aes(x = Year, y = value)) +
    scale_y_continuous("Prevalence, per 100 000", labels = function(x) x * 1E5) +
    scale_x_continuous("Year") + 
    facet_grid(Stage~Sex, labeller = labeller(Stage = c(A = "Asymptomatic TB", S = "Symptomatic TB"))) +
    expand_limits(y = 0, x = 2010)

  
  save(gs, file = paste0("out/g_Fitted_", iso, ".rdata"))
  
  ggsave(plot = gs$g_Inc, paste0("docs/figs/fit/Inc_", iso, ext), width = 5.5, height = 3.5)
  ggsave(plot = gs$g_CNR, paste0("docs/figs/fit/CNR_", iso, ext), width = 6.5, height = 3.5)
  ggsave(plot = gs$g_Prv, paste0("docs/figs/fit/Prv_", iso, ext), width = 6.5, height = 5.5)
}







