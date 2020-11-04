library(tidyverse)
library(rstan)

### Country list -----
countries <- c(
  KHM = "Cambodia",
  KEN = "Kenya",
  LAO = "Lao People's Democratic Republic", 
  MWI = "Malawi", 
  PAK = "Pakistan", 
  PHL = "Philippines", 
  TZA = "United Republic of Tanzania", 
  UGA = "Uganda",
  VNM = "Viet Nam", 
  ZMB = "Zambia"
)



incs <- list()


for (i in 1:length(countries)) {
  iso <- glue::as_glue(names(countries)[i])
  country <- countries[i]
  
  ### Load data ----
  load(paste0("data/Input_", iso, ".rdata"))
  load(paste0("out/ASC/Post_", iso, ".rdata"))
  
  wts <- extract(fitted_as_uni, pars = "inc_a")$inc_a
  wts <- wts[, , dat_as$n_t] * matrix(dat_as$Pop[, dat_as$n_t], dim(wts)[1], 2, byrow = T)
  wts <- wts / rowSums(wts)

  notification <- notification %>% mutate(ISO = iso, Country = country)
  dat_inc <- incidence %>% 
    left_join(rbind(
      notification %>%
        filter(Year == 2019) %>%
        group_by(Country, ISO, Year) %>%
        summarise(Pop = sum(Pop), Sex = "Total"),
      notification %>%
        filter(Year == 2019) %>%
        select(Year, Pop, Sex)
    )) %>%
    mutate(m = m / Pop, u = u / Pop, l = l / Pop, Index = "WHO") %>%
    select(Year, Sex, Index, m, l, u)
  
  
  fore_inc <- local({
    inc_a <- extract(fitted_as_uni, c("inc_a"))$inc_a[, , dat_as$n_t]
    inc_a <- cbind(inc_a, rowSums(inc_a * wts))
    colnames(inc_a) <- c("A_Female", "A_Male", "A_Total")
  
    inc_a <- data.table::data.table(inc_a) %>% 
      pivot_longer(starts_with("A")) %>%
      separate(name, c("Index", "Sex"), "_") %>%
      group_by(Index, Sex) %>%
      summarise(m = mean(value), l = quantile(value, 0.025), u = quantile(value, 0.975), Year = 2018)
    
    inc_s <- extract(fitted_as_uni, c("inc_s"))$inc_s[, , dat_as$n_t]
    inc_s <- cbind(inc_s, rowSums(inc_s * wts))
    colnames(inc_s) <- c("S_Female", "S_Male", "S_Total")
    
    inc_s <- data.table::data.table(inc_s) %>% 
      pivot_longer(starts_with("S")) %>%
      separate(name, c("Index", "Sex"), "_") %>%
      group_by(Index, Sex) %>%
      summarise(m = mean(value), l = quantile(value, 0.025), u = quantile(value, 0.975), Year = 2018)
    
    rbind(inc_s, inc_a)
  })
  
  
  inc <- rbind(fore_inc, dat_inc)
  inc$Country <- country

  save(inc, file = "out/Incidence_" + iso + ".rdata")

  incs[[i]] <- inc  
}

incs <- data.table::rbindlist(incs)

save(incs, file = "out/Incidence_All.rdata")

