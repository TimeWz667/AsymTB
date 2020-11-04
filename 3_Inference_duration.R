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


countries_cs <- c(
  KEN = "Kenya",
  MWI = "Malawi", 
  PHL = "Philippines", 
  TZA = "United Republic of Tanzania", 
  UGA = "Uganda", #
  ZMB = "Zambia"
)


### Load data ----

TTE_all <- list()


for (i in 1:length(countries)) {
  iso <- names(countries)[i]
  country <- countries[i]
  
  
  ### Load data ----
  load(paste0("out/ASC/Post_", iso, ".rdata"))
  
  wts <- extract(fitted_as_uni, pars = "inc_a")$inc_a
  wts <- wts[, , dim(wts)[3]]
  wts <- wts / rowSums(wts)
  
  
  ### Duration ----
  if (country %in% countries_cs) {
    ext <- extract(fitted_asc_uni, pars = c("ra", "r_sym", "rs", "r_aware", "rc", "r_det"))
    
    TTE <- with(ext, {
      A <- 1 / (ra + r_sym)
      S <- A + 1 / (rs + r_aware)
      C <- S + 1 / (rc + r_det)
      A <- A[1:500, ]
      S <- S[1:500, ]
      C <- C[1:500, ]
      
      rbind(
        data.table::data.table(
          Country = country,
          Sex = "Female",
          A = A[, 1],
          S = S[, 1],
          C = C[, 1]
        ),
        data.table::data.table(
          Country = country,
          Sex = "Male",
          A = A[, 2],
          S = S[, 2],
          C = C[, 2]
        ),
        data.table::data.table(
          Country = country,
          Sex = "Total",
          A = rowSums(A * wts[1:500, ]),
          S = rowSums(S * wts[1:500, ]),
          C = rowSums(C * wts[1:500, ])
        )
      )
    }) %>%
      mutate(DA = A, DT = C) %>%
      pivot_longer(c(A, S, C), names_to = "Phase", values_to = "Duration")
  } else {
    ext <- extract(fitted_as_uni, pars = c("ra", "r_sym", "rs", "r_det"))
    
    TTE <- with(ext, {
      A <- 1 / (ra + r_sym)
      C <- A + 1 / (rs + r_det)
      A <- A[1:500, ]
      C <- C[1:500, ]
      
      rbind(
        data.table::data.table(
          Country = country,
          Sex = "Female",
          A = A[, 1],
          C = C[, 1]
        ),
        data.table::data.table(
          Country = country,
          Sex = "Male",
          A = A[, 2],
          C = C[, 2]
        ),
        data.table::data.table(
          Country = country,
          Sex = "Total",
          A = rowSums(A * wts[1:500, ]),
          C = rowSums(C * wts[1:500, ])
        )
      )
    }) %>%
      mutate(DA = A, DT = C) %>%
      pivot_longer(c(A, C), names_to = "Phase", values_to = "Duration")
  }
  
  TTE <- TTE %>% mutate(Phase = factor(Phase, levels = c("A", "S", "C")))
  
  TTE_all[[i]] <- TTE
  save(TTE, file = paste0("out/TTE_", iso, ".rdata"))
}

TTE_all <- data.table::rbindlist(TTE_all)

save(TTE_all, file = "out/TTE_All.rdata")



### Compare to P:N ratio


Dur_sim_all <- list()
Dur_data_all <- list()

for (i in 1:length(countries)) {
  iso <- names(countries)[i]
  country <- countries[i]
  
  
  ### Load data ----
  load(paste0("out/ASC/Post_", iso, ".rdata"))
  
  wts <- extract(fitted_as_uni, pars = "inc_a")$inc_a
  wts <- wts[, , dim(wts)[3]]
  wts <- wts / rowSums(wts)
  
  
  ### Duration ----
  
  Dur_data <- with(dat_as, {
    prev1 <- (sum(Asym) + sum(Sym)) / sum(N)
    prev2 <- sum(Sym) / sum(N)
    
    
    noti <- sum(Noti[, Years == YearSurveyed]) / sum(Pop[, Years == YearSurveyed])
    
    data.table::data.table(
      Country = country,
      ISO = iso,
      Def = c("D1", "D2"),
      Type = "Data",
      m = c(prev1 / noti, prev2 / noti)
    )
  })

  
  Dur_sim <- with(extract(fitted_as_uni, pars = c("nr", "prv", "pr_s", "dur_s")), {
    
    ti <- dat_as$Years == dat_as$YearSurveyed
    ti <- ifelse(any(ti), which(ti), 1)
    
    prev1 <- rowSums(prv[, , ti] * wts)
    prev2 <- rowSums(prv[, , ti] * wts * pr_s)
    
    noti <- rowSums(nr[, , ti] * wts)
    
    data.table::data.table(
      Country = country,
      ISO = iso,
      Type = "Fitted",
      D1 = prev1 / noti,
      D2 = mean(prev2 / noti),
      D3 = mean(dur_s)
    )
  }) %>%
    pivot_longer(c(D1, D2, D3), names_to = "Def", values_to = "Duration") %>%
    group_by(Country, ISO, Type, Def) %>%
    summarise(m = mean(Duration), l = quantile(Duration, 0.025), u = quantile(Duration, 0.975))
  
  
  Dur_data_all[[i]] <- Dur_data
  Dur_sim_all[[i]] <- Dur_sim
  
  save(Dur_data, Dur_sim, file = paste0("out/PN_", iso, ".rdata"))
  
}

Dur_data_all <- data.table::rbindlist(Dur_data_all)
Dur_sim_all <- data.table::rbindlist(Dur_sim_all)


save(Dur_data_all, Dur_sim_all, file = paste0("out/PN_All.rdata"))

