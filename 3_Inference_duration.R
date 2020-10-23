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
  MWI = "Malawi", 
  PHL = "Philippines", 
  TZA = "United Republic of Tanzania", 
  UGA = "Uganda", #
  ZMB = "Zambia"
)



for (i in 1:length(countries)) {
  iso <- names(countries)[i]
  country <- countries[i]
  
  ### Load data ----
  load(paste0("out/ASC/Post_", iso, ".rdata"))
  
  wts <- extract(fitted_as_uni, pars = "inc_a")$inc_a
  wts <- wts[, , dim(wts)[3]]
  wts <- wts / rowSums(wts)
  

  ### Duration ----
  ext <- extract(fitted_as_uni, pars = c("ra", "r_sym", "rs", "r_det"))
  
  Durations_sex <- with(ext, {
    tab <- data.frame(
        D_sym = 1 / (ra + r_sym) , 
        D_care = 1 / (rs + r_det)
      )
    colnames(tab) <- c("D_sym_Female", "D_sym_Male", "D_care_Female", "D_care_Male")
    tab
  }) %>% 
    pivot_longer(everything()) %>%
    separate(name, c("Index", "Stage", "Sex"), "_") %>%
    mutate(Stage = factor(Stage, levels = c("sym", "care")))
  
  
  
  Durations_total <- with(ext, {
    data.frame(D_sym = rowSums(1 / (ra + r_sym) * wts), 
               D_care = rowSums(1 / (rs + r_det) * wts))
  }) %>%
      pivot_longer(starts_with(c("D_")), values_to = "value") %>%
      separate(name, c("Index", "Stage"), "_") %>%
      mutate(Stage = factor(Stage, levels = c("sym", "care")), Sex = "Total")
  
  Durations <- rbind(Durations_sex, Durations_total) %>% arrange(Stage, Sex)
  
  
  TTE_sex <- with(ext, {
    tab <- data.frame(
      T_sym = 1 / (ra + r_sym), 
      T_care = 1 / (ra + r_sym) + 1 / (rs + r_det)
    )
    colnames(tab) <- c("T_sym_Female", "T_sym_Male", "T_care_Female", "T_care_Male")
    tab
  }) %>% 
    pivot_longer(everything()) %>%
    separate(name, c("Index", "Stage", "Sex"), "_") %>%
    mutate(Stage = factor(Stage, levels = c("sym", "care")))
  
  
  
  TTE_total <- with(ext, {
    data.frame(
      T_sym = rowSums(1 / (ra + r_sym) * wts), 
      T_care = rowSums((1 / (ra + r_sym) + 1 / (rs + r_det)) * wts)
    )
  }) %>%
    pivot_longer(starts_with(c("T_")), values_to = "value") %>%
    separate(name, c("Index", "Stage"), "_") %>%
    mutate(Stage = factor(Stage, levels = c("sym", "care")), Sex = "Total")
  
  TTE <- rbind(TTE_sex, TTE_total) %>% arrange(Stage, Sex)
  

  save(Durations, TTE, file = paste0("out/Duration_", iso, ".rdata"))
}





### Load data ----


sim_dur <- data.table::rbindlist(lapply(1:length(countries), function(i) {
  iso <- names(countries)[i]
  country <- countries[i]
  
  load(paste0("out/ASC/Post_", iso, ".rdata"))
  
  wts <- extract(fitted_as_uni, pars = "inc_a")$inc_a
  wts <- wts[, , dim(wts)[3]]
  wts <- wts / rowSums(wts)
  
  if (country %in% countries_cs) {
    ext <- extract(fitted_asc_uni, pars = c("ra", "r_sym", "rs", "r_aware", "rc", "r_det"))
    
    dur <- with(ext, {
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
    
    dur <- with(ext, {
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
  
})) %>%
  mutate(Phase = factor(Phase, levels = c("A", "S", "C")))


save(sim_dur, file = "out/Duration_All.rdata")

