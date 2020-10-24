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

