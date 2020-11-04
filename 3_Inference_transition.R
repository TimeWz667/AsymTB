rm(list = ls())

library(rstan)


get_transition_matrix <- function(post) {
  mci <- function(xs) {
    c(
      m = mean(xs),
      l = quantile(xs, 0.025),
      u = quantile(xs, 0.975)
    )
  }
  
  pars <- extract(post, c("p_sp", "r_tr", "r_det_sn", "r_det_sp", "r_sym", 
                          "r_sc", "ra", "rn", "rp"))
    
  res <- data.table::rbindlist(list(
    data.table::data.table(
      From = "A",
      To = c("Sn", "Sp", "SC", "Death"),
      with(pars, {
        temp <- cbind((1 - p_sp) * r_sym, p_sp * r_sym, r_sc, ra - r_sc)
        temp <- temp / rowSums(temp)
        temp <- t(apply(temp, 2, mci))
        colnames(temp) <- c("m", "l", "u")
        temp
      })
    ),
    data.table::data.table(
      From = "Sn",
      To = c("Sp", "SC", "Death", "Notification"),
      with(pars, {
        temp <- cbind(r_tr, r_sc, rn - r_sc, r_det_sn)
        temp <- temp / rowSums(temp)
        temp <- t(apply(temp, 2, mci))
        colnames(temp) <- c("m", "l", "u")
        temp
      })
    ),
    data.table::data.table(
      From = "Sp",
      To = c("SC", "Death", "Notification"),
      with(pars, {
        temp <- cbind(r_sc, rp - r_sc, r_det_sp)
        temp <- temp / rowSums(temp)
        temp <- t(apply(temp, 2, mci))
        colnames(temp) <- c("m", "l", "u")
        temp
      })
    )
  ))
  res$To <- factor(res$To, levels = c("Sn", "Sp", "Notification", "SC", "Death"))
  res$From <- factor(res$From, levels = c("A", "Sn", "Sp"))
  res
}


load("out/Full/KEN/Total.rdata")

trmap_ken <- get_transition_matrix(fitted_total)


load("out/Full/BLT/Total.rdata")

trmap_blt <- get_transition_matrix(fitted_total)


save(trmap_blt, trmap_ken, file = "out/TransMap.rdata")
