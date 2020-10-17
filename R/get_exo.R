get_exo <- function(hiv, dr, untr = c("no_untr", "as_sn", "full"), bg_death = T,
                    r_sc_l = 0.1, r_sc_u = 0.3, scale_dur = 1, 
                    dr_sn = 0.022, dr_sp = 0.12, dr_hiv_sn = 0.022, dr_hiv_sp = 0.7) {
  
  untr <- match.arg(untr)
  
  
  exo <- list(
    r_sc_l = 0.1,
    r_sc_u = 0.3,
    scale_dur = 1,
    r_death_bg = ifelse(rep(bg_death, length(dr)), dr, 0),
    r_death_sn = switch(untr, 
                        no_untr = rep(0, length(hiv)),
                        ifelse(hiv == "HIV", dr_hiv_sn, dr_sn)),
    r_death_sp = switch(untr, 
                        no_untr = rep(0, length(hiv)),
                        as_sn = ifelse(hiv == "HIV", dr_hiv_sn, dr_sn),
                        full = ifelse(hiv == "HIV", dr_hiv_sp, dr_sp))
  )
  exo
}


get_exo_np <- function(dr, untr = c("no_untr", "as_sn", "full"), bg_death = T,
                        r_sc_l = 0.1, r_sc_u = 0.3, scale_dur = 1, 
                        dr_sn = 0.022, dr_sp = 0.12) {
  
  untr <- match.arg(untr)
  
  n <- length(dr)
  
  exo <- list(
    r_sc_l = 0.1,
    r_sc_u = 0.3,
    scale_dur = 1,
    r_death_bg = ifelse(rep(bg_death, n), dr, 0),
    r_death_sn = switch(untr, 
                        no_untr = rep(0, n),
                        as_sn = rep(dr_sn, n),
                        full = rep(dr_sn, n)),
    r_death_sp = switch(untr, 
                        no_untr = rep(0, n),
                        as_sn = rep(dr_sn, n),
                        full = rep(dr_sp, n))
  )
  exo
}

get_exo_sym <- function(dr, untr = c("no_untr", "as_sn", "avg"), bg_death = T,
                        r_sc_l = 0.1, r_sc_u = 0.3, scale_dur = 1, 
                        dr_sn = 0.022, dr_sp = 0.12) {
  
  untr <- match.arg(untr)
  
  n <- length(dr)
  
  exo <- list(
    r_sc_l = 0.1,
    r_sc_u = 0.3,
    scale_dur = 1,
    r_death_bg = ifelse(rep(bg_death, n), dr, 0),
    r_death_tb = switch(untr, 
                        no_untr = rep(0, n),
                        as_sn = rep(dr_sn, n),
                        avg = rep((dr_sn + dr_sp) / 2, n))
  )
  exo
}
