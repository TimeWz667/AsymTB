get_exo_asc <- function(mor, 
                       untr_a = c("no_untr", "as_sn"), 
                       untr_s = c("no_untr", "as_sn", "avg"), bg_death = T, 
                       r_sc_l = 0.1, r_sc_u = 0.3, scale_dur = 1, 
                       dr_sn = 0.022, dr_sp = 0.12, pr_fp = 0) {
  
  
  untr_a <- match.arg(untr_a)
  untr_s <- match.arg(untr_s)
  
  n <- nrow(mor)
  
  if (bg_death) {
    bg <- mor$DeaR
  } else {
    bg <- rep(0, n)
  }
  
  list(
    r_sc_l = r_sc_l,
    r_sc_u = r_sc_u,
    scale_dur = scale_dur,
    pr_fp = pr_fp,
    r_death_a = switch(untr_a, 
                       no_untr = rep(0, n),
                       as_sn = rep(dr_sn, n)) + bg,
    r_death_s = switch(untr_s, 
                       no_untr = rep(0, n),
                       as_sn = rep(dr_sn, n),
                       avg = mor$pr_sp * dr_sp + (1 - mor$pr_sp) * dr_sn) + bg
  )
}
