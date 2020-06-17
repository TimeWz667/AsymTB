
model_none <- function(pars) {
  with(as.list(pars), {
    rn <- r_sc + r_death_sn + 1 / del_sn - adr
    rp <- r_sc + r_death_sp + 1 / del_sp - adr
    
    
    sn0 <- (1 - p_sp) * rp 
    sp0 <- p_sp * rn 
    
    sn <- sn0 / (sn0 + sp0)
    sp <- sp0 / (sn0 + sp0)
    
    c(A = 0, Sn =sn, Sp = sp, NotiSn = sn/del_sn, NotiSp = sp/del_sp)
  })
}


model_asym <- function(pars) {
  with(as.list(pars), {
    ra <- r_sc
    rn <- r_sc + r_death_sn + 1 / del_sn - adr
    rp <- r_sc + r_death_sp + 1 / del_sp - adr
    
    sn0 <- (1 - p_sp) * r_sym / rn    
    sp0 <- p_sp * r_sym / rp
    
    a <- 1 / (1 + sn0 + sp0)
    
    c(A = a, Sn = sn0 * a, Sp = sp0 * a, NotiSn = sn0 * a/del_sn, NotiSp = sp0 * a/del_sp)
  })
}


model_transition <- function(pars) {
  with(as.list(pars), {
    rn <- r_sc + r_death_sn + 1 / del_sn - adr
    rp <- r_sc + r_death_sp + 1 / del_sp - adr
    
    sn0 <- (1 - p_sp) * rp 
    sp0 <- p_sp * rn + r_tr 
    
    sn <- sn0 / (sn0 + sp0)
    sp <- sp0 / (sn0 + sp0)
    
    c(A = 0, Sn =sn, Sp = sp, NotiSn = sn/del_sn, NotiSp = sp/del_sp)
  })
}


model_asym_transition <- function(pars) {
  with(as.list(pars), {
    ra <- r_sc
    rn <- r_sc + r_death_sn + 1 / del_sn - adr
    rp <- r_sc + r_death_sp + 1 / del_sp - adr
    
    sn0 <- (1 - p_sp) * r_sym / (r_tr + rn)    
    sp0 <- p_sp * r_sym / (rp) + r_tr * sn0 / (rp)
    
    a <- 1 / (1 + sn0 + sp0)
    
    c(A = a, Sn = sn0 * a, Sp = sp0 * a, NotiSn = sn0 * a/del_sn, NotiSp = sp0 * a/del_sp)
  })
}


get_parameters <- function(r_sym = 4, 
                           p_sp = 0.3,
                           r_tr = 0.312, # qexp(pexp(0.006) * 52)
                           r_sc = 0.15, 
                           r_death_sp = 0.127, 
                           r_death_sn = 0.024,
                           del_sp = 1.2, 
                           del_sn = 1,
                           adr = 0) {
  c(
    r_sym = r_sym, 
    p_sp = p_sp,
    r_tr = r_tr,
    r_sc = r_sc, 
    r_death_sp = r_death_sp, 
    r_death_sn = r_death_sn,
    del_sp = del_sp, 
    del_sn = del_sn,
    adr = adr
  )
}

