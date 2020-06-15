deriv(A) <- inc - r_sym * A - ga_a * A
deriv(Sn) <- (1 - p_sp) * r_sym * A - r_tr * Sn - ga_sn * Sn
deriv(Sp) <- p_sp * r_sym * A + r_tr * Sn - ga_sp * Sp


initial(A) <- a * N0
initial(Sn) <- sn0 * a * N0
initial(Sp) <- sp0 * a * N0


sn0 <- (1 - p_sp) * r_sym / (r_tr + ga_sn - adr)    
sp0 <- (p_sp * r_sym + r_tr * sn0) / (ga_sp - adr)
a <- 1 / (1 + sn0 + sp0)


output(N) <- N
output(Noti_Sn) <- Sn / del_sn
output(Noti_Sp) <- Sp / del_sp
output(PropA) <- A / N
output(PropSn) <- Sn / N
output(PropSp) <- Sp / N


N <- A + Sn + Sp

inc <- ga_a * A + ga_sn * Sn + ga_sp * Sp - adr * N

ga_a <- r_sc
ga_sn <- r_sc + r_death_sn + 1/del_sn
ga_sp <- r_sc + r_death_sp + 1/del_sp


N0 <- user(1)

adr <- user(0) # Annual decline rate

# Parameters
r_sym <- user(4) # Symptom activation
p_sp <- user(0.3) # sp at symptom onset
r_tr <- user(0.3) # Conversion
r_sc <- user(0.15)# Self-cure
r_death_sp <- user(0.127)
r_death_sn <- user(0.024)
del_sp <- user(0.5) # Delay sp
del_sn <- user(0.8) # Delay sn
