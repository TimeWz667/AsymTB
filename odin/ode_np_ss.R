deriv(Sn) <- (1 - p_sp) * inc - r_tr * Sn - ga_sn * Sn
deriv(Sp) <- p_sp * inc + r_tr * Sn - ga_sp * Sp


initial(Sn) <- sn0 / (sn0 + sp0) * N0
initial(Sp) <- sp0 / (sn0 + sp0) * N0


sn0 <- (1 - p_sp) * (ga_sp - adr) 
sp0 <- p_sp * (ga_sn - adr)  + r_tr 


output(N) <- N
output(Noti_Sn) <- Sn / del_sn
output(Noti_Sp) <- Sp / del_sp
output(PropSn) <- Sn / N
output(PropSp) <- Sp / N


N <- Sn + Sp

inc <- ga_sn * Sn + ga_sp * Sp - adr * N

ga_sn <- r_sc + r_death_sn + 1/del_sn
ga_sp <- r_sc + r_death_sp + 1/del_sp

N0 <- user(1)


adr <- user(0)

r_sym <- user(4)
p_sp <- user(0.3)
r_tr <- user(0.3)
r_sc <- user(0.15) 
r_death_sp <- user(0.127)
r_death_sn <- user(0.024)
del_sp <- user(0.15) 
del_sn <- user(0.15)
