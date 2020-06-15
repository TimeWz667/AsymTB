deriv(A) <- - r_sym * A - ga_a * A
deriv(Sn) <- (1 - p_sp) * r_sym * A - r_tr * Sn - ga_sn * Sn
deriv(Sp) <- p_sp * r_sym * A + r_tr * Sn - ga_sp * Sp
deriv(Cured) <- (Sn + Sp + A) * r_sc
deriv(Death) <- r_death_sn * Sn + r_death_sp * Sp
deriv(NotiSn) <- Sn/del_sn 
deriv(NotiSp) <- Sp/del_sp 


initial(A) <- A0
initial(Sn) <- Sn0
initial(Sp) <- Sp0
initial(Cured) <- Cured0
initial(Death) <- Death0
initial(NotiSn) <- NotiSn0
initial(NotiSp) <- NotiSp0

ga_a <- r_sc
ga_sn <- r_sc + r_death_sn + 1/del_sn
ga_sp <- r_sc + r_death_sp + 1/del_sp


A0 <- user(1)
Sn0 <- user(0)
Sp0 <- user(0)
Cured0 <- user(0)
Death0 <- user(0)
NotiSn0 <- user(0)
NotiSp0 <- user(0)

# Parameters
r_sym <- user(4) # Symptom activation
p_sp <- user(0.3) # sp at symptom onset
r_tr <- user(0.3) # Conversion
r_sc <- user(0.15) # Self-cure
r_death_sp <- user(0.127)
r_death_sn <- user(0.024)
del_sp <- user(0.5) # Delay sp
del_sn <- user(0.8) # Delay sn
