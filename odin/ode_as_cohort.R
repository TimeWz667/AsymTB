deriv(A) <- - r_sym * A - ga_a * A
deriv(S) <- r_sym * A - ga_s * S
deriv(Cured) <- r_sc * (S + A)
deriv(Death) <- r_death * S
deriv(Noti) <- S/del


initial(A) <- A0
initial(S) <- S0
initial(Cured) <- Cured0
initial(Death) <- Death0
initial(Noti) <- Noti0


ga_a <- r_sc
ga_s <- r_sc + r_death + 1/del


A0 <- user(1)
S0 <- user(0)
Cured0 <- user(0)
Death0 <- user(0)
Noti0 <- user(0)


# Parameters
r_sym <- user(4) # Symptom activation
r_sc <- user(0.15) # Self-cure
r_death <- user(0.127)
del <- user(0.5) # Delay sp
