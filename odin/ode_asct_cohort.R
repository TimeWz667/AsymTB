deriv(A) <- - r_sym * A - ga_a * A
deriv(S) <- r_sym * A - r_aware * S - ga_s * S
deriv(C) <- r_aware * S - r_det * C - ga_c * C
deriv(Tr) <- r_det * C - (r_death_ontr + r_cure + r_ltfu) * Tr
deriv(LTFU) <- r_ltfu * Tr
deriv(Death) <- r_death_a * A + r_death_s * (S + C) + r_death_ontr * Tr
deriv(SelfCured) <- r_sc * (A + S + C)
deriv(Cured) <- r_cure * Tr



initial(A) <- N0
initial(S) <- 0
initial(C) <- 0
initial(Tr) <- 0
initial(LTFU) <- 0
initial(Death) <- 0
initial(SelfCured) <- 0
initial(Cured) <- 0



ga_a <- r_sc + r_death_a
ga_s <- r_sc + r_death_s
ga_c <- r_sc + r_death_s

N0 <- user(1)



# Parameters
r_sym <- user(4) # Symptom activation
r_aware <- user(4)
r_sc <- user(0.15) # Self-cure

r_det <- user(1) # Delay sp

r_cure <- user(2)
r_ltfu <- user(0.3)

r_death_a <- user(0)
r_death_s <- user(0.071)
r_death_ontr <- user(0.1)
