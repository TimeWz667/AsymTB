deriv(A) <- - r_sym * A - ga_a * A
deriv(S) <- r_sym * A - r_det * S - ga_s * S
deriv(Tr) <- r_det * S - (r_death_ontr + r_cure + r_ltfu) * Tr
deriv(LTFU) <- r_ltfu * Tr
deriv(Death) <- r_death_a * A + r_death_s * S + r_death_ontr * Tr
deriv(SelfCured) <- r_sc * (A + S)
deriv(Cured) <- r_cure * Tr



initial(A) <- N0
initial(S) <- 0
initial(Tr) <- 0
initial(LTFU) <- 0
initial(Death) <- 0
initial(SelfCured) <- 0
initial(Cured) <- 0


ga_a <- r_sc + r_death_a
ga_s <- r_sc + r_death_s

N0 <- user(1)



# Parameters
r_sym <- user(4) # Symptom activation

r_det <- user(0.5) # rate to case-detection


r_sc <- user(0.15) # Self-curer_death_bg <- user(0)
r_death_a <- user(0)
r_death_s <- user(0)

r_cure <- user(2)
r_ltfu <- user(0.3)
r_death_ontr <- user(0.1)
