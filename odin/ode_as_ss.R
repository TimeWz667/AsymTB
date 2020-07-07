deriv(A) <- inc - r_sym * A - ga_a * A
deriv(S) <- r_sym * A - ga_s * S

initial(A) <- N0 * a0 / (s0 + a0)
initial(S) <- N0 * s0 / (s0 + a0)


s0 <- r_sym    
a0 <- ga_s - adr


output(N) <- N
output(Noti) <- S / del

output(PropA) <- A / N
output(PropS) <- S / N


N <- A + S

inc <- ga_a * A + ga_s * S - adr * N

ga_a <- r_sc
ga_s <- r_sc + r_death + 1/del


N0 <- user(1000)

adr <- user(0.01) # Annual decline rate

# Parameters
r_sym <- user(4) # Symptom activation
r_sc <- user(0.15)# Self-cure
r_death <- user(0.127)

del <- user(0.5) # Delay sp

