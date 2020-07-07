deriv(A) <- inc - r_sym * A - ga_a * A
deriv(S) <- r_sym * A - r_aware * S - ga_s * S
deriv(C) <- r_aware * S - ga_c * C

initial(A) <- N0 * a0 / (a0 + s0 + c0)
initial(S) <- N0 * s0 / (a0 + s0 + c0)
initial(C) <- N0 * c0 / (a0 + s0 + c0)

   
a0 <- (ga_s + r_aware - adr) / r_sym
s0 <- 1
c0 <- r_aware / (ga_c - adr)

output(N) <- N
output(Noti) <- C / del

output(PropA) <- A / N
output(PropS) <- S / N
output(PropC) <- C / N


N <- A + S + C

inc <- ga_a * A + ga_s * S + ga_c * C - adr * N

ga_a <- r_sc
ga_s <- r_sc + r_death
ga_c <- r_sc + r_death + 1 / del


N0 <- user(1000)

adr <- user(0.01) # Annual decline rate

# Parameters
r_sym <- user(4) # Symptom activation
r_aware <- user(4)
r_sc <- user(0.15)# Self-cure
r_death <- user(0.071)

del <- user(0.5) # Delay sp

