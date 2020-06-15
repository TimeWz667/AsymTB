inp <- 
  list(
    r_sym = 12/5, 
    adr = 0.05,
    N0 = 1E5,
    r_tr = 0
  )


f <- "odin/ode_anp_ss.R"
model <- odin::odin(f)

cm <- model(user = inp)
cm$run(seq(0, 5, 0.5))


f <- "odin/ode_np_ss.R"
model <- odin::odin(f)

cm <- model(user = inp)
cm$run(seq(0, 5, 0.5))
