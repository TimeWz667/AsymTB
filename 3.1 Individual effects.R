source("R/eval_i_pacf.R")

options(odin.verbose = F,
        odin.target = ifelse(odin::can_compile(), "c", "r"),
        odin.compiler_warnings = F,
        odin.no_check_unused_equations = T
)



n_iter <- 500
beta <- c(0.22, 0.22, 1)


countries <- c("Malawi", "Kenya")

for (country in countries) {
  load(file = paste0("Case", country, "/Fitted.rdata"))
  
  i_pacf <- list()
  for(i in 1:12) {
    i_pacf[[i]] <- eval_i_pacf(fit_asym[[i]], i, beta, n_iter)
    cat("Duration", i, "Completed\n")
  }
  
  save(i_pacf, file = paste0("Case", country, "/PACF_individual.rdata"))
}
