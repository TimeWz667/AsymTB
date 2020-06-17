source("Source/CohortPACF.R")

options(odin.verbose = F,
        odin.target = ifelse(odin::can_compile(), "c", "r"),
        odin.compiler_warnings = F,
        odin.no_check_unused_equations = T
)

load(file = "output/Malawi/Fitted.rdata")


n_iter <- 500
beta <- c(0.22, 0.22, 1)


sen <- 0.05

c_pacf_05 <- list()
for(i in 1:12) {
  c_pacf_05[[i]] <- eval_c_pacf(fit_asym[[i]], i, sen, beta, n_iter)
  cat("Duration", i, "Completed\n")
}


sen <- 0.2

c_pacf_20 <- list()
for(i in 1:12) {
  c_pacf_20[[i]] <- eval_c_pacf(fit_asym[[i]], i, sen, beta, n_iter)
  cat("Duration", i, "Completed\n")
}


save(c_pacf_05, c_pacf_20, file = "Output/Malawi/CohortEffects.rdata")
