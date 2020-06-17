source("R/eval_p_pacf.R")

options(odin.verbose = F,
        odin.target = ifelse(odin::can_compile(), "c", "r"),
        odin.compiler_warnings = F,
        odin.no_check_unused_equations = T
)



n_foreward = 10
n_acf = 3
n_iter <- 500


countries <- c("Malawi", "Kenya")

for (country in countries) {
  load(paste0("data/Input_", country, ".rdata"))
  load(file = paste0("Case", country, "/Fitted.rdata"))
  
  p_pacf_05 <- list()
  for(i in 1:12) {
    p_pacf_05[[i]] <- eval_p_pacf(fit_asym[[i]], pop_foreword = pop_foreword, 0.05, n_foreward, n_acf, n_iter)
    cat("Duration", i, "Completed\n")
  }
  
  
  p_pacf_20 <- list()
  for(i in 1:12) {
    p_pacf_20[[i]] <- eval_p_pacf(fit_asym[[i]], pop_foreword = pop_foreword, 0.2, n_foreward, n_acf, n_iter)
    cat("Duration", i, "Completed\n")
  }
  
  save(p_pacf_05, p_pacf_20, file = paste0("Case", country, "/PACF_population.rdata"))
}
