library(rstan)


source("R/fit.R")

options(mc.cores = parallel::detectCores())
options(odin.verbose = F,
        odin.target = ifelse(odin::can_compile(), "c", "r"),
        odin.compiler_warnings = F,
        odin.no_check_unused_equations = T
)

rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = "-march=corei7 -mtune=corei7")




n_iter <- 1E4
n_collect <- 1000

countries <- c("Malawi", "Kenya")

for (country in countries) {
  load(paste0("data/Data_", country, ".rdata"))
  
  fit_baseline <- fit_np(input_data$sn_sp, n_iter = n_iter, n_collect = n_collect)

  fit_asym <- lapply(1:12, function(dur_asym) {
    res <- suppressWarnings(fit_anp(input_data$a_sn_sp, dur_asym = dur_asym, n_iter = n_iter))
    cat("Asymptomatic period: ", dur_asym, "\n")
    return(res)
  }) 
  
  ts.plot(t(sapply(fit_asym, function(x) quantile(x$female$Duration[, "All"]))), xlim = c(0, 12), title = "Female")
  lines(c(0, 0, 0, 0, 0), quantile(fit_baseline$female$Duration[, "All"]))
  ts.plot(t(sapply(fit_asym, function(x) quantile(x$male$Duration[, "All"]))), xlim = c(0, 12), title = "Male")
  lines(c(0, 0, 0, 0, 0), quantile(fit_baseline$male$Duration[, "All"]))
  
  save(fit_baseline, fit_asym, file = paste0("Case", country, "/Fitted.rdata"))
} 


## refit if needed -----
country <- "Kenya"
refit <- c(5)
load(file = paste0("Case", country, "/Fitted.rdata"))
load(paste0("data/Data_", country, ".rdata"))

for (dur_asym in refit) {
  fit_asym[[dur_asym]] <- suppressWarnings(fit_anp(input_data$a_sn_sp, dur_asym = dur_asym, n_iter = n_iter))
  cat("Asymptomatic period: ", dur_asym, "\n")
}


ts.plot(ts(t(sapply(fit_asym, function(x) quantile(x$female$Duration[, "All"])))), xlim = c(0, 12), title = "Female")
lines(c(0, 0, 0, 0, 0), quantile(fit_baseline$female$Duration[, "All"]))
ts.plot(ts(t(sapply(fit_asym, function(x) quantile(x$male$Duration[, "All"])))), xlim = c(0, 12), title = "Male")
lines(c(0, 0, 0, 0, 0), quantile(fit_baseline$male$Duration[, "All"]))

