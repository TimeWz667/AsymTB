library(tidyverse)
library(rstan)
source("R/inference.R")


### Country list -----
source("data/country_list.R")




### Collect duration ----
Durations_list <- lapply(names(countries), function(iso) {
  iso <- glue::as_glue(iso)
  country <- countries[names(countries) == iso]
  infer_duration(iso, country)
})
names(Durations_list) <- names(countries)

for (iso in names(countries)) {
  iso <- glue::as_glue(iso)
  Durations <- Durations_list[[iso]]
  save(Durations, file = "out/Durations_" + iso + ".rdata")
}



Durations <- list(
  Durations_Sex = bind_rows(lapply(Durations_list, function(x) x$Durations_Sex)),
  Durations_All = bind_rows(lapply(Durations_list, function(x) x$Durations_All))
)

save(Durations, file = "out/Durations_All.rdata")



### Collect conversion rate ----

Conversion_list <- lapply(names(countries), function(iso) {
  iso <- glue::as_glue(iso)
  country <- countries[names(countries) == iso]
  infer_conversion(iso, country)
})
names(Conversion_list) <- names(countries)


for (iso in names(countries)) {
  iso <- glue::as_glue(iso)
  Conversion <- Conversion_list[[iso]]
  save(Conversion, file = "out/Conversion_" + iso + ".rdata")
}



Conversion <- list(
  Reg = lapply(Conversion_list, function(x) x$Reg),
  Pars = bind_rows(lapply(Conversion_list, function(x) x$Pars)),
  Dens = bind_rows(lapply(Conversion_list, function(x) x$Dens)),
  Wts = bind_rows(lapply(Conversion_list, function(x) x$Wts)) %>% mutate(ISO = names(countries))
)

save(Conversion, file = "out/Conversion_All.rdata")


### Collect cascade ----
Cascade_list <- lapply(names(countries), function(iso) {
  iso <- glue::as_glue(iso)
  country <- countries[names(countries) == iso]
  infer_cascade(iso, country)
})
names(Cascade_list) <- names(countries)


Cohort_list <- lapply(names(countries), function(iso) {
  country <- countries[names(countries) == iso]
  infer_cohort(iso, country, 200)
})
names(Cohort_list) <- names(countries)


for (iso in names(countries)) {
  iso <- glue::as_glue(iso)
  Cascade <- Cascade_list[[iso]]
  Cohort <- Cohort_list[[iso]]
  save(Cascade, Cohort, file = "out/Cascade_" + iso + ".rdata")
}



Cascade <- list(
  Cascade_Sex = bind_rows(lapply(Cascade_list, function(x) x$Cascade_Sex)),
  Cascade_All = bind_rows(lapply(Cascade_list, function(x) x$Cascade_All))
)

Cohort <- list(
  Cohort = bind_rows(lapply(Cohort_list, function(x) x$Cohort)),
  End = bind_rows(lapply(Cohort_list, function(x) x$End))
)

save(Cascade, Cohort, file = "out/Cascade_All.rdata")


