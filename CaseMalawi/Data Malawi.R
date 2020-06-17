## Population/Notification data ----
load("data/Input_Malawi.rdata")
pop <- pop_data

## Prevalence data -----
mat_asym <- rbind( # asymptom + no care seeking intention
  Sym = c(Sn = 70 - 14 - 19, Sp = 62 - 10 - 21),
  Asym = c(Sn = 14 + 19, Sp = 10 + 21)
)

#mat_asym <- rbind( # asymptom
#  Sym = c(Sn = 70 - 19, Sp = 62-21),
#  Asym = c(Sn = 19, Sp = 21)
#)

p_asym <- sum(mat_asym[2, ]) / sum(mat_asym)

p_sp_asym <- sum(mat_asym[2, 2]) / sum(mat_asym[2, ])

p_male_asym <- 16609 / (16609 + 9435)

p_asym_f <- p_asym * (1 - p_male_asym) * c((1 - p_sp_asym), p_sp_asym)
p_asym_m <- p_asym * p_male_asym * c((1 - p_sp_asym), p_sp_asym)

mat_sex <- rbind(
  female = c(Sn = 69 - 28, Sp = 28),
  male = c(Sn = 63 - 34, Sp = 34)
)

mat_female <- round(mat_sex["female", ] * cbind(p_asym_f, 1 - p_asym_f))
mat_male <- round(mat_sex["male", ] * cbind(p_asym_m, 1 - p_asym_m))




## Collect data
data_malawi <- list()
# 
# data_malawi$single_time <- list(
#   female = list(N = 18490, Pop = pop$pop_f$Pop[1], Years = cases$year, year_survey = 2013, 
#              Asym = sum(mat_female[, 1]), Sn = mat_female[1, 2], Sp = mat_female[2, 2], 
#              Noti_Sn = cases$noti_sn_f[1], 
#              Noti_Sp = cases$noti_sp_f[1]),
#   male = list(N = 13089, Pop = pop$pop_m$Pop[1], Years = cases$year, year_survey = 2013,
#            Asym = sum(mat_male[, 1]), Sn = mat_male[1, 2], Sp = mat_male[2, 2], 
#            Noti_Sn = cases$noti_sn_m[1], 
#            Noti_Sp = cases$noti_sp_m[1])
# )


data_malawi$a_sn_sp <- list(
  female = list(N = 18490, Pop = pop$pop_f$Pop, Years = cases$year, year_survey = 2013, 
             Asym = sum(mat_female[, 1]), Sn = mat_female[1, 2], Sp = mat_female[2, 2], 
             Noti_Sn = cases$noti_sn_f, 
             Noti_Sp = cases$noti_sp_f,
             n_t = length(pop$pop_f$Pop)),
  male = list(N = 13089, Pop = pop$pop_m$Pop,Years = cases$year, year_survey = 2013,
           Asym = sum(mat_male[, 1]), Sn = mat_male[1, 2], Sp = mat_male[2, 2], 
           Noti_Sn = cases$noti_sn_m, 
           Noti_Sp = cases$noti_sp_m,
           n_t = length(pop$pop_f$Pop))
)


data_malawi$sn_sp <- list(
  female = list(N = 18490, Pop = pop$pop_f$Pop, Years = cases$year, year_survey = 2013,
                Sn = mat_sex[1, 1], Sp = mat_sex[1, 2], 
                Noti_Sn = cases$noti_sn_f, 
                Noti_Sp = cases$noti_sp_f,
                n_t = length(pop$pop_f$Pop)),
  male = list(N = 13089, Pop = pop$pop_m$Pop, Years = cases$year, year_survey = 2013,
              Sn = mat_sex[2, 1], Sp = mat_sex[2, 2], 
              Noti_Sn = cases$noti_sn_m, 
              Noti_Sp = cases$noti_sp_m,
              n_t = length(pop$pop_f$Pop))
)

input_data <- data_malawi

save(input_data, file = "data/Data_Malawi.rdata")
