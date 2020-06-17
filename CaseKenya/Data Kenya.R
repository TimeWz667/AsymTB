## Population/Notification data ----
load("data/Input_Kenya.rdata")
pop <- pop_data

## Prevalence data -----

prev <- read.csv("data/kenya_prev.csv")

table(prev$Smearpositive, prev$Age >= 15)
table(prev$Culturepositive, prev$Age >= 15)

table(prev$HIVStattus)
table(prev$OtherSymptoms, prev$Culturepositive)


enrolled <- subset(prev, OtherSymptomsWeeksNumber == 0)


enrolled$sym <- (rowSums(enrolled[, c("Cough", "ChestPains", "Coughing", "NightSweats", "Fatigue", 
                                      "Sputum", "Fever", "BreatheShortness", "WeightLoss", 
                                      "OtherSymptoms", "BloodCough")], na.rm=T) > 0) + 0

enrolled$sp <- (enrolled$Smearpositive == "POS" & enrolled$Culturepositive != "CON" & 
  enrolled$Xpertpositive == "POS") + 0

enrolled$bc <- (enrolled$Culturepositive == "POS" | enrolled$Xpertpositive == "POS") + 0

enrolled$sn <- enrolled$bc - enrolled$sp

temp <- subset(enrolled, sp == 1)
mat_sp <- table(temp$sym, temp$Sex)

temp <- subset(enrolled, sn == 1)
mat_sn <- table(temp$sym, temp$Sex)

prv_f_a_n_p <- list(
  N = sum(enrolled$Sex == "F"),
  Asym = sum(mat_sp["0", "F"]) + sum(mat_sn["0", "F"]),
  Sn = sum(mat_sn["1", "F"]), 
  Sp = sum(mat_sp["1", "F"])
)

prv_m_a_n_p <- list(
  N = sum(enrolled$Sex == "M"),
  Asym = sum(mat_sp["0", "M"]) + sum(mat_sn["0", "M"]),
  Sn = sum(mat_sn["1", "M"]), 
  Sp = sum(mat_sp["1", "M"])
)

prv_f_n_p <- list(
  N = sum(enrolled$Sex == "F"),
  Sn = sum(mat_sn[, "F"]), 
  Sp = sum(mat_sp[, "F"])
)

prv_m_n_p <- list(
  N = sum(enrolled$Sex == "M"),
  Sn = sum(mat_sn[, "M"]), 
  Sp = sum(mat_sp[, "M"])
)



## Collect data
data_kenya <- list()
# 
# data_kenya$single_time <- list(
#   female = c(prv_f_a_n_p,
#              Pop = pop$pop_f$Pop[1],
#              Noti_Sn = cases$noti_sn_f[1], 
#              Noti_Sp = cases$noti_sp_f[1]
#              ),
#   male = c(prv_m_a_n_p,
#            Pop = pop$pop_m$Pop[1],
#            Noti_Sn = cases$noti_sn_m[1], 
#            Noti_Sp = cases$noti_sp_m[1]
#   )
# )


data_kenya$a_sn_sp <- list(
  female = c(prv_f_a_n_p,
             list(
               Years = cases$year, year_survey = 2016,
               Pop = pop$pop_f$Pop,
               Noti_Sn = cases$noti_sn_f, 
               Noti_Sp = cases$noti_sp_f,
               n_t = length(pop$pop_f$Pop)
             )
  ),
  male = c(prv_m_a_n_p,
           list(
             Years = cases$year, year_survey = 2016,
             Pop = pop$pop_m$Pop,
             Noti_Sn = cases$noti_sn_m, 
             Noti_Sp = cases$noti_sp_m,
             n_t = length(pop$pop_m$Pop)
           )
  )
)


data_kenya$sn_sp <- list(
  female = c(prv_f_n_p,
             list(
               Years = cases$year, year_survey = 2016,
               Pop = pop$pop_f$Pop,
               Noti_Sn = cases$noti_sn_f, 
               Noti_Sp = cases$noti_sp_f,
               n_t = length(pop$pop_f$Pop)
             )
  ),
  male = c(prv_m_n_p,
           list(
             Years = cases$year, year_survey = 2016,
             Pop = pop$pop_m$Pop,
             Noti_Sn = cases$noti_sn_m, 
             Noti_Sp = cases$noti_sp_m,
             n_t = length(pop$pop_m$Pop)
           )
  )
)

input_data <- data_kenya
save(input_data, file = "data/Data_Kenya.rdata")
