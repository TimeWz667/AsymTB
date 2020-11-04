source("R/steady_states.R")


pars <- get_parameters()
m1 <- model_none(pars)
m2 <- model_asym(pars)

m3 <- model_asym_transition(pars)

res <- data.table::rbindlist(lapply(seq(0, 1, 0.01), function(dur_asym) {
  data.table::rbindlist(lapply(seq(0, 1, 0.01), function(p_sp){
    pars <- get_parameters(p_sp = p_sp, r_sym = 1 / dur_asym, r_tr = 0.2)
    
    if (dur_asym == 0) {
      m1 <- model_none(pars)
      m2 <- model_transition(pars)
    } else {
      m1 <- model_asym(pars)
      m2 <- model_asym_transition(pars)
    }
    names(m1) <- names(m2) <- NULL
    return(data.frame(
      DurAsym = dur_asym, PrSp = p_sp, 
      Conversion = c("No conversion", paste0("Conversion: ", pars["r_tr"])),
      SpSn = c(sum(m1[2:3])/sum(m1[1]), sum(m2[2:3])/sum(m2[1])), 
      ErrDelaySn = c(m1[4]/sum(m1[2]), m2[4]/sum(m2[2])) - pars["delay_sn"], 
      ErrDelaySp = c(m1[5]/m1[3], m2[5]/m2[3]) - pars["delay_sp"],
      ErrDelay = c(sum(m1[4:5])/sum(m1[1:3]), sum(m2[4:5])/sum(m2[1:3])) - pars["delay_sp"],
      ErrDelaySym = c(sum(m1[4:5])/sum(m1[2:3]), sum(m2[4:5])/sum(m2[2:3])) - pars["delay_sp"],
      ErrDur = c(sum(m1[4:5])/sum(m1[1:3]), sum(m2[4:5])/sum(m2[1:3])) - pars["delay_sp"] - dur_asym
      ))
  }))
}))

res$SpSn <- pmin(res$SpSn, 2)


##### Visualisation -----

library(ggplot2)


g <- ggplot(res, aes(x = DurAsym, y = PrSp, z = SpSn)) + 
  geom_tile(aes(fill = SpSn)) + 
  stat_contour(breaks = 1, colour = "black") +
  scale_fill_distiller("Sym/Asym", palette = "Spectral", direction = -1) +
  scale_y_continuous("Proportion of Smear+ve at symptom onset (%)") +
  scale_x_continuous("Duration of asymptomatic phase, year") +
  facet_grid(.~Conversion)

ggsave("output/Sym2Asym.jpg", plot = g, width = 7, height = 4)


ggplot(res, aes(x = DurAsym, y = PrSp, z = ErrDelaySn)) + 
  geom_tile(aes(fill = ErrDelaySn)) + 
  scale_fill_distiller("Error", direction = -1) +
  facet_grid(.~Conversion)


ggplot(res, aes(x = DurAsym, y = PrSp, z = ErrDelaySp)) + 
  geom_tile(aes(fill = ErrDelaySp)) + 
  scale_fill_distiller("Error", direction = -1) +
  facet_grid(.~Conversion)


ggplot(res, aes(x = DurAsym, y = PrSp, z = ErrDelay)) + 
  geom_tile(aes(fill = ErrDelay)) + 
  scale_fill_distiller("Error", direction = -1) +
  facet_grid(.~Conversion)


ggplot(res, aes(x = DurAsym, y = PrSp, z = ErrDelaySym)) + 
  geom_tile(aes(fill = ErrDelaySym)) + 
  scale_fill_distiller("Error", direction = -1) +
  facet_grid(.~Conversion)


g <- ggplot(res, aes(x = DurAsym, y = PrSp, z = ErrDur)) + 
  geom_tile(aes(fill = ErrDur)) + 
  scale_fill_distiller("Error", direction = -1) +
  scale_y_continuous("Proportion of Smear+ve at symptom onset (%)") +
  scale_x_continuous("Duration of asymptomatic phase (year)") +
  facet_grid(.~Conversion)

# ggsave("Output/ErrDur.jpg", plot = g, width = 7, height = 4)

