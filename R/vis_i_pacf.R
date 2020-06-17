extract_idx <- function(ip, dur_asym = 1) {
  data.table::rbindlist(list(
    data.table::data.table(ip$female$None$Indices, ACF = "None", Sex = "Female"),
    data.table::data.table(ip$male$None$Indices, ACF = "None", Sex = "Male"),
    data.table::data.table(ip$female$Sp$Indices, ACF = "Smear+", Sex = "Female"),
    data.table::data.table(ip$male$Sp$Indices, ACF = "Smear+", Sex = "Male"),
    data.table::data.table(ip$female$Sym$Indices, ACF = "Symptom+", Sex = "Female"),
    data.table::data.table(ip$male$Sym$Indices, ACF = "Symptom+", Sex = "Male"),
    data.table::data.table(ip$female$All$Indices, ACF = "All", Sex = "Female"),
    data.table::data.table(ip$male$All$Indices, ACF = "All", Sex = "Male")
  )) %>% 
    mutate(ACF = factor(ACF, levels = c("None", "Smear+", "Symptom+", "All")), DurAsym = dur_asym)
}


visualise_i_pacf_bi <- function(ip1, ip2, dur1, dur2) {
  require(ggplot2)
  require(tidyverse)
  
  fn_ss <- function(ys, ACF, Sex) {
    as.tibble(ys) %>% 
      mutate(t = (1:nrow(.) - 1) / 12, Noti = NotiSp + NotiSn, Sym = Sp + Sn, Asym = A,
             ACF = 1 - Noti - Sym - Asym - Cured - Death) %>%
      filter(t < 2) %>% select(t, ACF, Asym, Sym, Noti, Cured, Death) %>%
      pivot_longer(-t, values_to = "n", names_to = "name") %>%
      mutate(ACF = ACF, Sex = Sex)
  }
  
  
  extract_idx <- function(ip, dur_asym = 1) {
    res <- list()
    
    res$Indices <- data.table::rbindlist(list(
      data.table::data.table(ip$female$None$Indices, ACF = "None", Sex = "Female"),
      data.table::data.table(ip$male$None$Indices, ACF = "None", Sex = "Male"),
      data.table::data.table(ip$female$Sp$Indices, ACF = "Smear+", Sex = "Female"),
      data.table::data.table(ip$male$Sp$Indices, ACF = "Smear+", Sex = "Male"),
      data.table::data.table(ip$female$Sym$Indices, ACF = "Symptom+", Sex = "Female"),
      data.table::data.table(ip$male$Sym$Indices, ACF = "Symptom+", Sex = "Male"),
      data.table::data.table(ip$female$All$Indices, ACF = "All", Sex = "Female"),
      data.table::data.table(ip$male$All$Indices, ACF = "All", Sex = "Male")
    )) %>% 
      mutate(ACF = factor(ACF, levels = c("None", "Smear+", "Symptom+", "All")), DurAsym = dur_asym)
    
    res$Ys <- data.table::rbindlist(list(
      fn_ss(ip$female$None$Ys, "None", "Female"),
      fn_ss(ip$female$Sp$Ys, "Smear+", "Female"),
      fn_ss(ip$female$Sym$Ys, "Symptom+", "Female"),
      fn_ss(ip$female$All$Ys, "All", "Female"),
      fn_ss(ip$male$None$Ys, "None", "Male"),
      fn_ss(ip$male$Sp$Ys, "Smear+", "Male"),
      fn_ss(ip$male$Sym$Ys, "Symptom+", "Male"),
      fn_ss(ip$male$All$Ys, "All", "Male")
    )) %>%
      mutate(name = factor(name, levels = c("Cured", "Death", "ACF", "Noti", "Sym", "Asym")),
             ACF = factor(ACF, levels = c("None", "Smear+", "Symptom+", "All")), DurAsym = dur_asym)
    
    res
  }
  
  ## Bind data ---
  idx1 <- extract_idx(ip1, dur1)
  idx2 <- extract_idx(ip2, dur2)
  indices_bind <- rbind(idx1$Indices, idx2$Indices)
  
  ## Visualisation ---
  gs <- list()
  
  gs$g_state_dist_1 <- idx1$Ys %>%  
    ggplot() +
    geom_bar(aes(x = t, y = n, fill = name), position = "stack", stat = "identity") +
    facet_grid(ACF~Sex)
  
  gs$g_state_dist_2 <- idx2$Ys %>%  
    ggplot() +
    geom_bar(aes(x = t, y = n, fill = name), position = "stack", stat = "identity") +
    facet_grid(ACF~Sex)
  
  gs$g_endpoints_1 <- idx1$Indices %>% 
    filter(name %in% c("Cured", "Death",  "NotiSp", "NotiSn", "ACF_A", "ACF_Sn", "ACF_Sp")) %>%
    mutate(name = factor(name, levels = c("Cured", "Death",  "NotiSp", "ACF_Sp", "NotiSn", "ACF_Sn", "ACF_A"))) %>%
    ggplot() +
    geom_bar(aes(x = ACF, y = mean, fill = name), position = "stack", stat = "identity") +
    facet_grid(.~Sex)
  
  gs$g_endpoints_2 <- idx2$Indices %>% 
    filter(name %in% c("Cured", "Death",  "NotiSp", "NotiSn", "ACF_A", "ACF_Sn", "ACF_Sp")) %>%
    mutate(name = factor(name, levels = c("Cured", "Death",  "NotiSp", "ACF_Sp", "NotiSn", "ACF_Sn", "ACF_A"))) %>%
    ggplot() +
    geom_bar(aes(x = ACF, y = mean, fill = name), position = "stack", stat = "identity") +
    facet_grid(.~Sex)
  
  gs$g_endpoints <- indices_bind %>% 
    filter(name %in% c("Cured", "Death",  "NotiSp", "NotiSn", "ACF_A", "ACF_Sn", "ACF_Sp")) %>%
    mutate(name = factor(name, levels = c("Cured", "Death",  "NotiSp", "ACF_Sp", "NotiSn", "ACF_Sn", "ACF_A"))) %>%
    ggplot() +
    geom_bar(aes(x = ACF, y = mean, fill = name), position = "stack", stat = "identity") +
    facet_grid(DurAsym~Sex)
  
  
  gs$g_dur_inf_1 <- idx1$Indices %>% filter(name %in% c("Dur", "Inf")) %>%
    mutate(name = factor(name, levels = c("Dur", "Inf"))) %>%
    ggplot() +
    geom_pointrange(aes(x=ACF, y=m, ymin=l, ymax=u, colour = Sex), position = position_dodge(0.2)) +
    expand_limits(y = 0) +
    facet_grid(name~.)
  
  gs$g_dur_inf_2 <- idx2$Indices %>% filter(name %in% c("Dur", "Inf")) %>%
    mutate(name = factor(name, levels = c("Dur", "Inf"))) %>%
    ggplot() +
    geom_pointrange(aes(x=ACF, y=m, ymin=l, ymax=u, colour = Sex), position = position_dodge(0.2)) +
    expand_limits(y = 0) +
    facet_grid(name~.)
  
  gs$g_dur_inf <- indices_bind %>% filter(name %in% c("Dur", "Inf")) %>%
    mutate(name = factor(name, levels = c("Dur", "Inf"))) %>%
    ggplot() +
    geom_pointrange(aes(x=ACF, y=m, ymin=l, ymax=u, colour = Sex), position = position_dodge(0.2)) +
    expand_limits(y = 0) +
    facet_grid(name~DurAsym)
  
  gs$g_dur_1 <- idx1$Indices %>% filter(startsWith(name, "Dur_")) %>%
    ggplot() +
    geom_bar(aes(x = ACF, y = mean, fill = name), position = "stack", stat = "identity") +
    geom_pointrange(data = idx1$Indices %>% filter(name == "Dur"), aes(x=ACF, y=m, ymin=l, ymax=u)) +
    facet_grid(.~Sex)
  
  gs$g_dur_2 <- idx2$Indices %>% filter(startsWith(name, "Dur_")) %>%
    ggplot() +
    geom_bar(aes(x = ACF, y = mean, fill = name), position = "stack", stat = "identity") +
    geom_pointrange(data = idx2$Indices %>% filter(name == "Dur"), aes(x=ACF, y=m, ymin=l, ymax=u)) +
    facet_grid(.~Sex)
  
  gs$g_dur <- indices_bind %>% filter(startsWith(name, "Dur_")) %>%
    ggplot() +
    geom_bar(aes(x = ACF, y = mean, fill = name), position = "stack", stat = "identity") +
    geom_pointrange(data = indices_bind %>% filter(name == "Dur"), aes(x=ACF, y=m, ymin=l, ymax=u)) +
    facet_grid(DurAsym~Sex)
  
  gs$g_inf_1 <- idx1$Indices %>% filter(startsWith(name, "Inf_")) %>%
    ggplot() +
    geom_bar(aes(x = ACF, y = mean, fill = name), position = "stack", stat = "identity") +
    geom_pointrange(data = idx1$Indices %>% filter(name == "Inf"), aes(x=ACF, y=m, ymin=l, ymax=u)) +
    facet_grid(.~Sex)
  
  gs$g_inf_2 <- idx2$Indices %>% filter(startsWith(name, "Inf_")) %>%
    ggplot() +
    geom_bar(aes(x = ACF, y = mean, fill = name), position = "stack", stat = "identity") +
    geom_pointrange(data = idx2$Indices %>% filter(name == "Inf"), aes(x=ACF, y=m, ymin=l, ymax=u)) +
    facet_grid(.~Sex)
  
  gs$g_inf <- indices_bind %>% filter(startsWith(name, "Inf_")) %>%
    ggplot() +
    geom_bar(aes(x = ACF, y = mean, fill = name), position = "stack", stat = "identity") +
    geom_pointrange(data = indices_bind %>% filter(name == "Inf"), aes(x=ACF, y=m, ymin=l, ymax=u)) +
    facet_grid(DurAsym~Sex)
  
  gs
}


visualise_i_pacf_12 <- function(i_pacf) {
  idx <- data.table::rbindlist(lapply(1:12, function(dur_asym) {
    ip <- i_pacf[[dur_asym]]
    
    data.table::rbindlist(list(
      data.table::data.table(ip$female$None$Indices, ACF = "None", Sex = "Female"),
      data.table::data.table(ip$male$None$Indices, ACF = "None", Sex = "Male"),
      data.table::data.table(ip$female$Sp$Indices, ACF = "Smear+", Sex = "Female"),
      data.table::data.table(ip$male$Sp$Indices, ACF = "Smear+", Sex = "Male"),
      data.table::data.table(ip$female$Sym$Indices, ACF = "Symptom+", Sex = "Female"),
      data.table::data.table(ip$male$Sym$Indices, ACF = "Symptom+", Sex = "Male"),
      data.table::data.table(ip$female$All$Indices, ACF = "All", Sex = "Female"),
      data.table::data.table(ip$male$All$Indices, ACF = "All", Sex = "Male")
    )) %>% 
      mutate(ACF = factor(ACF, levels = c("None", "Smear+", "Symptom+", "All")), 
             DurAsym = dur_asym)
  }))
  gs <- list()
  
  
  gs$g_dur <- idx %>% filter(name == "Dur") %>%
    ggplot() +
    geom_pointrange(aes(x = DurAsym, y = m, ymin = l, ymax = u, colour = ACF)) +
    facet_grid(.~Sex) +
    scale_y_continuous("Duration [TB activation, notification)") +
    scale_x_continuous("Duration of asymptomatic phase, month", breaks = c(0, 4, 8, 12)) +
    expand_limits(y = 0, x = 0)
  
  gs$g_death <- idx %>% filter(name == "Death") %>%
    ggplot() +
    geom_pointrange(aes(x = DurAsym, y = m, ymin = l, ymax = u, colour = ACF)) +
    facet_grid(.~Sex) +
    scale_y_continuous("Untreated deaths, %") +
    scale_x_continuous("Duration of asymptomatic phase, month", breaks = c(0, 4, 8, 12)) +
    expand_limits(y = 0, x = 0)
  
  gs$g_inf <- idx %>% filter(name == "Inf") %>%
    ggplot() +
    geom_pointrange(aes(x = DurAsym, y = m, ymin = l, ymax = u, colour = ACF)) +
    facet_grid(.~Sex) +
    scale_y_continuous("Transmissibility, per contact-year") +
    scale_x_continuous("Duration of asymptomatic phase, month", breaks = c(0, 4, 8, 12)) +
    expand_limits(y = 0, x = 0)
  
  gs
}

