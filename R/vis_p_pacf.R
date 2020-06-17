visualise_p_pacf <- function(p_pacf) {
  require(ggplot2)
  require(tidyverse)
  
  flows <- data.table::rbindlist(lapply(names(p_pacf$female$None$Flow), function(key) {
    data.table::rbindlist(list(
      data.table::data.table(p_pacf$female$None$Flow[[key]], ACF = "None", Sex = "Female"),
      data.table::data.table(p_pacf$male$None$Flow[[key]], ACF = "None", Sex = "Male"),
      data.table::data.table(p_pacf$female$Sp$Flow[[key]], ACF = "Smear+", Sex = "Female"),
      data.table::data.table(p_pacf$male$Sp$Flow[[key]], ACF = "Smear+", Sex = "Male"),
      data.table::data.table(p_pacf$female$Sym$Flow[[key]], ACF = "Symptom+", Sex = "Female"),
      data.table::data.table(p_pacf$male$Sym$Flow[[key]], ACF = "Symptom+", Sex = "Male"),
      data.table::data.table(p_pacf$female$All$Flow[[key]], ACF = "All", Sex = "Female"),
      data.table::data.table(p_pacf$male$All$Flow[[key]], ACF = "All", Sex = "Male")
    ))
  }))
  flows$ACF <- factor(flows$ACF, levels = c("None", "Smear+", "Symptom+", "All"))
  
  
  aocs <- data.table::rbindlist(lapply(names(p_pacf$female$None$AOC), function(key) {
    data.table::rbindlist(list(
      data.table::data.table(p_pacf$female$None$AOC[[key]], ACF = "None", Sex = "Female"),
      data.table::data.table(p_pacf$male$None$AOC[[key]], ACF = "None", Sex = "Male"),
      data.table::data.table(p_pacf$female$Sp$AOC[[key]], ACF = "Smear+", Sex = "Female"),
      data.table::data.table(p_pacf$male$Sp$AOC[[key]], ACF = "Smear+", Sex = "Male"),
      data.table::data.table(p_pacf$female$Sym$AOC[[key]], ACF = "Symptom+", Sex = "Female"),
      data.table::data.table(p_pacf$male$Sym$AOC[[key]], ACF = "Symptom+", Sex = "Male"),
      data.table::data.table(p_pacf$female$All$AOC[[key]], ACF = "All", Sex = "Female"),
      data.table::data.table(p_pacf$male$All$AOC[[key]], ACF = "All", Sex = "Male")
    ))
  }))
  aocs$ACF <- factor(aocs$ACF, levels = c("None", "Smear+", "Symptom+", "All"))
  
  
  avers <- data.table::rbindlist(lapply(names(p_pacf$female$d10), function(key) {
    data.table::rbindlist(list(
      data.table::data.table(p_pacf$female$d10[[key]], ACF = "Smear+", Sex = "Female"),
      data.table::data.table(p_pacf$male$d10[[key]], ACF = "Smear+", Sex = "Male"),
      data.table::data.table(p_pacf$female$d20[[key]], ACF = "Symptom+", Sex = "Female"),
      data.table::data.table(p_pacf$male$d20[[key]], ACF = "Symptom+", Sex = "Male"),
      data.table::data.table(p_pacf$female$d30[[key]], ACF = "All", Sex = "Female"),
      data.table::data.table(p_pacf$male$d30[[key]], ACF = "All", Sex = "Male")
    ))
  }))
  avers <- tibble::as_tibble(avers)
  avers$ACF <- factor(avers$ACF, levels = c("Smear+", "Symptom+", "All"))
  
  
  
  gs_flow <- lapply(unique(flows$name), function(k) {
    flows %>% filter(name == k) %>%
      ggplot(aes(x = t)) +
      geom_line(aes(y = mean, colour = ACF)) +
      geom_ribbon(aes(ymin = l, ymax = u, fill = ACF), alpha = 0.2) +
      scale_y_continuous("Number") +
      scale_x_continuous("Year") +
      facet_grid(.~Sex) +
      expand_limits(y = 0)
  })
  names(gs_flow) <- paste0("g_ts_", unique(flows$name))
  
  
  gs_aoc <- lapply(unique(aocs$name), function(k) {
    aocs %>% filter(name == k) %>%
      ggplot(aes(x = t)) +
      geom_line(aes(y = mean, colour = ACF)) +
      geom_ribbon(aes(ymin = l, ymax = u, fill = ACF), alpha = 0.2) +
      scale_y_continuous("Number") +
      scale_x_continuous("Year") +
      facet_grid(.~Sex) +
      expand_limits(y = 0)
  })
  names(gs_aoc) <- paste0("g_aoc_", unique(aocs$name))
  
  
  gs <- c(gs_flow, gs_aoc)
  
  gs$g_avert <- avers %>% 
    filter(name %in% c("Death", "Act")) %>%
    ggplot(aes(x = t)) +
    geom_line(aes(y = mean, colour = ACF)) +
    geom_ribbon(aes(ymin = l, ymax = u, fill = ACF), alpha = 0.2) +
    scale_y_continuous("Averted events, per 100 000", labels = function(x) x * 1E5) +
    scale_x_continuous("Year") +
    facet_grid(name~Sex, scale = "free_y") +
    expand_limits(y = 0)
  
  
  gs$g_occur <- flows %>%
    filter(name %in% c("Act", "Noti", "Mor", "ACF")) %>%
    mutate(name = ifelse(name == "ACF", "Found", name)) %>%
    select(t, name, mean, ACF, Sex) %>%
    pivot_wider(names_from = name, values_from = mean) %>%
    mutate(Notification = Noti + Found) %>%
    rename(Incidence = Act, Mortality= Mor, `Care Sought` = Noti) %>%
    select(-Found) %>% 
    pivot_longer(c(Incidence, Mortality, Notification, `Care Sought`), values_to = "m") %>%
    mutate(name = factor(name, levels = c("Incidence", "Care Sought", "Notification", "Mortality"))) %>%
    ggplot() +
    geom_line(aes(x = t, y = m, colour = name)) +
    scale_color_discrete("Index") +
    scale_y_continuous("Event rate, per 100 000 capita-year", labels = function(x) x * 1E5) +
    scale_x_continuous("Year", labels = function(x) x + 2020) + 
    facet_grid(Sex~ACF, scale = "free_y") +
    theme(legend.position = "bottom", axis.text.x.bottom = element_text(angle = 45, hjust = 1))
  
  gs
}


visualise_p_pacf <- function(p_pacf_1, p_pacf_2, lab_1, lab_2) {
  gs <- list()
  
  
  flows <- data.table::rbindlist(lapply(c("Act", "Noti", "Mor", "ACF"), function(key) {
    data.table::rbindlist(list(
      data.table::data.table(p_pacf_1$total$None$Flow[[key]], ACF = "None", GP = lab_1),
      data.table::data.table(p_pacf_2$total$None$Flow[[key]], ACF = "None", GP = lab_2),
      data.table::data.table(p_pacf_1$total$Sp$Flow[[key]], ACF = "Smear+", GP = lab_1),
      data.table::data.table(p_pacf_2$total$Sp$Flow[[key]], ACF = "Smear+", GP = lab_2),
      data.table::data.table(p_pacf_1$total$Sym$Flow[[key]], ACF = "Symptom+", GP = lab_1),
      data.table::data.table(p_pacf_2$total$Sym$Flow[[key]], ACF = "Symptom+", GP = lab_2),
      data.table::data.table(p_pacf_1$total$All$Flow[[key]], ACF = "All", GP = lab_1),
      data.table::data.table(p_pacf_2$total$All$Flow[[key]], ACF = "All", GP = lab_2)
    ))
  }))
  flows$ACF <- factor(flows$ACF, levels = c("None", "Smear+", "Symptom+", "All"))
  
  gs$g_occur <- flows %>%
    mutate(name = ifelse(name == "ACF", "Found", name)) %>%
    select(t, name, mean, ACF, GP) %>%
    pivot_wider(names_from = name, values_from = mean) %>%
    mutate(Notification = Noti + Found) %>%
    rename(Incidence = Act, `Untreated mortality`= Mor, `Care Sought` = Noti) %>%
    select(-Found) %>% 
    pivot_longer(c(Incidence, `Untreated mortality`, Notification, `Care Sought`), values_to = "m") %>%
    mutate(name = factor(name, levels = c("Incidence", "Care Sought", "Notification", "Untreated mortality"))) %>%
    ggplot() +
    geom_line(aes(x = t, y = m, colour = name)) +
    scale_color_discrete("Index") +
    scale_y_continuous("Event rate, per 100 000 capita-year", labels = function(x) x * 1E5) +
    scale_x_continuous("Year", labels = function(x) x + 2020, breaks = c(2020, 2025, 2030)) + 
    facet_grid(GP~ACF, scale = "free_y") +
    theme(legend.position = "bottom", axis.text.x.bottom = element_text(angle = 45, hjust = 1))
  
  
  
  avers <- data.table::rbindlist(lapply(c("Mor", "Act"), function(key) {
    data.table::rbindlist(list(
      data.table::data.table(p_pacf_1$total$d10[[key]], ACF = "Smear+", GP = lab_1),
      data.table::data.table(p_pacf_2$total$d10[[key]], ACF = "Smear+", GP = lab_2),
      data.table::data.table(p_pacf_1$total$d20[[key]], ACF = "Symptom+", GP = lab_1),
      data.table::data.table(p_pacf_2$total$d20[[key]], ACF = "Symptom+", GP = lab_2),
      data.table::data.table(p_pacf_1$total$d30[[key]], ACF = "All", GP = lab_1),
      data.table::data.table(p_pacf_2$total$d30[[key]], ACF = "All", GP = lab_2)
    ))
  }))
  avers$ACF <- factor(avers$ACF, levels = c("Smear+", "Symptom+", "All"))
  
  
  gs$g_avert <- avers %>% 
    ggplot(aes(x = t)) +
    geom_line(aes(y = mean, colour = ACF)) +
    geom_ribbon(aes(ymin = l, ymax = u, fill = ACF), alpha = 0.2) +
    scale_y_continuous("Averted events, thousands", labels = function(x) x/1E3) +
    scale_x_continuous("Year", labels = function(x) x + 2020, breaks = c(2020, 2025, 2030)) +
    facet_grid(name~GP, scale = "free_y", 
               labeller = labeller(name=c(Act="Incident cases", Mor="Untreated deaths"))) +
    expand_limits(y = 0)
  
  
  gs
}

