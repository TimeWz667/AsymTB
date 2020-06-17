
eval_p_pacf <- function(fitted, pop_foreword, sen = 0.1, n_foreward = 15, n_acf = 3, n_iter = 300) {
  
  an_acf <- function(y0, pars, sens = c(0, 0, 1)) {
    found <- y0 * c(sens, 0, 0, 0, 0)
    keep <- y0 - found
    # names(found) <- paste0(names(found), 0)
    names(keep) <- paste0(names(keep), 0)
    
    list(
      found = found,
      keep = c(pars, keep)
    )
  }
  
  
  find_ss <- function(pars) {
    with(as.list(pars), {
      ra <- r_sc
      rn <- r_sc + r_death_sn + 1 / del_sn - adr
      rp <- r_sc + r_death_sp + 1 / del_sp - adr
      
      sn0 <- (1 - p_sp) * r_sym / (r_tr + rn)    
      sp0 <- p_sp * r_sym / (rp) + r_tr * sn0 / (rp)
      
      a <- 1 / (1 + sn0 + sp0)
      
      c(A0 = a, Sn0 = sn0 * a, Sp0 = sp0 * a)
    })
  }
  
  
  do_pacf <- function(pars, cm, prv0, pop, n_acf, n_foreward, sens = c(0, 0, 1)) {
    inp <- c(pars, find_ss(pars), c(Cured0 = 0, Death0 = 0, NotiSn0 = 0, NotiSp0 = 0))
    
    cm$set_user(user = inp)
    
    # move to the first acf timing
    t0 <- 0
    ys <- cm$run(seq(t0, t0 + 0.5, 1/10))
    y0 <- ys[nrow(ys), 2:8]
    found <- c()
    t0 <- 0.5
    while(t0 < n_acf) {
      yield <- an_acf(y0, pars, sens)
      found <- rbind(found, c(t = t0, yield$found))
      cm$set_user(user = yield$keep)
      ys <- rbind(ys, cm$run(seq(t0, t0 + 0.5, 1/10))[-1, ])
      y0 <- ys[nrow(ys), 2:8]
      t0 <- t0 + 0.5
    }
    yield <- an_acf(y0, pars, sens)
    found <- rbind(found, c(t = t0, yield$found))
    cm$set_user(user = yield$keep)
    
    ys <- rbind(ys, cm$run(seq(t0, n_foreward, 1/10))[-1, ])
    
    found <- found[, 1:4]
    colnames(found) <- c("t", "ACF_A", "ACF_Sn", "ACF_Sp")
    found <- cbind(found, ACF = rowSums(found[, -1]))
    
    found_ts <- apply(found / 0.5, 2, rep, each = 5)
    found_ts[, "t"] <- ys[1 + 1:(10 * n_acf), "t"]
    
    
    ys <- merge(ys, found_ts, all.x = T)
    ys[is.na(ys)] <- 0
    
    cols <- c("A", "Sn", "Sp", "Cured", "Death", "NotiSn", "NotiSp", 
              "Act", "Sym", "Noti", "Mor", "ACF", "ACF_A", "ACF_Sn", "ACF_Sp")
    
    ys[, cols] <- ys[, cols] * prv0
    ysn <- ys
    ysn[, cols] <- ysn[, cols] * pop
    
    aoc <- cbind(
      t = ysn[-1, "t"],
      Act = cumsum((ysn[-1, "Act"] + ysn[-nrow(ysn), "Act"]) / 20),
      Sym = cumsum((ysn[-1, "Sym"] + ysn[-nrow(ysn), "Sym"]) / 20),
      Mor = cumsum((ysn[-1, "Mor"] + ysn[-nrow(ysn), "Mor"]) / 20),
      A = cumsum((ysn[-1, "A"] + ysn[-nrow(ysn), "A"]) / 20),
      Sn = cumsum((ysn[-1, "Sn"] + ysn[-nrow(ysn), "Sn"]) / 20),
      Sp = cumsum((ysn[-1, "Sp"] + ysn[-nrow(ysn), "Sp"]) / 20)
    )
    aoc <- rbind(c(t = 0, Act = 0, Sym = 0, Mor = 0, A = 0, Sn = 0, Sp = 0), aoc)
    
    list(Ys = as.matrix(ys), YsN =as.matrix(ysn), Found = found, AOC = aoc)
  }
  
  do_sel_pacf <- function(sel, exo, prv0, pop, n_acf, n_foreward, sens = c(0, 0, 1)) {
    temp <- lapply(sel, function(pars) {
      do_pacf(c(pars, exo), cm, prv0, pop, n_acf, n_foreward, sens)
    })
    
    yss <- array(0, c(dim(temp[[1]]$Ys), n_iter))
    ysns <- array(0, c(dim(temp[[1]]$YsN), n_iter))
    founds <- array(0, c(dim(temp[[1]]$Found), n_iter))
    aocs <- array(0, c(dim(temp[[1]]$AOC), n_iter))
    for (j in 1:n_iter) {
      yss[, , j] <- temp[[j]]$Ys
      ysns[, , j] <- temp[[j]]$YsN
      founds[, , j] <- temp[[j]]$Found
      aocs[, , j] <- temp[[j]]$AOC
    }
    
    dimnames(yss)[[2]] <- colnames(temp[[1]]$Ys)
    dimnames(ysns)[[2]] <- colnames(temp[[1]]$YsN)
    dimnames(founds)[[2]] <- colnames(temp[[1]]$Found)
    dimnames(aocs)[[2]] <- colnames(temp[[1]]$AOC)

    list(
      Ys = yss, 
      YsN = ysns,
      Found = founds,
      AOC = aocs
    )
  }
  
  sum_array <- function(arr, key, new_name = key) {
    ti <- arr[, "t", 1]
    arr <- arr[, key, ]
    
    data.frame(
      t = ti,
      name = new_name,
      mean = apply(arr, 1, mean),
      l = apply(arr, 1, quantile, p = 0.025),
      m = apply(arr, 1, quantile, p = 0.5),
      u = apply(arr, 1, quantile, p = 0.975)
    )
  }
  
  sum_avert_array <- function(arr1, arr0, key, new_name = key) {
    arr <- (arr0 - arr1)[, key, ]
    
    data.frame(
      t = arr0[, "t", 1],
      name = new_name,
      mean = apply(arr, 1, mean),
      l = apply(arr, 1, quantile, p = 0.025),
      m = apply(arr, 1, quantile, p = 0.5),
      u = apply(arr, 1, quantile, p = 0.975)
    )
  }
  
  diff_pacf <- function(res1, res0) {
    list(
      Act = sum_avert_array(res1$AOC, res0$AOC, "Act"),
      Sym = sum_avert_array(res1$AOC, res0$AOC, "Sym"),
      Mor = sum_avert_array(res1$AOC, res0$AOC, "Mor"),
      PrvA = sum_avert_array(res1$AOC, res0$AOC, "A"),
      PrvSn = sum_avert_array(res1$AOC, res0$AOC, "Sn"),
      PrvSp = sum_avert_array(res1$AOC, res0$AOC, "Sp")
    )
  }
  
  sum_pacf <- function(res) {
    list(
      Flow = list(
        Act = sum_array(res$Ys, "Act"),
        Sym = sum_array(res$Ys, "Sym"),
        SC = sum_array(res$Ys, "SC"),
        AC = sum_array(res$Ys, "AC"),
        Noti = sum_array(res$Ys, "Noti"),
        Mor = sum_array(res$Ys, "Mor"),
        ACF = sum_array(res$Ys, "ACF")
      ),
      AOC = list(
        Act = sum_array(res$AOC, "Act"),
        Sym = sum_array(res$AOC, "Sym"),
        Mor = sum_array(res$AOC, "Mor"),
        A = sum_array(res$AOC, "A"),
        Sn = sum_array(res$AOC, "Sn"),
        Sp = sum_array(res$AOC, "Sp"),
        NotiSn = sum_array(res$YsN, "NotiSn"),
        NotiSp = sum_array(res$YsN, "NotiSp")
      )
    )
   
  }
  
  ## Main -----
  f <- "odin/ode_anp_ts.R"
  model <- odin::odin(f)
  cm <- model()
  
  exo <- fitted$female$Exo
  fitted$female$Pop <- pop_foreword$pop_f
  fitted$male$Pop <- pop_foreword$pop_m
  
  res <- lapply(fitted, function(x) {
    prv0 <- x$Prevalence[nrow(x$Prevalence), 1]
    mc <- as.data.frame(x$MC)[c("del_sn", "del_sp", "r_tr", "p_sp", "adr")]
    pop <- x$Pop[1:n_foreward, "Pop"]
    pop <- c(pop[1], rep(pop, each = n_foreward))
    
    if(nrow(mc) < n_iter) {
      sel <- sample.int(nrow(mc), n_iter, replace = T)
    } else {
      sel <- sample.int(nrow(mc), n_iter)
    }
    
    sel <- lapply(sel, function(i) {
      pars <- as.list(mc[i, ])
      pars
    })
    
    temp <- list(
      None = do_sel_pacf(sel, exo, prv0, pop, n_acf, n_foreward, c(0, 0, 0)),
      Sp = do_sel_pacf(sel, exo, prv0, pop, n_acf, n_foreward, c(0, 0, sen)),
      Sym = do_sel_pacf(sel, exo, prv0, pop, n_acf, n_foreward, c(0, sen, sen)),
      All = do_sel_pacf(sel, exo, prv0, pop, n_acf, n_foreward, c(sen, sen, sen))
    )
    
    res <- lapply(temp, sum_pacf)
    
    res$d10 <- diff_pacf(temp$Sp, temp$None)
    res$d20 <- diff_pacf(temp$Sym, temp$None)
    res$d30 <- diff_pacf(temp$All, temp$None)
    res$d21 <- diff_pacf(temp$Sym, temp$Sp)
    res$d31 <- diff_pacf(temp$All, temp$Sp)
    res$d32 <- diff_pacf(temp$All, temp$Sym)
    
    res$raw <- temp
    res$pop <- pop
    res
  })
  
  ns  <- array(res$female$pop + res$male$pop, dim(res$male$raw$None$Ys)) 
  wf <- array(res$female$pop, dim(res$male$raw$None$Ys)) / ns
  wm <- array(res$male$pop, dim(res$male$raw$None$Ys)) / ns
  
  temp <- 
    list(
    None = list(Ys = res$female$raw$None$Ys * wf + res$male$raw$None$Ys * wm, 
                YsN = res$female$raw$None$YsN + res$male$raw$None$YsN, 
                AOC = res$female$raw$None$AOC + res$male$raw$None$AOC),
    Sp = list(Ys = res$female$raw$Sp$Ys * wf + res$male$raw$Sp$Ys * wm, 
              YsN = res$female$raw$Sp$YsN + res$male$raw$Sp$YsN, 
              AOC = res$female$raw$Sp$AOC + res$male$raw$Sp$AOC),
    Sym = list(Ys = res$female$raw$Sym$Ys * wf + res$male$raw$Sym$Ys * wm, 
               YsN = res$female$raw$Sym$YsN + res$male$raw$Sym$YsN, 
               AOC = res$female$raw$Sym$AOC + res$male$raw$Sym$AOC),
    All = list(Ys = res$female$raw$All$Ys * wf + res$male$raw$All$Ys * wm, 
               YsN = res$female$raw$All$YsN + res$male$raw$All$YsN, 
               AOC = res$female$raw$All$AOC + res$male$raw$All$AOC) 
  )
  
  total <- lapply(temp, sum_pacf)
  
  total$d10 <- diff_pacf(temp$Sp, temp$None)
  total$d20 <- diff_pacf(temp$Sym, temp$None)
  total$d30 <- diff_pacf(temp$All, temp$None)
  total$d21 <- diff_pacf(temp$Sym, temp$Sp)
  total$d31 <- diff_pacf(temp$All, temp$Sp)
  total$d32 <- diff_pacf(temp$All, temp$Sym)
  
  res$total <- total
  
  res$female$raw <- NULL
  res$male$raw <- NULL
  
  res
}

