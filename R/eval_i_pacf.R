
eval_i_pacf <- function(fitted, dur_asym = 3, beta = c(0.22, 0.22, 1), n_iter = 300) {
  
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
  
  
  do_pacf <- function(pars, cm, start_mo = 3, sens = c(0, 0, 1), beta = c(0.22, 0.22, 1)) {
    inp <- c(pars, c(A0 = 1, Sn0 = 0, Sp0 = 0, Cured0 = 0, Death0 = 0, NotiSn0 = 0, NotiSp0 = 0))
    
    cm$set_user(user = inp)
    
    # move to the first acf timing
    if (start_mo > 0) {
      ys <- cm$run(seq(0, start_mo / 12, 1/12))
      y0 <- ys[nrow(ys), -1]
    } else {
      ys <- cm$run(seq(0, 10, 1))
      y0 <- ys[1, -1]
      ys <- t(as.matrix(ys[1, ]))
    }
    
    found <- c()
    while (sum(y0[1:3]) > 1E-5 | nrow(ys) < 100) {
      yield <- an_acf(y0, pars, sens)
      
      found <- rbind(found, yield$found)
      
      cm$set_user(user = yield$keep)
      ys <- rbind(ys, cm$run(seq(0, 0.5, 1/12))[-1, ])
      y0 <- ys[nrow(ys), -1]
    }
    
    found <- colSums(found)[1:3]
    names(found) <- paste0("ACF_", names(found))
    
    
    duration <- c(
      Dur_A = sum((ys[-1, 2] + ys[-nrow(ys), 2]) / 24),
      Dur_Sn = sum((ys[-1, 3] + ys[-nrow(ys), 3]) / 24),
      Dur_Sp = sum((ys[-1, 4] + ys[-nrow(ys), 4]) / 24)
    )
    
    inf <- duration * beta
    names(inf) <- c("Inf_A", "Inf_Sn", "Inf_Sp")
    
    res <- c(ys[nrow(ys), c("Cured", "Death", "NotiSp", "NotiSn")], found, 
             duration, Dur=sum(duration), 
             inf, "Inf"=sum(inf))
    list(Ys = ys, Indices = res)
  }
  
  sum_idx <- function(idx) {
    data.frame(
      name = rownames(idx),
      mean = rowMeans(idx, na.rm = T),
      l = apply(idx, 1, quantile, p = 0.025, na.rm = T),
      m = apply(idx, 1, median, na.rm = T),
      u = apply(idx, 1, quantile, p = 0.975, na.rm = T),
      row.names = NULL
    )
  }
  
  do_sel_pacf <- function(sel, sens = c(0, 0, 1), beta = c(0.22, 0.22, 1)) {
    temp <- lapply(sel, function(pars) {
      temp <- lapply(0:6, function(mo) do_pacf(pars, cm, start_mo = mo, sens, beta))
      
      yss <- array(0, c(61, 8, 7))
      for (j in 1:7) {
        yss[, , j] <- temp[[j]]$Ys[1:61, ]
      }
      yss <- apply(yss, c(1, 2), mean, na.rm = T)
      colnames(yss) <- colnames(temp[[1]]$Ys)
      
      list(
        Ys = yss, 
        Indices = rowMeans(sapply(temp, function(y) y$Indices), na.rm = T)
      )
    })
    
    yss <- array(0, c(61, 8, n_iter))
    for (j in 1:n_iter) {
      yss[, , j] <- temp[[j]]$Ys[1:61, ]
    }
    yss <- apply(yss, c(1, 2), mean, na.rm = T)
    colnames(yss) <- colnames(temp[[1]]$Ys)
    
    idx <- sapply(temp, function(y) y$Indices)
    list(
      Ys = yss, 
      Indices = sum_idx(idx),
      IndicesRaw = idx
    )
  }

  sum_diff <- function(res1, res0, res_base = res0) {
    dif <- res0$IndicesRaw - res1$IndicesRaw
    red <- dif / res_base$IndicesRaw
    
    list(
      Difference = sum_idx(dif),
      Reduction = sum_idx(dif / res0$IndicesRaw),
      IncreRed = sum_idx(dif / res_base$IndicesRaw)
    )
  }
  
  ## Main -----
  f <- "odin/ode_anp_acf_cohort.R"
  model <- odin::odin(f)
  cm <- model()

  res <- lapply(fitted, function(x) {
    mc <- as.data.frame(x$MC)[c("del_sn", "del_sp", "r_tr", "p_sp")]

    if (nrow(mc) < n_iter) {
      sel <- sample.int(nrow(mc), n_iter, replace = T)
    } else {
      sel <- sample.int(nrow(mc), n_iter)
    }
    
    sel <- lapply(sel, function(i) {
      pars <- as.list(mc[i, ])
      pars$r_sym <- 12 / dur_asym
      pars
    })
    
    temp <- list(
      None = do_sel_pacf(sel, c(0, 0, 0), beta),
      Sp = do_sel_pacf(sel, c(0, 0, 1), beta),
      Sym = do_sel_pacf(sel, c(0, 1, 1), beta),
      All = do_sel_pacf(sel, c(1, 1, 1), beta)
    )
    
    temp$d10 <- sum_diff(temp$Sp, temp$None)
    temp$d20 <- sum_diff(temp$Sym, temp$None)
    temp$d30 <- sum_diff(temp$All, temp$None)
    temp$d21 <- sum_diff(temp$Sym, temp$Sp, temp$None)
    temp$d31 <- sum_diff(temp$All, temp$Sp, temp$None)
    temp$d32 <- sum_diff(temp$All, temp$Sym, temp$None)
    temp
  })
  res
}
