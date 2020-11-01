

get_transition_matrix <- function(post) {
  
  exo <- post$Total$Exo
  trmat <- array(0, c(6, 6, length(post_a_sn_sp$Total$MC$r_sym)))
  
  for (i in 1:length(post_a_sn_sp$Total$MC$r_sym)) {
  
    pars <- c(exo, lapply(post_a_sn_sp$Total$MC, function(p) p[i]))
    pars
  
  

  
    trmat[, , i] <- with(pars, {
      mat <- rbind(
        c(0, r_sym * (1 - p_sp), r_sym * p_sp, r_sc, 0, 0),
        c(0, 0, r_tr, r_sc, r_death_sn, 1/del_sn),
        c(0, 0, 0, r_sc, r_death_sp, 1/del_sp),
        c(0, 0, 0, 1, 0, 0),
        c(0, 0, 0, 0, 1, 0),
        c(0, 0, 0, 0, 0, 1)
      )
      mat / rowSums(mat)
    })
  
  }
  
  
  apply(trmat, c(1, 2), mean)

}


load("output/Post_Kenya.rdata")

get_transition_matrix(post_a_sn_sp)


load("output/Post_Blantyre.rdata")

get_transition_matrix(post_a_sn_sp)


