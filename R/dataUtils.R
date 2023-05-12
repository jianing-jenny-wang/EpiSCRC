
calcBasicSetUp <- function(dtBasicSetUp_one_simu){
  res = list()
  res$names_of_k <- dtBasicSetUp_one_simu$names_of_k
  res$W_val <- dtBasicSetUp_one_simu$W
  # res$D_w_val <- dtBasicSetUp_one_simu$D_w
  res$K <- dtBasicSetUp_one_simu$K # number of subregions
  res$J <- dtBasicSetUp_one_simu$J # number of lists
  res$D <- dtBasicSetUp_one_simu$D # number of latent factors
  res$P <- dtBasicSetUp_one_simu$P # number of normal people by region
  res$true_logitlisteffect <- dtBasicSetUp_one_simu$list_effects
  return(res)
}


calSimulationParams <- function(dtYfullTable_one_simu, dtNtarget_one_simu, dtyObs_by_loc_one_simu, P, K, J, prev_ub = 0.2){
  single_dt <- dtYfullTable_one_simu # temporarily use one simulated dataset
  # Create the observed dataset for one single simulated dataset
  dtYObs <- subset(dtYfullTable_one_simu, obs_label == 1)

  N_target <- dtNtarget_one_simu # number of target people by region

  N_obs <- dtyObs_by_loc_one_simu$x # number of observed people by region
  TotalN_target <- sum(N_target)
  TotalN_obs <- sum(N_obs)
  Prev <- N_target/P

  #### Data Augmentation ####
  prev_ub <- 0.2
  M_aug <- round(P*prev_ub,digits = 0)
  aug_size <- sum(M_aug) # The augmented size

  pseudo_prev <- N_target/M_aug # Effectively what is the true value of prevalence after DA
  logit_pseudo_prev <- logit(pseudo_prev) # Effectively what is the true value of lambda after DA

  all_zeros_size <- M_aug - N_obs # N of additional zero rows by location
  dt_all_zeros <- list()
  for(i in 1:K){
    dt_all_zeros_loc_i <- array(0, dim = c(all_zeros_size[i],ncol(dtYObs)-1))
    dt_all_zeros_loc_i <- cbind(dt_all_zeros_loc_i, rep(i, nrow(dt_all_zeros_loc_i)))
    dt_all_zeros_loc_i <- as.data.frame(dt_all_zeros_loc_i)
    colnames(dt_all_zeros_loc_i) <- colnames(dtYObs)
    dt_all_zeros[[i]] <- dt_all_zeros_loc_i
  }
  dt_all_zeros <- do.call(rbind,dt_all_zeros)
  # Augmentation to make all locations have untarget + target people
  dtYObs_aug <- rbind(dtYObs, dt_all_zeros)
  # Order by observed first, and then unobserved, nrow = sum(P)
  dtYObs_aug <- arrange(dtYObs_aug, desc(obs_label))
  loc_label = dtYObs_aug$loc_label
  dtYObs_yaug <- dtYObs_aug[,seq(J)] # this is the dataset to be passed into the model
  ########## FINISH DATA MANIPULATION ##########
  Mod_GenMth_ICAR_Data <- list(y = dtYObs_yaug)

  out = list(aug_size=aug_size,
             logit_pseudo_prev = logit_pseudo_prev,
             TotalN_target = TotalN_target,
             TotalN_obs = TotalN_obs,
             Prev = Prev,
             loc_label = loc_label,
             Mod_GenMth_ICAR_Data = Mod_GenMth_ICAR_Data
            )

  return(out)
}

