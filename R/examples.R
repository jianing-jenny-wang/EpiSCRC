

path_to_input <- "/projectnb2/jw-thesis/Proj2-CRC-SAE/Simulation/Simu5/Simu5_GenDt/Simu5_2_OHcty_woNug/Simu5.2.3_J5VaryNoExtremePDetect/"
task_id <- "task_id"
dtBasicSetUp_one_simu = paste0(path_to_input, "BasicSetUp.RData")
dtYfullTable_one_simu = paste0(path_to_input, "dtYfullTable_simu", task_id, ".RData")
dtNtarget_one_simu = paste0(path_to_input, "dtNtarget_simu", task_id, ".RData")
dtyObs_by_loc_one_simu = paste0(path_to_input, "dtyObs_by_loc_simu", task_id, ".RData")

# experiment = fit(dtBasicSetUp_one_simu,
#                  dtYfullTable_one_simu,
#                  dtNtarget_one_simu,
#                  dtyObs_by_loc_one_simu,
#                  alpha = c(0.05, 0.03, 0.06),
#                  tau2 = c(0.5, 0.1, 0.2),
#                  prec_theta = c(10, 20, 50),
#                  prev_ub = 0.2,
#                  sigma_mu = 10,
#                  sigma_loading = 10,
#                  SD_HalfCauchy_Prior = 10,
#                  seedLoadingMat = c(123, 234, 345),
#                  seedTheta = c(123, 234, 345),
#                  seedCovTheta = NaN,
#                  monitors = c("N","pesudo_logit_prev", "sigma2.phi", "resPhi", "beta0","loading_mat_upper_tri"),
#                  nchain = 3,
#                  iter = 68000,
#                  nburnin = 50000,
#                  nthin = 100,
#                  setSeed = 9)
