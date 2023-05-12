#' Fit end-to-end ICAR model
#'
#' @param dtBasicSetUp_one_simu
#' Example of data - writing some descriptions here
#' @param dtYfullTable_one_simu
#' @param dtNtarget_one_simu
#' @param dtyObs_by_loc_one_simu
#' @param tau_prior
#' @param sigma_theta_prior
#' @param theta_prior
#' @param alpha
#' @param tau2
#' @param prec_theta
#' @param prev_ub
#' @param sigma_mu
#' @param sigma_loading
#' @param SD_HalfCauchy_Prior
#' @param seedLoadingMat
#' @param seedTheta
#' @param seedCovTheta
#' @param monitors
#' @param nchain
#' @param iter
#' @param nburnin
#' @param nthin
#' @param setSeed
#' @return
#' Posterior samples of ICAR model paramters
#' @export
#'
#' @examples
#' fit(dtBasicSetUp_one_simu,
#'     dtYfullTable_one_simu,
#'     dtNtarget_one_simu,
#'     dtyObs_by_loc_one_simu,
#'     tau_prior = c(inv_tau2  ~ T(dt(mu = 0, tau = sd_spRE_prior, df = 1), 0, ),"tau2 <- 1/inv_tau2"),
#'     igma_theta_prior = c(sigma_theta ~ T(dt(mu = 0, tau = sd_theta_prior, df = 1), 0, )),
#'     theta_prior = c(theta[i] ~ dnorm(mu_theta, sd = sigma_theta),"loading_latent_factor[i,1:J] <- theta[i] * t(loading_mat_upper_tri[1:J,1])"),
#'     alpha = c(0.05, 0.03, 0.06),
#'     tau2 = c(0.5, 0.1, 0.2),
#'     prec_theta = c(10, 20, 50),
#'     prev_ub = 0.2,
#'     sigma_mu = 10,
#'     sigma_loading = 10,
#'     SD_HalfCauchy_Prior = 10,
#'     seedLoadingMat = c(123, 234, 345),
#'     seedTheta = c(123, 234, 345),
#'     seedCovTheta=NaN,
#'     monitors = c("N","pesudo_logit_prev", "sigma2.phi", "resPhi", "beta0","loading_mat_upper_tri"),
#'     nchain = 3,
#'     iter = 68000,
#'     nburnin = 50000,
#'     nthin = 100,
#'     setSeed=9)
fit <- function(dtBasicSetUp_one_simu, dtYfullTable_one_simu, dtNtarget_one_simu, dtyObs_by_loc_one_simu, tau_prior, sigma_theta_prior, theta_prior, alpha, tau2, prec_theta, prev_ub = 0.2, sigma_mu = 10, sigma_loading = 10, SD_HalfCauchy_Prior = 10, seedLoadingMat=NaN, seedTheta=NaN, seedCovTheta=NaN, monitors = c("N","pesudo_logit_prev", "sigma2.phi", "resPhi", "beta0","loading_mat_upper_tri"), nchain = 3, iter = 68000, nburnin = 50000, nthin = 100, setSeed=FALSE){

  basicSetup <- calcBasicSetUp(dtBasicSetUp_one_simu)
  P <- basicSetup$P
  K <- basicSetup$K

  J <- basicSetup$J
  D <- basicSetup$D

  names_of_k <- res$names_of_k
  W_val      <- res$W_val

  mask_loading_factor_mat <- getMask(J, D)

  modelConst  <- prepareICARmodelConst(names_of_k, W_val)
  Adj_vec     <- modelConst$Adj_vec
  Weights_vec <- modelConst$Weights_vec
  D_w_ls      <- modelConst$D_w_ls


  simParams <- calSimulationParams(dtYfullTable_one_simu, dtNtarget_one_simu, dtyObs_by_loc_one_simu, P, K, J, prev_ub = prev_ub)
  aug_size             <- simParams$aug_size
  loc_label            <- simParams$loc_label

  Mod_GenMth_ICAR_Data <- simParams$Mod_GenMth_ICAR_Data
  Mod_GenMth_ICAR_Consts <- getModelConsts(aug_size, J, K, D, Adj_vec, Weights_vec, D_w_ls, SD_HalfCauchy_Prior, mask_loading_factor_mat, loc_label, sigma_mu = sigma_mu, sigma_loading = sigma_loading)
  Mod_GenMth_ICAR_Inits <- getModelInits(nchain, alpha, tau2, prec_theta, J, D, K, aug_size, seedLoadingMat=seedLoadingMat, seedTheta=seedTheta, seedCovTheta=seedCovTheta)
  Mod_GenMth_ICAR_Code <- hyper2Likelihood(tau_prior, sigma_theta_prior, theta_prior)


  C_mcmc_Mod_GenMth_ICAR <- prepareNimbleModel(scrNimbleCode = Mod_GenMth_ICAR_Code,
                                               constants = Mod_GenMth_ICAR_Consts,
                                               data = Mod_GenMth_ICAR_Data,
                                               inits = Mod_GenMth_ICAR_Inits,
                                               monitors=monitors)

  samples_Mod_GenMth_ICAR <- getPosteriorByMCMC(compiledNimbleModel = C_mcmc_Mod_GenMth_ICAR,
                                                nchains = nchains,
                                                iter = iter,
                                                nburnin = nburnin,
                                                nthin = nthin,
                                                setSeed=setSeed)

  return(samples_Mod_GenMth_ICAR)
}
