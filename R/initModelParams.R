genLoading_mat_init <- function(J, D, seed=NaN){
  ## Create initial values for loading factor matrix ##
  if(!is.nan(seed)){
    set.seed(seed)
  }
  loading_mat_init <- matrix(rnorm(n = J*D, mean = 0, sd = 1), nrow = J, ncol = D)
  loading_mat_init <- round(loading_mat_init, digits = 3)
  loading_mat_init[upper.tri(loading_mat_init)] <- 0
  return(loading_mat_init)
}

genTheta_init <- function(aug_size, D, seed=NaN){
  ## Create initial values for individual latent factors ##
  if(!is.nan(seed)){
    set.seed(seed)
  }

  thetas = rnorm(n = aug_size*D, mean = 0, sd = 0.1)
  Theta = matrix(thetas, nrow = aug_size, ncol = D )
  return(Theta)
}


genCovTheta_init <- function(D, seed = NaN){
  if(!is.nan(seed)){
    set.seed(seed)
  }
  corrs = matrixsampling::rinvwishart(n = 1, nu = D, Omega = diag(D))
  cov_theta = matrix(corrs, nrow = D, ncol = D)
  return(cov_theta)
}

putInitPerChain <- function(alpha, tau2, prec_theta, J, D, K, aug_size, loading_mat, theta, cov_theta){
  init = list(alpha = alpha,
              tau2 = tau2,  # if use invgamma for spRE precision parameter
              inv_tau2 = 1.0/tau2,  # if use half-cauchy for spRE variance parameter
              logit_list_effects = rep(1, J),
              loading_mat = loading_mat,
              prec_theta = prec_theta, # if use invgamma for theta precision parameter
              sigma_theta = 1.0/prec_theta, # if use half-cauchy for theta variance parameter
              cov_theta = cov_theta,
              theta = theta,
              phi = rep(0.1, K),
              z = rep(1, aug_size)
  )
  return(init)
}

getModelInits <- function(nchain, alpha, tau2, prec_theta, J, D, K, aug_size, seedLoadingMat=NaN, seedTheta=NaN, seedCovTheta=NaN){
  # inits = Mod_GenMth_ICAR_Inits
  if(length(seedLoadingMat)<nchain){
    seedLoadingMat = rep(NaN, nchain)
  }
  if(length(seedTheta)<nchain){
    seedTheta = rep(NaN, nchain)
  }
  if(length(seedCovTheta)<nchain){
    seedCovTheta = rep(NaN, nchain)
  }

  loading_mat = lapply(seedLoadingMat, FUN = function(x) genLoading_mat_init(J, D, seed=x))
  theta       = lapply(seedTheta,      FUN = function(x) genTheta_init(aug_size, D, seed=x))
  cov_theta   = lapply(seedCovTheta,   FUN = function(x) genCovTheta_init(D, seed=x))

  inits = mapply(putInitPerChain, alpha, tau2, prec_theta, J, D, K, aug_size, loading_mat, theta, cov_theta, SIMPLIFY = FALSE)
  return(inits)
}
