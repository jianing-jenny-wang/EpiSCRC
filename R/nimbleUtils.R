
prepareNimbleModel <- function(scrNimbleCode, constants, data, inits, monitors = c("N","pesudo_logit_prev", "sigma2.phi", "resPhi", "beta0","loading_mat_upper_tri")) {
  nimbleModelUncompiled <- nimble::nimbleModel(code = scrNimbleCode,
                                               constants = constants,
                                               data = data,
                                               inits = inits)
  # ## (1) Configuration
  modelConfig           <- nimble::configureMCMC(nimbleModelUncompiled,
                                                 monitors = monitors,
                                                 enableWAIC = nimbleOptions('enableWAIC' = TRUE),
                                                 print = TRUE)
  # ## (2) Build MCMC algorithm for the model
  mcmcAlgorithm         <- nimble::buildMCMC(modelConfig)
  # ## (3) Compile C++ after configuration
  compiledNimbleModel   <- nimble::compileNimble(mcmcAlgorithm,
                                                 project = nimbleModelUncompiled)

  return(compiledNimbleModel)
}


getPosteriorByMCMC <- function(compiledNimbleModel, nchains = 3, iter = 68000, nburnin = 50000, nthin = 100, setSeed=9) {
  posteriorSamples <- nimble::runMCMC(compiledNimbleModel,
                                      niter = niter,
                                      nburnin = nburnin,
                                      nchains = nchains,
                                      thin = nthin,
                                      setSeed = setSeed,
                                      progressBar = TRUE,
                                      summary = TRUE,
                                      WAIC = TRUE,
                                      samplesAsCodaMCMC = TRUE)
  return(posteriorSamples)
}


prepareNimbleCode <- function(){
  # "scrNimbleCode = Mod_GenMth_ICAR_Code"
  ### Create a NIMBLE model ###
  ## Define the model code, all nodes, variables and relationship, and all constants
  scrNimbleCode <- nimble::nimbleCode( {
    ############
    ## Priors ##
    ############

    # Spatial RE [ICAR] #
    alpha ~ dflat() # vague uniform prior due to centering constraint at 0

    # variance parameter of the phi #
    ## Option 1: Variance ~ Inverse Gamma, Precision tau2 ~ Gamma
    # tau2 ~ dgamma(shape = 0.5, scale = 1/0.5) # inverse variance of spatial component (precision of the ICAR Component, equivalent variance ~ IG(0.5,0.5))
    ## Option 2: Variance ~ Half-Cauchy (t distribution with df = 1)
    inv_tau2  ~ T(dt(mu = 0, tau = sd_spRE_prior, df = 1), 0, ) # half-t truncated at 0 and set no upper limit
    tau2 <- 1/inv_tau2 # compute precision  phi[1:K] ~ dcar_normal(adj[1:L], weights[1:L], num[1:K], tau2, zero_mean = 1) # L = length of Adj_vec, K = n of regions, zero_mean = 1 put constraint centers at 0

    # phi from ICAR #
    phi[1:K] ~ dcar_normal(adj[1:L], weights[1:L], num[1:K], tau2, zero_mean = 1) # L = length of Adj_vec, K = n of regions, zero_mean = 1 put constraint centers at 0

    # List_Effects FE #
    for(j in 1:J){
      logit_list_effects[j] ~ dnorm(0, sd = sigma_mu)
    }

    # variance parameter of the theta #
    ## Option 1: Variance ~ Inverse Gamma, precision ~ Gamma
    # prec_theta ~ dgamma(shape = 0.5, scale = 1/0.5) # inverse variance of theta prior
    # sigma_theta <- 1/prec_theta
    ## Option 2: Variance ~ Half-Cauchy (t distribution with df = 1)
    sigma_theta ~ T(dt(mu = 0, tau = sd_theta_prior, df = 1), 0, )  # half-t truncated at 0 and set no upper limit

    # Factor Loading Matrix - alpha FE#
    for(r in 1:J){
      for(c in 1:D){
        loading_mat[r,c] ~ dnorm(0, sd = sigma_loading)
        loading_mat_upper_tri[r,c] <- loading_mat[r,c]*mask_loading_mat[r,c]
      }
    }

    # Latent Factors - theta covariance matrix#
    # if(D > 1){
    #   # LKJ prior for covariance matrix of theta #
    #   # https://r-nimble.org/html_manual/cha-writing-models.html #
    #   # Ustar[1:D,1:D] ~ dlkj_corr_cholesky(eta = 1, p = D) # eta = 1, assign uniform distribution (non-informative); eta > 1 less correlation (~ identity matrix); eta < 1 stronger correlation
    #   # U[1:D,1:D] <- uppertri_mult_diag(Ustar[1:D, 1:D], theta_sds[1:D])
    #   # Inverse-Wishart prior for covariance matrix of theta #
    #   cov_theta ~ dinvwish(S = S, df = D)
    # }


    ################
    ## Likelihood ##
    ################

    # Prevalence Level #
    # for(k in 1:K){
    #  eps[k] ~ dnorm(0, sd = sigma_eps) # location wise dispersion error term
    # }

    for (i in 1:M){ # M = augmented large population

      logit(lambda[i]) <- alpha + phi[prev_re[i]] # pseudo logit prev via nested indicators, dim = length of M

      z[i] ~ dbern(lambda[i]) # inclusion indicator to be a target subject, same location's subjects share the same pseudo prev

      # if(D == 1){
      theta[i] ~ dnorm(mu_theta, sd = sigma_theta) # vector of individual heterogeneity effect is a random (latent) effect
      loading_latent_factor[i,1:J] <- theta[i] * t(loading_mat_upper_tri[1:J,1]) # dim = 1 * trans(J*D)
      # }
      # if(D > 1){
      #   # option1: LKJ prior for covariance matrix of theta #
      #   # t(theta[1:M, 1:D]) ~ dmnorm(mu_theta[1:D], cholesky = U[1:D, 1:D], prec_param = 0)
      #
      #   # option 2: Inverse-Wishart prior for covariance matrix of theta #
      #   # t(theta[1:M, 1:D]) ~  dmnorm(mean = mu_theta, cholesky = chol(cov_theta), prec_param = 0)

      #   loading_latent_factor[i,1:J] <- t(loading_mat_upper_tri[1:J,1:D]%*%t(theta[i,1:D])) # trans(J*D * D*1)
      # }

      for(j in 1:J){
        # Detection prob model #
        logit(p.detect[i,j]) <- logit_list_effects[j] + loading_latent_factor[i,j] # detection probability

        # Compute effective p #
        p.eff[i,j] <- z[i] * p.detect[i,j]  # being a potential target * being observed
        y[i,j] ~ dbern(p.eff[i,j]) # individual capture histories
      } # j
    } # i

    #########################
    ## Derive the Quantity ##
    #########################
    # <1> Number of target people
    N <- sum(z[1:M])
    # <2> Overall baseline value of lambda (to compute baseline prevalence)
    beta0 <- alpha # compute overall baseline value of lambda
    # <3> Variance of spatial component
    # sigma2.phi <- 1/tau2 # variance of spatial component using inverse gamma for variance parameter
    sigma2.phi <- inv_tau2 # variance of spatial component using half-cauchy for variance parameter
    # <4> Area-specific residual
    resPhi[1:K] <- phi[1:K]
    # <5> Overall residuals at individual level
    # res[1:M] <- (z[1:M]-lambda[1:M])/sqrt(lambda[1:M])
    # <6> Lambda
    for(k in 1:K){
      pesudo_logit_prev[k] <- alpha + phi[k]
    }

  })
  ########## FINISH BAYESIAN MODEL WRITTING ##########
  return(scrNimbleCode)

}


