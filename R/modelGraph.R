
tau_prior = c(inv_tau2  ~ T(dt(mu = 0, tau = sd_spRE_prior, df = 1), 0, ),
              "tau2 <- 1/inv_tau2"
              )
sigma_theta_prior = c(sigma_theta ~ T(dt(mu = 0, tau = sd_theta_prior, df = 1), 0, ))

theta_prior = c(theta[i] ~ dnorm(mu_theta, sd = sigma_theta),# vector of individual heterogeneity effect is a random (latent) effect
                "loading_latent_factor[i,1:J] <- theta[i] * t(loading_mat_upper_tri[1:J,1])"# dim = 1 * trans(J*D)
                )

helper <- function(codeline){
  if(is.character(codeline)){
    codeline <- parse(text=codeline)
  }
  codeline <- as.expression(codeline)
  return(codeline)
}


list2expr <- function(codeList){
  parsedList <- lapply(codeList, helper)
  output_call <- as.call(c(as.symbol("{"), unlist(parsedList)))
  return(output_call)
}


listLikelihoodIndividual <- function(theta_prior){
  likelihood_i = c(
    "logit(lambda[i]) <- alpha + phi[prev_re[i]]", # pseudo logit prev via nested indicators, dim = length of M
    z[i] ~ dbern(lambda[i]), # inclusion indicator to be a target subject, same location's subjects share the same pseudo prev

    #Theta prior
    theta_prior,

    expression(for(j in 1:J){
      # Detection prob model #
      logit(p.detect[i,j]) <- logit_list_effects[j] + loading_latent_factor[i,j] # detection probability

      # Compute effective p #
      p.eff[i,j] <- z[i] * p.detect[i,j]  # being a potential target * being observed
      y[i,j] ~ dbern(p.eff[i,j]) # individual capture histories
    } # j
    )
  )
  return(likelihood_i)
}

exprLikelihoodAllData <- function(expr_likelihood_i){
  parsedLikelihood <- as.expression(
    as.call(
      c(as.symbol("for"),
        as.symbol("i"),
        as.call(c(as.symbol(":"), 1, as.symbol("M"))),
        expr_likelihood_i
      )
    )
  )

  return(parsedLikelihood)
}


listFullLikelihoodwithHyper <- function(tau_prior, sigma_theta_prior, parsedLikelihood){
  model <- c(
    # Spatial RE [ICAR] #
    alpha ~ dflat(), # vague uniform prior due to centering constraint at 0

    # variance parameter of the phi #
    tau_prior,

    # phi from ICAR #
    phi[1:K] ~ dcar_normal(adj[1:L], weights[1:L], num[1:K], tau2, zero_mean = 1), # L = length of Adj_vec, K = n of regions, zero_mean = 1 put constraint centers at 0

    # List_Effects FE #
    expression(for(j in 1:J){
      logit_list_effects[j] ~ dnorm(0, sd = sigma_mu)
    }),

    # variance parameter of the theta #
    sigma_theta_prior,

    # Factor Loading Matrix - alpha FE#
    expression(for(r in 1:J){
      for(c in 1:D){
        loading_mat[r,c] ~ dnorm(0, sd = sigma_loading)
        loading_mat_upper_tri[r,c] <- loading_mat[r,c]*mask_loading_mat[r,c]
      }
    }),

    ################
    ## Likelihood ##
    ################

    parsedLikelihood,

    #########################
    ## Derive the Quantity ##
    #########################
    # <1> Number of target people
    "N <- sum(z[1:M])",
    # <2> Overall baseline value of lambda (to compute baseline prevalence)
    "beta0 <- alpha", # compute overall baseline value of lambda
    # <3> Variance of spatial component
    # sigma2.phi <- 1/tau2 # variance of spatial component using inverse gamma for variance parameter
    "sigma2.phi <- inv_tau2", # variance of spatial component using half-cauchy for variance parameter
    # <4> Area-specific residual
    "resPhi[1:K] <- phi[1:K]",
    # <5> Overall residuals at individual level
    # res[1:M] <- (z[1:M]-lambda[1:M])/sqrt(lambda[1:M])
    # <6> Lambda
    expression(for(k in 1:K){
      pesudo_logit_prev[k] <- alpha + phi[k]
    })
  )

  return(model)
}


#' Full likelihood of data and prior and hyper
#'
#' @param tau_prior
#' hyper of phi
#' @param sigma_theta_prior
#' hyper of theta
#' @param theta_prior
#' prior theta
#'
#' @return
#' return the full likelihood of ICAR model
#' @export
#'
#' @examples
#' tau_prior <- c(inv_tau2  ~ T(dt(mu = 0, tau = sd_spRE_prior, df = 1), 0, ),"tau2 <- 1/inv_tau2")
#' sigma_theta_prior <- c(sigma_theta ~ T(dt(mu = 0, tau = sd_theta_prior, df = 1), 0, ))
#' theta_prior <- c(theta[i] ~ dnorm(mu_theta, sd = sigma_theta), "loading_latent_factor[i,1:J] <- theta[i] * t(loading_mat_upper_tri[1:J,1])")
#' model_ll <- hyper2Likelihood(tau_prior, sigma_theta_prior, theta_prior)
hyper2Likelihood <- function(tau_prior, sigma_theta_prior, theta_prior){
  # Return nimble code of ICAR Likelihood using hyper and prior
  likelihood_i      <- listLikelihoodIndividual(theta_prior)
  expr_likelihood_i <- list2expr(likelihood_i)
  parsedLikelihood  <- exprLikelihoodAllData(expr_likelihood_i)
  model             <- listFullLikelihoodwithHyper(tau_prior, sigma_theta_prior, parsedLikelihood)
  modelGraph        <- list2expr(model)
  return(modelGraph)
}

