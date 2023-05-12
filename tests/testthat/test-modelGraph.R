test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("full likelihood equal nimbleCode",{
  nimble_code = prepareNimbleCode()
  tau_prior = c(inv_tau2  ~ T(dt(mu = 0, tau = sd_spRE_prior, df = 1), 0, ),
                "tau2 <- 1/inv_tau2"
  )
  sigma_theta_prior = c(sigma_theta ~ T(dt(mu = 0, tau = sd_theta_prior, df = 1), 0, ))

  theta_prior = c(theta[i] ~ dnorm(mu_theta, sd = sigma_theta),# vector of individual heterogeneity effect is a random (latent) effect
                  "loading_latent_factor[i,1:J] <- theta[i] * t(loading_mat_upper_tri[1:J,1])"# dim = 1 * trans(J*D)
  )
  model_ll <- hyper2Likelihood(tau_prior, sigma_theta_prior, theta_prior)

  out_char <- as.character(model_ll)
  exp_char <- as.character(nimble_code)
  testthat::expect_equal(length(out_char), length(exp_char))

  for(i in 1:length(out_char)){
    testthat::expect_equal(out_char[i], exp_char[i])
  }
})

