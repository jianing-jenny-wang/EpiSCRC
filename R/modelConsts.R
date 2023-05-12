prepareICARmodelConst <- function(names_of_k, W_val){
  ###################################################################
  #### MODEL 3 - Assuming ICAR Model with Beta0, spRE, no Nugget ####
  ###################################################################
  ### Prepare INTRINSIC CAR (ICAR) model specification ###
  ## Reshape the adjacency matrix to get adjacency list
  ## (i.e. A vector of indices indicating which regions are neighbors of which.)
  ## Change the character name of Adj matrix to numeric name
  W_numeric_name <- data.frame(name = names_of_k, id = seq(1:length(names_of_k)))
  row.names(W_val) <- ifelse(row.names(W_val) %in% W_numeric_name$name, W_numeric_name$id,NA)
  colnames(W_val) <- ifelse(colnames(W_val) %in% W_numeric_name$name, W_numeric_name$id,NA)
  W <- W_val
  W_ls <- reshape2::melt(W)
  Adj_vec <- W_ls$Var1[W_ls$value==1]
  ## A vector of weights. In this case, all weights are 1.
  Weights_vec <- W_ls$value[W_ls$value==1]
  ## A vector of length N. num[n] indicates how many neighbors region n contains.
  ## This helps map the adj vector to the starting region.
  D_w <- diag(apply(W, 2, sum))
  D_w_ls <- reshape2::melt(D_w)
  D_w_ls <- D_w_ls$value[which(D_w_ls$value>0)]

  consts = list(Adj_vec = Adj_vec,
                Weights_vec = Weights_vec,
                D_w_ls = D_w_ls
  )
  return(consts)
}

getMask <- function(J, D){
  ## Create a mask matrix to obtain loading factor matrix with upper tri = 0 ##
  mask_loading_factor_mat <- matrix(1, nrow = J, ncol = D)
  mask_loading_factor_mat[upper.tri(mask_loading_factor_mat)] <- 0
  return(mask_loading_factor_mat)
}

getModelConsts <- function(aug_size, J, K, D, Adj_vec, Weights_vec, D_w_ls, SD_HalfCauchy_Prior, mask_loading_factor_mat, loc_label, sigma_mu = 10, sigma_loading = 10){
  # constants = Mod_GenMth_ICAR_Consts
  constants <- list(M = aug_size,
                    J = J,
                    K = K,
                    D = D, # length of latent factors in theta
                    S = diag(D), # Identity matrix of inverse wishart distribution for theta
                    L = length(Adj_vec),
                    adj = Adj_vec,
                    weights = Weights_vec,
                    num = D_w_ls, # number of neighbors
                    sigma_mu = sigma_mu, # list effects sigma
                    sd_spRE_prior = SD_HalfCauchy_Prior, # standard deviation of half-cauchy prior distribution for spRE variance parameter
                    mu_theta = rep(0, D), # mean of latent factors theta
                    sd_theta_prior = SD_HalfCauchy_Prior, # standard deviation of half-cauchy prior distribution for theta variance parameter
                    sigma_loading = sigma_loading,  # standard deviation of single loading element
                    mask_loading_mat = mask_loading_factor_mat, # TRUE indicates upper triangulation of the factor loading matrix
                    prev_re = loc_label # location label for each individual
  )
  return(constants)
}

