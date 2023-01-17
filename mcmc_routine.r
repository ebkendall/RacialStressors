library(mvtnorm, quietly=T)
library(foreach, quietly=T)
library(doParallel, quietly=T)
library(deSolve, quietly=T)
library(LaplacesDemon, quietly=T)

# Rcpp packages
library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)
sourceCpp("mcmc_routine_c.cpp")

# Needed for OpenMP C++ parallel
Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
Sys.setenv("PKG_LIBS" = "-fopenmp")

# Gibbs update of the mu_tilde
update_mu_tilde = function(pars, par_index, n_sub) {
    
    mu_i = matrix(pars[par_index$mu_i], ncol = 3)
    
    # Prior mean and variance for mu_tilde
    big_sigma = matrix(c(1.097200, 1.014710, 1.134521,
                         1.014710, 1.332313, 1.475067,
                         1.134521, 1.475067, 2.024380), ncol = 3, byrow = T)
    big_sigma_inv = solve(big_sigma)
    mu_0 = c(6.411967, 6.481880, 6.335972)
    
    # variance for mu^(i)
    upsilon = matrix(pars[par_index$upsilon], nrow = 3, ncol = 3)
    upsilon_inv = solve(upsilon)
    
    # variance of Gibbs update
    V = solve(big_sigma_inv + n_sub * upsilon_inv)
    
    # Mean of Gibbs update
    M = V %*% (big_sigma_inv %*% mu_0 + n_sub * (upsilon_inv %*% colMeans(mu_i)) )
    
    return(rmvnorm(1, mean = M, sigma = V))
}

# Gibbs update of the mu_tilde
update_upsilon = function(pars, par_index, n_sub) {
    
    mu_i = matrix(pars[par_index$mu_i], ncol = 3)
    
    sum_mu_i = (mu_i[1, ] - pars[par_index$mu_tilde]) %*% t(mu_i[1, ] - pars[par_index$mu_tilde])
    for(i in 2:n_sub) {
        sum_mu_i = sum_mu_i + (mu_i[i, ] - pars[par_index$mu_tilde]) %*% t(mu_i[i, ] - pars[par_index$mu_tilde])
    }
    
    # Prior for Upsilon
    psi = diag(3)
    nu = 3 + 2
    
    # Gibbs update
    new_psi = psi + sum_mu_i
    new_nu  = nu + n_sub
    
    return(c(rinvwishart(new_nu, new_psi)))
}

# Gibbs update of the mu_i
update_mu_i = function(y_2, pars, par_index, n_sub, V_i, EIDs, id) {
    mu_i = foreach(i=1:n_sub, .combine = 'rbind', .packages = "mvtnorm") %dopar% {
        
        D_small = V_i[[i]]
        upsilon = matrix(pars[par_index$upsilon], nrow = 3, ncol = 3)
        up_solve = solve(upsilon)
        
        y_sub = matrix(y_2[id == EIDs[i]], ncol=1)
        tau2 = pars[par_index$tau2]
        
        V = solve((1/tau2) * (t(D_small) %*% D_small) + up_solve)
        M = V %*% ((1/tau2) * (t(D_small) %*% y_sub) + up_solve %*% matrix(pars[par_index$mu_tilde],ncol=1))
        
        mu_i_small = c(rmvnorm(n = 1, mean = M, sigma = V))
        
        return(mu_i_small)
    }
    
    return(mu_i)
}

# Gibbs update of the tau2
update_tau2 = function(y_2, pars, par_index, V_i, n_sub, EIDs, id) {
    
    a = b = 1
    mu_i = matrix(pars[par_index$mu_i], ncol = 3)
    
    temp = 0
    for(i in 1:n_sub) {
        y_sub = matrix(y_2[id == EIDs[i]], ncol=1)
        mu_temp = t(mu_i[i,,drop=F])
        temp = temp + t(y_sub - V_i[[i]] %*% mu_temp) %*% (y_sub - V_i[[i]] %*% mu_temp)
    }
    
    a_new = a + 0.5 * length(y_2)
    b_new = b + 0.5 * temp
    
    return(rinvgamma(n=1, shape = a_new, scale = b_new))
}

update_V_i = function(B) {
  V_i = vector(mode = 'list', length = length(B))
  for(i in 1:length(B)) {
    state_sub = c(B[[i]])
    V_i[[i]] = matrix(nrow = length(state_sub), ncol = 3)
    V_i[[i]][,1] = as.integer(state_sub == 1)
    V_i[[i]][,2] = as.integer(state_sub == 2)
    V_i[[i]][,3] = as.integer(state_sub == 3)
  }
  return(V_i)
}

# -----------------------------------------------------------------------------
# The mcmc routine for samping the parameters
# -----------------------------------------------------------------------------
mcmc_routine = function( y_1, y_2, t, id, init_par, prior_par, par_index,
                         steps, burnin, n_cores, n_sub){

  cl <- makeCluster(n_cores, outfile="")
  registerDoParallel(cl)

  pars = init_par
  n = length(y_1)
  n_par = length(pars)
  chain = matrix( 0, steps, n_par - length(par_index$mu_i))
  B_chain = matrix( 0, steps - burnin, length(y_1))

  group = list(c(par_index$beta), c(par_index$misclass))
  n_group = length(group)

  # proposal covariance and scale parameter for Metropolis step
  pcov = list();	for(j in 1:n_group)  pcov[[j]] = diag(length(group[[j]]))*0.001
  pscale = rep( 1, n_group)
  # load(paste0('Model_out/mcmc_out_4_13.rda'))
  # pcov = mcmc_out$pcov
  # pscale = mcmc_out$pscale
  # rm(mcmc_out)

  accept = rep( 0, n_group)

  EIDs = unique(id)
  
  # Initializing the state space list B
  B = list()
  for(i in 1:length(EIDs)) {
    state_sub = y_1[id == EIDs[i]]
    b_temp = matrix(state_sub, ncol = 1)
    B[[i]] = b_temp
  }

  # Initializing the V_i matrix
  V_i = update_V_i(B)

  # Vector to store the values of mu_i
  big_mu_i = vector(mode = 'list', length = 10)

  
  # Evaluate the log_post of the initial parameters
  log_post_prev = log_f_i_cpp_total(EIDs, pars, prior_par, par_index, y_1, t, id, B)

  if(!is.finite(log_post_prev)){
    print("Infinite log-posterior; choose better initial parameters")
    break
  }
 
  # Begin the MCMC algorithm --------------------------------------------------
  chain[1,] = pars[-par_index$mu_i]
  for(ttt in 2:steps){
      
    # mu_tilde: Gibbs update
    pars[par_index$mu_tilde] = update_mu_tilde(pars, par_index, n_sub)
    chain[ttt, par_index$mu_tilde] = pars[par_index$mu_tilde]
    
    # upsilon: Gibbs update
    pars[par_index$upsilon] = update_upsilon(pars, par_index, n_sub)
    chain[ttt, par_index$upsilon] = pars[par_index$upsilon]
    
    # mu_i: Gibbs update
    pars[par_index$mu_i] = update_mu_i_cpp(y_2, pars, par_index, V_i, EIDs, id)
    if(ttt %% 1000 == 0) M[[ttt/1000]] = pars[par_index$mu_i]

    # tau2: Gibbs update
    pars[par_index$tau2] = update_tau2(y_2, pars, par_index, V_i, n_sub, EIDs, id)
    chain[ttt, par_index$tau2] = pars[par_index$tau2]
    
    # S_chain: Metropolis-within-Gibbs update
    B_V = update_b_i_cpp(8, EIDs, pars, prior_par, par_index, y_1, t, id, B, V_i)
    B = B_V[[1]]
    V_i = B_V[[2]]
      
    for(j in 1:n_group){

      # Propose an update
      ind_j = group[[j]]
      proposal = pars
      if(length(ind_j) > 1) {
          proposal[ind_j] = rmvnorm( n=1, mean=pars[ind_j],sigma=pcov[[j]]*pscale[j])
      } else {
          proposal[ind_j] = rnorm( n=1, mean=pars[ind_j],sd=sqrt(pcov[[j]]*pscale[j]))
      }

      # Compute the log density for the proposal
      log_post = log_f_i_cpp_total(EIDs, proposal, prior_par, par_index, y_1, t, id, B)

      # Only propose valid parameters during the burnin period
      if(ttt < burnin){
        while(!is.finite(log_post)){
          print('bad proposal')
          proposal = pars
          if(length(ind_j > 1)) {
              proposal[ind_j] = rmvnorm( n=1, mean=pars[ind_j],sigma=pcov[[j]]*pscale[j])
          } else {
              proposal[ind_j] = rnorm( n=1, mean=pars[ind_j],sd=sqrt(pcov[[j]]*pscale[j]))
          }
          
          log_post = log_f_i_cpp_total(EIDs, proposal, prior_par, par_index, y_1, t, id, B)
        }
      }

      # Evaluate the Metropolis-Hastings ratio
      if( log_post - log_post_prev > log(runif(1,0,1)) ){
        log_post_prev = log_post
        pars[ind_j] = proposal[ind_j]
        accept[j] = accept[j] +1
        # print("accept"); print(accept)
      }
      chain[ttt,ind_j] = pars[ind_j]

      # Proposal tuning scheme ------------------------------------------------
      if(ttt < burnin){
        # During the burnin period, update the proposal covariance in each step
        # to capture the relationships within the parameters vectors for each
        # transition.  This helps with mixing.
        if(ttt == 100)  pscale[j] = 1
        
        if(ttt %% 50 == 0) {
            print("accept ratio")
            print(accept[j])
        }
        
        if (length(ind_j) > 1) {
            if(100 <= ttt & ttt <= 2000){
              temp_chain = chain[1:ttt,ind_j]
              pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])
    
            } else if(2000 < ttt){
              temp_chain = chain[(ttt-2000):ttt,ind_j]
              pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])
            }
        } else {
            if(100 <= ttt & ttt <= 2000){
                temp_chain = chain[1:ttt,ind_j]
                pcov[[j]] = matrix(var(temp_chain[ !duplicated(temp_chain)]))
                
            } else if(2000 < ttt){
                temp_chain = chain[(ttt-2000):ttt,ind_j]
                pcov[[j]] = matrix(var(temp_chain[ !duplicated(temp_chain)]))
            }
        }
        
        if( sum( is.na(pcov[[j]]) ) > 0)  pcov[[j]] = diag( length(ind_j) )

        # Tune the proposal covariance for each transition to achieve
        # reasonable acceptance ratios.
        if(ttt %% 30 == 0){
          if(ttt %% 480 == 0){
            accept[j] = 0

          } else if( accept[j] / (ttt %% 480) < .4 ){ 
            pscale[j] = (.75^2)*pscale[j]

          } else if( accept[j] / (ttt %% 480) > .5 ){ 
            pscale[j] = (1.25^2)*pscale[j]
          }
        }
      }
      # -----------------------------------------------------------------------
    }

    # Restart the acceptance ratio at burnin.
    if(ttt == burnin)  accept = rep( 0, n_group)
    if(ttt > burnin) {
      B_chain[ttt - burnin, ] = do.call( 'c', B)
    }

    if(ttt%%1==0)  cat('--->',ttt,'\n')
  }
  # ---------------------------------------------------------------------------

  stopCluster(cl)
  print(accept/(steps-burnin))

  return(list( chain=chain[burnin:steps,], B_chain = B_chain,
               accept=accept/(steps-burnin),
               pscale=pscale, pcov = pcov, M = M))
}
# -----------------------------------------------------------------------------
