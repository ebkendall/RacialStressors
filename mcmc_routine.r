library(mvtnorm, quietly=T)
library(foreach, quietly=T)
library(doParallel, quietly=T)
library(deSolve, quietly=T)
library(LaplacesDemon, quietly=T)

# Construct the transition rate matrix
Q <- function(t,beta){

    betaMat = matrix(beta, ncol = 2, byrow = F) # determine the covariates
  
    q1  = exp( c(1,t) %*% betaMat[1,] )  # Transition from Base   to  Stress
    q2  = exp( c(1,t) %*% betaMat[2,] )  # Transition from Base   to  Recov
    q3  = exp( c(1,t) %*% betaMat[3,] )  # Transition from Stress to  Base
    q4  = exp( c(1,t) %*% betaMat[4,] )  # Transition from Stress to  Recov
    q5  = exp( c(1,t) %*% betaMat[5,] )  # Transition from Recov  to  Base
    q6  = exp( c(1,t) %*% betaMat[6,] )  # Transition from Recov  to  Stress
    
    qmat = matrix(c(  0,  q1,  q2,
                     q3,   0,  q4,
                     q5,  q6,   0),
                nrow = 3, byrow = T)
    diag(qmat) = -rowSums(qmat)

  return(qmat)
}

# Setting up the differential equations for deSolve
model_t <- function(t,p,parms) {
    qmat = Q(t, parms$b)
    pmat = matrix(c(  p[1],  p[2], p[3],
                      p[4],  p[5], p[6],
                      p[7],  p[8], p[9]),
                nrow = 3, byrow = T)
    
    # Vectorizing the matrix multiplication row-wise
    dP = c(t(pmat %*% qmat))
    return(list(dP))
}

# Gibbs update of the mu_tilde
update_mu_tilde = function(pars, par_index, mu_i, n_sub) {
    
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
update_upsilon = function(pars, par_index, mu_i, n_sub) {
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

# Evaluating the log posterior
fn_log_post_continuous <- function(pars, prior_par, par_index, y_1, y_2, t, id, mu_i) {

    # Order: Base, Stress, Recov
    init_logit = c( 1, exp(pars[par_index$pi_logit][1]), exp(pars[par_index$pi_logit][2]))

    # Initial state probabilities
    init = init_logit / sum(init_logit)

    # Misclassification response matrix
    resp_fnc = matrix(c(1, exp(pars[par_index$misclass][1]), exp(pars[par_index$misclass][2]),
                        exp(pars[par_index$misclass][3]), 1, exp(pars[par_index$misclass][4]),
                        exp(pars[par_index$misclass][5]), exp(pars[par_index$misclass][6]), 1),
                        ncol=3, byrow=TRUE)

    resp_fnc = resp_fnc / rowSums(resp_fnc)
    
    tau2 <- exp(pars[par_index$log_tau2])
    
    beta <- pars[par_index$beta]

    # initial condition for deSolve
    p_ic <- c( p1=1, p2=0, p3=0,
               p4=0, p5=1, p6=0,
               p7=0, p8=0, p9=1)
    
    eids = unique(id)
  
    # Parallelized computation of the log-likelihood
    log_total_val = foreach(i=unique(id), .combine='+', 
                            .export = c("model_t", "Q"), 
                            .packages = c("deSolve"),
                            .inorder=FALSE) %dopar% {
        
        f_i = val = 1
        y_1_i = y_1[id == i]    # the observed state
        # y_2_i = y_2[id == i]    # the rsa measurements
        t_i = t[id == i]        # time
        
        mu_1 = 1 # dnorm(x = y_2_i[1], mean = mu_i[which(eids == i), 1], sd = sqrt(tau2))
        mu_2 = 1 # dnorm(x = y_2_i[1], mean = mu_i[which(eids == i), 2], sd = sqrt(tau2))
        mu_3 = 1 # dnorm(x = y_2_i[1], mean = mu_i[which(eids == i), 3], sd = sqrt(tau2))
        
        if(y_1_i[1] <= 3) { # observed state
          f_i = init %*% diag(c(mu_1,mu_2,mu_3) * resp_fnc[, y_1_i[1]])
        } else { # un-observed state (99)
            print("un-observed state")
            f_i = init %*% diag(c(mu_1,mu_2,mu_3))
        }

        log_norm = 0
        
        for(k in 2:length(t_i)) {
            out <- deSolve::ode(p_ic, times = t_i[(k-1):k], 
                                      func = model_t, 
                                      parms = list(b=beta))
            
            P <- matrix(c( out[2,"p1"],  out[2,"p2"],  out[2,"p3"],  
                           out[2,"p4"],  out[2,"p5"],  out[2,"p6"],  
                           out[2,"p7"],  out[2,"p8"],  out[2,"p9"]),
                        nrow = 3, byrow = T)
            
            mu_1 = 1 # dnorm(x = y_2_i[k], mean = mu_i[which(eids == i), 1], sd = sqrt(tau2))
            mu_2 = 1 # dnorm(x = y_2_i[k], mean = mu_i[which(eids == i), 2], sd = sqrt(tau2))
            mu_3 = 1 # dnorm(x = y_2_i[k], mean = mu_i[which(eids == i), 3], sd = sqrt(tau2))

            if(y_1_i[k] <= 3) { # observed state
              D_i = diag(c(mu_1,mu_2,mu_3) * resp_fnc[, y_1_i[k]])
            } else { # unknown state (99)
              D_i = diag(c(mu_1,mu_2,mu_3))
            }

            val = f_i %*% P %*% D_i

            norm_val = sqrt(sum(val^2))
            f_i = val / norm_val
            log_norm = log_norm + log(norm_val)
        }

        return(log(sum(f_i)) + log_norm)
    }

    mean = prior_par$prior_mean
    sd = diag(prior_par$prior_sd)
    log_prior_dens = dmvnorm( x=pars[c(c(par_index$beta), c(par_index$misclass), 
                                       c(par_index$pi_logit))], #, c(par_index$log_tau2)
                                       mean=mean, sigma=sd, log=T)
    return(log_total_val + log_prior_dens)

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
  chain = matrix( 0, steps, n_par)

  # group = as.list(unlist(par_index))
  # names(group) = NULL
  # group = list(c(par_index$beta[c(1,7)]), c(par_index$beta[c(2,8)]), c(par_index$beta[c(3,9)]),
  #              c(par_index$beta[c(4,10)]), c(par_index$beta[c(5,11)]), c(par_index$beta[c(6,12)]),
  #              c(par_index$misclass), c(par_index$pi_logit), c(par_index$log_tau2))
  group = list(c(par_index$beta), c(par_index$misclass), c(par_index$pi_logit)) # , c(par_index$log_tau2)
  n_group = length(group)

  # proposal covariance and scale parameter for Metropolis step
  # pcov = list();	for(j in 1:n_group)  pcov[[j]] = diag(length(group[[j]]))*0.001
  # pscale = rep( 1, n_group)
  load(paste0('Model_out/mcmc_out_4_5.rda'))
  pcov = mcmc_out$pcov
  pscale = mcmc_out$pscale
  rm(mcmc_out)

  accept = rep( 0, n_group)
  
  # Initializing the random effects matrix mu_i
  # M = vector(mode = "list", length = 10)
  # mu_i = matrix(rep(c(8, 6, 7), n_sub), nrow = n_sub, ncol = 3, byrow = T)

  # Evaluate the log_post of the initial parameters
  log_post_prev = fn_log_post_continuous( pars, prior_par, par_index, y_1, y_2, t, id, mu_i)

  if(!is.finite(log_post_prev)){
    print("Infinite log-posterior; choose better initial parameters")
    break
  }

  # Begin the MCMC algorithm --------------------------------------------------
  chain[1,] = pars
  for(ttt in 2:steps){
      
    # # mu_tilde: Gibbs update
    # pars[par_index$mu_tilde] = update_mu_tilde(pars, par_index, mu_i, n_sub)
    chain[ttt, par_index$mu_tilde] = pars[par_index$mu_tilde]
    
    # # upsilon: Gibbs update
    # pars[par_index$upsilon] = update_upsilon(pars, par_index, mu_i, n_sub)
    chain[ttt, par_index$upsilon] = pars[par_index$upsilon]
    
    # # mu_i: update
    # var_mu_i = matrix(pars[par_index$upsilon], nrow = 3, ncol = 3)
    # mu_i = rmvnorm(n_sub, mean = pars[par_index$mu_tilde], sigma = var_mu_i)
      
    for(j in 1:n_group){

      # Propose an update
      ind_j = group[[j]]
      proposal = pars
      if(length(ind_j > 1)) {
          proposal[ind_j] = rmvnorm( n=1, mean=pars[ind_j],sigma=pcov[[j]]*pscale[j])
      } else {
          proposal[ind_j] = rnorm( n=1, mean=pars[ind_j],sd=sqrt(pcov[[j]]*pscale[j]))
      }

      # Compute the log density for the proposal
      log_post = fn_log_post_continuous(proposal, prior_par, par_index, y_1, y_2, t, id, mu_i)

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
          
          log_post = fn_log_post_continuous(proposal, prior_par, par_index, y_1, y_2, t, id, mu_i)
        }
      }

      # Evaluate the Metropolis-Hastings ratio
      if( log_post - log_post_prev > log(runif(1,0,1)) ){
        log_post_prev = log_post
        pars[ind_j] = proposal[ind_j]
        accept[j] = accept[j] +1
      }
      chain[ttt,ind_j] = pars[ind_j]

      # Proposal tuning scheme ------------------------------------------------
      if(ttt < burnin){
        # During the burnin period, update the proposal covariance in each step
        # to capture the relationships within the parameters vectors for each
        # transition.  This helps with mixing.
        if(ttt == 100)  pscale[j] = 1
        
        if(ttt %% 100 == 0) {
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

    if(ttt%%1==0)  cat('--->',ttt,'\n')
  }
  # ---------------------------------------------------------------------------

  stopCluster(cl)
  print(accept/(steps-burnin))

  return(list( chain=chain[burnin:steps,], accept=accept/(steps-burnin),
               pscale=pscale, pcov = pcov))
}
# -----------------------------------------------------------------------------
