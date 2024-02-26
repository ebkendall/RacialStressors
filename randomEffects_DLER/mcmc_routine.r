library(mvtnorm, quietly=T)
library(LaplacesDemon, quietly=T)

# Rcpp packages
library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)
sourceCpp("mcmc_routine_c.cpp")

# -----------------------------------------------------------------------------
# The mcmc routine for samping the parameters
# -----------------------------------------------------------------------------
mcmc_routine = function( y_1, y_2, t, id, init_par, prior_par, par_index, steps,
                         burnin, n_sub, case_b, cov_info, simulation, B, 
                         covariate_struct, big_steps, trial_num, ind){

    pars = init_par
    n = length(y_1)
    n_par = length(pars)
    chain = matrix( 0, steps, n_par)
    B_chain = matrix( 0, steps, length(y_1))

    if(case_b) {
        if(covariate_struct == 1) {
            group = list(c(par_index$zeta[1:6]), c(par_index$zeta[7:12]),
                         c(par_index$zeta[13:18]), c(par_index$zeta[19:24]),
                         c(par_index$delta), c(par_index$tau2,par_index$sigma2),
                         c(par_index$gamma))
        } else if(covariate_struct == 2) {
            group = list(c(par_index$zeta[1:6]), c(par_index$zeta[7:12]),
                         c(par_index$delta),c(par_index$tau2, par_index$sigma2),
                         c(par_index$delta_new))
        } else {
            group = list(c(par_index$zeta[1:6]), c(par_index$zeta[7:12]),
                         c(par_index$zeta[13:18]), c(par_index$zeta[19:24]),
                         c(par_index$zeta[25:30]),
                         c(par_index$delta), c(par_index$tau2,par_index$sigma2),
                         c(par_index$gamma), c(par_index$delta_new))
        }
    } else {
        if(covariate_struct == 1) {
            group = list(c(par_index$zeta[1:6]), c(par_index$zeta[7:12]),
                         c(par_index$zeta[13:18]), c(par_index$zeta[19:24]),
                         c(par_index$delta), c(par_index$tau2,par_index$sigma2),
                         c(par_index$gamma), c(par_index$misclass))
        } else if(covariate_struct == 2) {
            group = list(c(par_index$zeta[1:6]), c(par_index$zeta[7:12]),
                         c(par_index$delta),c(par_index$tau2, par_index$sigma2),
                         c(par_index$delta_new), c(par_index$misclass))
        } else {
            group = list(c(par_index$zeta[1:6]), c(par_index$zeta[7:12]),
                         c(par_index$zeta[13:18]), c(par_index$zeta[19:24]),
                         c(par_index$zeta[25:30]), c(par_index$delta), 
                         c(par_index$tau2,par_index$sigma2),
                         c(par_index$gamma), c(par_index$delta_new), 
                         c(par_index$misclass))
        }
    }

    names(group) = NULL
    n_group = length(group)

    # proposal covariance and scale parameter for Metropolis step
    # load(paste0("Model_out/mcmc_out_", ind, "_", trial_num - 3, "_30b.rda"))
    # pcov = mcmc_out$pcov
    # pscale = mcmc_out$pscale
    # rm(mcmc_out)
    pcov = list(); for(j in 1:n_group)  pcov[[j]] = diag(length(group[[j]]))*0.001
    pscale = rep( 1, n_group)
    
    accept = rep( 0, n_group)
    EIDs = unique(id)
    
    # Begin the MCMC algorithm -------------------------------------------------
    chain[1,] = pars
    for(ttt in 2:big_steps){

        chain_ind = ttt %% steps
        if(chain_ind == 0) chain_ind = steps
        
        chain[chain_ind,] = pars
        if(ttt %% 200 == 0) {
            print("pars")
            print(chain[chain_ind - 1,])
        }
        
        # Update the state space
        new_B = update_b_i(EIDs, pars, par_index, id, B, y_2, y_1, cov_info, 
                           case_b, covariate_struct)
        B = new_B
        
        # Evaluate the log_post of the initial parameters
        log_post_prev = fn_log_post_continuous(EIDs, pars, prior_par, par_index, 
                                               y_1, id, y_2, cov_info, case_b, 
                                               B, covariate_struct)
        
        if(!is.finite(log_post_prev)){
            print("Infinite log-posterior; choose better initial parameters")
            break
        }

        for(j in 1:n_group){

            # Propose an update
            ind_j = group[[j]]
            proposal = pars
            proposal[ind_j] = rmvnorm( n=1, mean=pars[ind_j],
                                        sigma=pcov[[j]]*pscale[j])
            # Guarantee alpha < beta
            if(sum(ind_j %in% par_index$delta) == 3) {
                while(proposal[par_index$delta][2] >= proposal[par_index$delta][3]) {
                    proposal[ind_j] = rmvnorm( n=1, mean=pars[ind_j],
                                                sigma=pcov[[j]]*pscale[j])
                }
            }

            # Compute the log density for the proposal
            log_post = fn_log_post_continuous(EIDs, proposal, prior_par,
                                              par_index, y_1, id, y_2, cov_info,
                                              case_b, B, covariate_struct)

            # Only propose valid parameters during the burnin period
            if(ttt < burnin){
                while(!is.finite(log_post)){
                    print(paste0("bad proposal, ", j, ":"))
                    print(proposal[ind_j])

                    proposal = pars
                    
                    proposal[ind_j] = rmvnorm( n=1, mean=pars[ind_j],
                                                sigma=pcov[[j]]*pscale[j])
                    # Guarantee alpha < beta
                    if(sum(ind_j %in% par_index$delta) == 3) {
                        while(proposal[par_index$delta][2] >= proposal[par_index$delta][3]) {
                            proposal[ind_j] = rmvnorm( n=1, mean=pars[ind_j],
                                                        sigma=pcov[[j]]*pscale[j])
                        }
                    }
                    
                    log_post = fn_log_post_continuous(EIDs, proposal, prior_par,
                                                    par_index, y_1, id, y_2, cov_info, 
                                                    case_b, B, covariate_struct)
                }
            }
            
            # Ensuring that we do not have problems from C++
            if(!is.finite(log_post) | is.nan(log_post)) {
                print(paste0("bad proposal post burnin: ", log_post))
            }
            
            if( log_post - log_post_prev > log(runif(1,0,1)) ){
                log_post_prev = log_post
                pars[ind_j] = proposal[ind_j]
                accept[j] = accept[j] +1
            }
            
            chain[chain_ind,ind_j] = pars[ind_j]

            # Proposal tuning scheme ------------------------------------------------
            if(ttt < burnin){
                # During the burnin period, update the proposal covariance in each step
                # to capture the relationships within the parameters vectors for each
                # transition.  This helps with mixing.
                if(ttt == 100)  pscale[j] = 1
                
                if(100 <= ttt & ttt <= 2000){
                    temp_chain = chain[1:ttt,ind_j]
                    pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])
                    
                } else if(2000 < ttt){
                    temp_chain = chain[(ttt-2000):ttt,ind_j]
                    pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])
                }
                
                if( sum( is.na(pcov[[j]]) ) > 0)  pcov[[j]] = diag( length(ind_j) )

                # Tune the proposal covariance for each transition to achieve
                # reasonable acceptance ratios.
                if(ttt %% 30 == 0){
                    if(ttt %% 480 == 0){
                        accept[j] = 0
                    } else if( accept[j] / (ttt %% 480) < .4 ){ 
                        pscale[j] = (.5^2)*pscale[j]

                    } else if( accept[j] / (ttt %% 480) > .5 ){ 
                        pscale[j] = (1.5^2)*pscale[j]
                    }
                }
            }
            # -----------------------------------------------------------------------
        }

        # Restart the acceptance ratio at burnin.
        if(ttt == burnin){ accept = rep( 0, n_group) }
        if(ttt > burnin){ B_chain[ chain_ind, ] = do.call( 'c', B) }
        
        if(ttt%%1==0)  cat('--->',ttt,'\n')

        if(ttt > burnin & ttt%%steps == 0) {
            if(ttt/steps == 1) {
                mcmc_out = list(chain=chain[burnin:steps,], B_chain=B_chain,
                                accept=accept/(steps-burnin), pscale=pscale, 
                                pcov = pcov)
                ind_keep = seq(burnin + 10, steps, by=10)
                chain_big   = chain[ind_keep,]
                B_chain_big = B_chain[ind_keep, ]
            } else {
                mcmc_out = list(chain=chain, B_chain=B_chain,
                                accept=accept/(steps-burnin), pscale=pscale, 
                                pcov = pcov)
                ind_keep = seq(10, steps, by=10)
                chain_big   = rbind(chain_big, chain[ind_keep,])
                B_chain_big = rbind(B_chain_big, B_chain[ind_keep, ])
            }

            # Save intermediary chain
            if(ttt/steps == 5) {
                if(simulation) {
                    save(mcmc_out, file = paste0('Model_out/mcmc_out_interm_',ind,'_', 
                                                trial_num, 'it', ttt/steps, '_sim_check.rda'))
                } else {
                    save(mcmc_out, file = paste0('Model_out/mcmc_out_interm_',ind,'_', 
                                                trial_num, 'it', ttt/steps, '_check.rda'))
                }
            }

            # Reset the chains
            chain = matrix( 0, steps, n_par)
            B_chain = matrix( 0, steps, length(y_1))
        }
    }
    # ---------------------------------------------------------------------------

    return(list(chain=chain_big, B_chain=B_chain_big,
                pscale=pscale, pcov = pcov))
}
# -----------------------------------------------------------------------------
