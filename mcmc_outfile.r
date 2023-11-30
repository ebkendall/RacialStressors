# This script file produces trace plots and histograms of the mcmc output files
library(latex2exp)

dir = 'Model_out/' 

# Information defining which approach to take ----------------------------------
trial_num = 1
simulation = F
thirty = T
case_b = F
# ------------------------------------------------------------------------------

# Size of posterior sample from mcmc chains
n_post = 45000
# Step number at which the adaptive tuning scheme was frozen
burnin = 5000
# Total number of steps the mcmc algorithm is computed for
steps = 50000
# Matrix row indices for the posterior sample to use
index_post = (steps - burnin - n_post + 1):(steps - burnin)

index_seeds = c(1:5)

if(trial_num >= 3) {
    par_index = list(zeta=1:30, misclass=0,delta = 31:33, tau2 = 34, sigma2 = 35:37,
                     gamma = 38:41)
    labels <- c(TeX(r'($\hat{\zeta}_{0,1}:$ Baseline: 1 $\to$ 2)'), 
                TeX(r'($\hat{\zeta}_{0,2}:$ Baseline: 1 $\to$ 3)'), 
                TeX(r'($\hat{\zeta}_{0,3}:$ Baseline: 2 $\to$ 1)'),
                TeX(r'($\hat{\zeta}_{0,4}:$ Baseline: 2 $\to$ 3)'),
                TeX(r'($\hat{\zeta}_{0,5}:$ Baseline: 3 $\to$ 1)'),
                TeX(r'($\hat{\zeta}_{0,6}:$ Baseline: 3 $\to$ 2)'),
                TeX(r'($\hat{\zeta}_{1,1}:$ age: 1 $\to$ 2)'), 
                TeX(r'($\hat{\zeta}_{1,2}:$ age: 1 $\to$ 3)'), 
                TeX(r'($\hat{\zeta}_{1,3}:$ age: 2 $\to$ 1)'),
                TeX(r'($\hat{\zeta}_{1,4}:$ age: 2 $\to$ 3)'),
                TeX(r'($\hat{\zeta}_{1,5}:$ age: 3 $\to$ 1)'),
                TeX(r'($\hat{\zeta}_{1,6}:$ age: 3 $\to$ 2)'),
                TeX(r'($\hat{\zeta}_{2,1}:$ sex1: 1 $\to$ 2)'), 
                TeX(r'($\hat{\zeta}_{2,2}:$ sex1: 1 $\to$ 3)'), 
                TeX(r'($\hat{\zeta}_{2,3}:$ sex1: 2 $\to$ 1)'),
                TeX(r'($\hat{\zeta}_{2,4}:$ sex1: 2 $\to$ 3)'),
                TeX(r'($\hat{\zeta}_{2,5}:$ sex1: 3 $\to$ 1)'),
                TeX(r'($\hat{\zeta}_{2,6}:$ sex1: 3 $\to$ 2)'),
                TeX(r'($\hat{\zeta}_{3,1}:$ yes edu: 1 $\to$ 2)'), 
                TeX(r'($\hat{\zeta}_{3,2}:$ yes edu: 1 $\to$ 3)'), 
                TeX(r'($\hat{\zeta}_{3,3}:$ yes edu: 2 $\to$ 1)'),
                TeX(r'($\hat{\zeta}_{3,4}:$ yes edu: 2 $\to$ 3)'),
                TeX(r'($\hat{\zeta}_{3,5}:$ yes edu: 3 $\to$ 1)'),
                TeX(r'($\hat{\zeta}_{3,6}:$ yes edu: 3 $\to$ 2)'),
                TeX(r'($\hat{\zeta}_{4,1}:$ DLER: 1 $\to$ 2)'), 
                TeX(r'($\hat{\zeta}_{4,2}:$ DLER: 1 $\to$ 3)'), 
                TeX(r'($\hat{\zeta}_{4,3}:$ DLER: 2 $\to$ 1)'),
                TeX(r'($\hat{\zeta}_{4,4}:$ DLER: 2 $\to$ 3)'),
                TeX(r'($\hat{\zeta}_{4,5}:$ DLER: 3 $\to$ 1)'),
                TeX(r'($\hat{\zeta}_{4,6}:$ DLER: 3 $\to$ 2)'),
                TeX(r'($\delta_1 = \mu$)'), TeX(r'($\delta_2 = \alpha$)'), TeX(r'($\delta_3 = \beta$)'),
                TeX(r'($\log(\tau^2)$)'), TeX(r'($\log(\sigma_1^2)$)'), TeX(r'($\log(\sigma_2^2)$)'), TeX(r'($\log(\sigma_3^2)$)'),
                TeX(r'($\hat{\gamma}_1:$ age)'), TeX(r'($\hat{\gamma}_2:$ sex1)'), 
                TeX(r'($\hat{\gamma}_3:$ yes edu)'), TeX(r'($\hat{\gamma}_4:$ DLER)'),
                TeX(r'($\tau^2 + \sigma_1^2$)'))
} else {
    par_index = list( zeta=1:25, misclass=0, delta = 26:28, tau2 = 29, 
                    sigma2 = 30:32, beta = 33:36)
    labels <- c(TeX(r'($\hat{\zeta}_{0,1}:$ Baseline: 1 $\to$ 2)'), 
                TeX(r'($\hat{\zeta}_{0,2}:$ Baseline: 2 $\to$ 1)'), 
                TeX(r'($\hat{\zeta}_{0,3}:$ Baseline: 2 $\to$ 3)'),
                TeX(r'($\hat{\zeta}_{0,4}:$ Baseline: 3 $\to$ 1)'),
                TeX(r'($\hat{\zeta}_{0,5}:$ Baseline: 3 $\to$ 2)'),
                TeX(r'($\hat{\zeta}_{1,1}:$ age: 1 $\to$ 2)'), 
                TeX(r'($\hat{\zeta}_{1,2}:$ age: 2 $\to$ 1)'), 
                TeX(r'($\hat{\zeta}_{1,3}:$ age: 2 $\to$ 3)'),
                TeX(r'($\hat{\zeta}_{1,4}:$ age: 3 $\to$ 1)'),
                TeX(r'($\hat{\zeta}_{1,5}:$ age: 3 $\to$ 2)'),
                TeX(r'($\hat{\zeta}_{2,1}:$ sex1: 1 $\to$ 2)'), 
                TeX(r'($\hat{\zeta}_{2,2}:$ sex1: 2 $\to$ 1)'), 
                TeX(r'($\hat{\zeta}_{2,3}:$ sex1: 2 $\to$ 3)'),
                TeX(r'($\hat{\zeta}_{2,4}:$ sex1: 3 $\to$ 1)'),
                TeX(r'($\hat{\zeta}_{2,5}:$ sex1: 3 $\to$ 2)'),
                TeX(r'($\hat{\zeta}_{3,1}:$ yes edu: 1 $\to$ 2)'), 
                TeX(r'($\hat{\zeta}_{3,2}:$ yes edu: 2 $\to$ 1)'), 
                TeX(r'($\hat{\zeta}_{3,3}:$ yes edu: 2 $\to$ 3)'),
                TeX(r'($\hat{\zeta}_{3,4}:$ yes edu: 3 $\to$ 1)'),
                TeX(r'($\hat{\zeta}_{3,5}:$ yes edu: 3 $\to$ 2)'),
                TeX(r'($\hat{\zeta}_{4,1}:$ DLER: 1 $\to$ 2)'), 
                TeX(r'($\hat{\zeta}_{4,2}:$ DLER: 2 $\to$ 1)'), 
                TeX(r'($\hat{\zeta}_{4,3}:$ DLER: 2 $\to$ 3)'),
                TeX(r'($\hat{\zeta}_{4,4}:$ DLER: 3 $\to$ 1)'),
                TeX(r'($\hat{\zeta}_{4,5}:$ DLER: 3 $\to$ 2)'),
                TeX(r'($\delta_1 = \mu$)'), TeX(r'($\delta_2 = \alpha$)'), TeX(r'($\delta_3 = \beta$)'),
                TeX(r'($\log(\tau^2)$)'), TeX(r'($\log(\sigma_1^2)$)'), TeX(r'($\log(\sigma_2^2)$)'), TeX(r'($\log(\sigma_3^2)$)'),
                TeX(r'($\hat{\gamma}_1:$ age)'), TeX(r'($\hat{\gamma}_2:$ sex1)'), 
                TeX(r'($\hat{\gamma}_3:$ yes edu)'), TeX(r'($\hat{\gamma}_4:$ DLER)'),
                TeX(r'($\tau^2 + \sigma_1^2$)'))
}

# -----------------------------------------------------------------------------
# Create mcmc trace plots and histograms
# -----------------------------------------------------------------------------
chain_list = vector(mode = "list", length = length(index_seeds))
post_means = matrix(nrow = length(index_seeds), ncol = length(labels))

ind = 0

for(seed in index_seeds){
    
    file_name = NULL
    
    if(simulation) {
        if(thirty) {
            if(case_b) {
                file_name = paste0(dir,'mcmc_out_',toString(seed), '_', trial_num, '_sim_30b.rda')   
            } else {
                file_name = paste0(dir,'mcmc_out_',toString(seed), '_', trial_num, '_sim_30.rda')      
            }
        } else {
            if(case_b) {
                file_name = paste0(dir,'mcmc_out_',toString(seed), '_', trial_num, '_sim_15b.rda')
            } else {
                file_name = paste0(dir,'mcmc_out_',toString(seed), '_', trial_num, '_sim_15.rda')   
            }
        }
    } else {
        if(thirty) {
            if(case_b) {
                file_name = paste0(dir,'mcmc_out_',toString(seed), '_', trial_num, '_30b.rda')   
            } else {
                file_name = paste0(dir,'mcmc_out_',toString(seed), '_', trial_num, '_30.rda')      
            }
        } else {
            if(case_b) {
                file_name = paste0(dir,'mcmc_out_',toString(seed), '_', trial_num, '_15b.rda')   
            } else {
                file_name = paste0(dir,'mcmc_out_',toString(seed), '_', trial_num, '_15.rda')      
            }
        }
    }
    
    if (file.exists(file_name)) {
        load(file_name)
        ind = ind + 1

        print(mcmc_out$accept)
	    # mcmc_out$chain[,c(par_index$tau2, par_index$sigma2)] = 
	    #     exp(mcmc_out$chain[,c(par_index$tau2, par_index$sigma2)])

        # Thinning the chain
        main_chain = mcmc_out$chain[index_post,]
        ind_keep = seq(1, nrow(main_chain), by=100)
        tau_sig_sum = main_chain[,par_index$tau2] + main_chain[,par_index$sigma2[1]]
        main_chain = cbind(main_chain, tau_sig_sum)

      	chain_list[[ind]] = main_chain[ind_keep, ]
    	post_means[ind,] <- colMeans(main_chain[ind_keep, ])
    }
}

# Plot and save the mcmc trace plots and histograms.
pdf_title = NULL
if(simulation) {
    if(thirty) {
        if(case_b) {
            pdf_title = paste0('Plots/mcmc_out_', trial_num, '_sim_30b.pdf')
        } else {
            pdf_title = paste0('Plots/mcmc_out_', trial_num, '_sim_30.pdf')   
        }
    } else {
        if(case_b) {
            pdf_title = paste0('Plots/mcmc_out_', trial_num, '_sim_15b.pdf')
        } else {
            pdf_title = paste0('Plots/mcmc_out_', trial_num, '_sim_15.pdf')   
        }
    }
} else {
    if(thirty) {
        if(case_b) {
            pdf_title = paste0('Plots/mcmc_out_', trial_num, '_30b.pdf')
        } else {
            pdf_title = paste0('Plots/mcmc_out_', trial_num, '_30.pdf')   
        }
    } else {
        if(case_b) {
            pdf_title = paste0('Plots/mcmc_out_', trial_num, '_15b.pdf')
        } else {
            pdf_title = paste0('Plots/mcmc_out_', trial_num, '_15.pdf')   
        }
    }
}
print(pdf_title)

pdf(pdf_title)
par(mfrow=c(4, 2))

stacked_chains = do.call( rbind, chain_list)
par_mean = par_median = upper = lower = rep( NA, length(labels))

if(thirty) {
    load('Data/true_par_2_30.rda')
} else {
    load('Data/true_par_15.rda')
}

# tau2_hat   = 0.3395152
# sigma2_hat = 1.227959
mle_ind = 1

for(r in 1:length(labels)){
    
    par_mean[r] = round( mean(stacked_chains[,r]), 4)
    par_median[r] = round( median(stacked_chains[,r]), 4)
    upper[r] = round( quantile( stacked_chains[,r], prob=.975), 4)
    lower[r] = round( quantile( stacked_chains[,r], prob=.025), 4)

    if (simulation) {
        plot( NULL, xlab=paste0("true val: ", round(true_par[r], 3)), ylab=NA, main=labels[r], xlim=c(1,nrow(chain_list[[1]])),
              ylim=range(stacked_chains[,r]) )
    } else {
        
        plot( NULL, xlab=paste0("[", lower[r], ", ", upper[r], "]"), ylab=NA, main=labels[r], xlim=c(1,nrow(chain_list[[1]])),
              ylim=range(stacked_chains[,r]) )
    }

    for(seed in 1:length(chain_list)) lines( chain_list[[seed]][,r], type='l', col=seed)

    hist( stacked_chains[,r], breaks=sqrt(nrow(stacked_chains)), ylab=NA, main=NA,
            freq=F, xlab=paste0('Mean = ',toString(par_mean[r]),
                                ' Median = ',toString(par_median[r])))
    abline( v=upper[r], col='red', lwd=2, lty=2)
    abline( v=lower[r], col='purple', lwd=2, lty=2)
    
    if(simulation) {
        abline( v=true_par[r], col='green', lwd=2, lty=2)
    } else {
        mle_vals = c(6.4557765, -0.2582171, -0.1166778, 
                    -1.0842574,  0.1973717, -0.0714276, 0.2777276)
        if(r %in% c(par_index$delta, par_index$tau2, par_index$sigma2)) {
            abline( v=mle_vals[mle_ind], col='blue', lwd=2, lty=2)
            mle_ind = mle_ind + 1
        }

        if(r == length(labels)) abline( v=log(1.227959) + log(0.3395152),  col='blue', lwd=2, lty=2)
    }
}

dev.off()
