# This script file produces trace plots and histograms of the mcmc output files
library(tidyverse)
library(gridExtra)

args = commandArgs(TRUE)
deSolve_or_expm = as.numeric(args[1])  # 1: deSolve, 2: expm

dir = 'real_ecog_analysis/Model_out/' 
model_name = c('deSolve', 'expm')

# Size of posterior sample from mcmc chains
n_post = 15000
# Step number at which the adaptive tuning scheme was frozen
burnin = 5000
# Total number of steps the mcmc algorithm is computed for
steps = 30000
# Matrix row indices for the posterior sample to use
index_post = (steps - burnin - n_post + 1):(steps - burnin)

par_index = list( beta=1:24, misclass = 25:30, pi_logit=31:33, l_delta = 34:37, 
                  l_theta=38:41, l_alpha=42:45, l_beta=46:49)

index_seeds = c(1:10)

labels <- c("Baseline: IS --> NREM", "Baseline: IS --> REM", "Baseline: IS --> LIMBO",   
            "Baseline: NREM --> IS", "Baseline: NREM --> REM", "Baseline: NREM --> LIMBO",
            "Baseline: REM --> IS", "Baseline: REM --> NREM", "Baseline: REM --> LIMBO",
            "Baseline: LIMBO --> IS", "Baseline: LIMBO --> NREM", "Baseline: LIMBO --> REM",
            "Time: IS --> NREM", "Time: IS --> REM", "Time: IS --> LIMBO",      
            "Time: NREM --> IS",   "Time: NREM --> REM", "Time: NREM --> LIMBO",
            "Time: REM --> IS", "Time: REM --> NREM", "Time: REM --> LIMBO",
            "Time: LIMBO --> IS", "Time: LIMBO --> NREM", "Time: LIMBO --> REM",
            "logit P( obs. NREM | true IS )", "logit P( obs. REM | true IS )",
            "logit P( obs. IS | true NREM )", "logit P( obs. REM | true NREM )", 
            "logit P( obs. IS | true REM )", "logit P( obs. NREM | true REM )", 
            "logit P( init NREM )", "logit P( init REM )", "logit P( init LIMBO )",
            "Delta (IS)", "Delta (NREM)", "Delta (REM)", "Delta (LIMBO)",
            "Theta (IS)", "Theta (NREM)", "Theta (REM)", "Theta (LIMBO)",
            "Alpha (IS)", "Alpha (NREM)", "Alpha (REM)", "Alpha (LIMBO)",
            "Beta (IS)", "Beta (NREM)", "Beta (REM)", "Beta (LIMBO)")

# -----------------------------------------------------------------------------
# Create mcmc trace plots and histograms
# -----------------------------------------------------------------------------
chain_list = vector(mode = "list", length = length(index_seeds))
post_means = matrix(nrow = length(index_seeds), ncol = length(labels))

ind = 0

for(seed in index_seeds){
    file_name = NULL
    if(deSolve_or_expm == 1) {
        file_name = paste0(dir,'mcmc_out_',toString(seed),'.rda')
    } else {
        file_name = paste0(dir,'mcmc_out_',toString(seed),'_expm.rda')
    }
    
    if (file.exists(file_name)) {
        load(file_name)
        ind = ind + 1

        # Thinning the chain
        main_chain = mcmc_out$chain[index_post,]
        ind_keep = seq(1, nrow(main_chain), by=10)

      	chain_list[[ind]] = main_chain[ind_keep, ]
    	post_means[ind,] <- colMeans(main_chain[ind_keep, ])
    }
}

# Plot and save the mcmc trace plots and histograms.

pdf(paste0('real_ecog_analysis/Plots/mcmc_', model_name[deSolve_or_expm], '.pdf'))
par(mfrow=c(4, 2))

stacked_chains = do.call( rbind, chain_list)
par_mean = par_median = upper = lower = rep( NA, length(labels))

for(r in 1:length(labels)){

    plot( NULL, xlab=NA, ylab=NA, main=labels[r], xlim=c(1,nrow(chain_list[[1]])),
            ylim=range(stacked_chains[,r]) )

    for(seed in 1:length(chain_list)) lines( chain_list[[seed]][,r], type='l', col=seed)

    par_mean[r] = round( mean(stacked_chains[,r]), 4)
    par_median[r] = round( median(stacked_chains[,r]), 4)
    upper[r] = round( quantile( stacked_chains[,r], prob=.975), 4)
    lower[r] = round( quantile( stacked_chains[,r], prob=.025), 4)
    print(paste0(r, ". ", labels[r],": [", lower[r], ", ", upper[r], "]"))

    hist( stacked_chains[,r], breaks=sqrt(nrow(stacked_chains)), ylab=NA, main=NA,
            freq=F, xlab=paste0('Mean = ',toString(par_mean[r]),
                                ' Median = ',toString(par_median[r])))
    abline( v=upper[r], col='red', lwd=2, lty=2)
    abline( v=lower[r], col='purple', lwd=2, lty=2)

}
cred_set_cumulative = cbind(lower, upper)
save(cred_set_cumulative, file = paste0('real_ecog_analysis/Plots/cred_set_cumulative_', 
                                            model_name[deSolve_or_expm], '.rda'))

dev.off()
