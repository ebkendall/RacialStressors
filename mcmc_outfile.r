# This script file produces trace plots and histograms of the mcmc output files
library(tidyverse)
library(gridExtra)
library(latex2exp)

dir = 'Model_out/' 
args <- commandArgs(TRUE)
trial_num = as.numeric(args[1])
simulation = as.numeric(args[2])

# Size of posterior sample from mcmc chains
n_post = 5000
# Step number at which the adaptive tuning scheme was frozen
burnin = 5000
# Total number of steps the mcmc algorithm is computed for
steps = 20000
# Matrix row indices for the posterior sample to use
index_post = (steps - burnin - n_post + 1):(steps - burnin)

# par_index = list( zeta=1:10, misclass=11:16,
#                   delta = 17:19, tau2 = 20, upsilon = 21:29,
#                   delta_i = 30:302)
# par_index = list( zeta=1:5, misclass=6:11,
#                   delta = 12:14, tau2 = 15, upsilon = 16:24,
#                   delta_i = 25:297)
par_index = list( zeta=1:8, misclass=9:14,
                  delta = 15:17, tau2 = 18, upsilon = 19:27,
                  delta_i = 28:300)
index_seeds = c(1)

labels <- c(TeX(r'($\hat{\zeta}_{0,1}:$ Baseline: 1 $\to$ 2)'),  
            # TeX(r'($\hat{\zeta}_{0,3}:$ Baseline: 2 $\to$ 1)'), 
            TeX(r'($\hat{\zeta}_{0,2}:$ Baseline: 2 $\to$ 3)'), 
            TeX(r'($\hat{\zeta}_{0,3}:$ Baseline: 3 $\to$ 1)'),
            TeX(r'($\hat{\zeta}_{0,4}:$ Baseline: 3 $\to$ 2)'),
            TeX(r'($\hat{\zeta}_{1,1}:$ Time: 1 $\to$ 2)'), 
            # TeX(r'($\hat{\zeta}_{1,3}:$ Time: 2 $\to$ 1)'), 
            TeX(r'($\hat{\zeta}_{1,2}:$ Time: 2 $\to$ 3)'), 
            TeX(r'($\hat{\zeta}_{1,3}:$ Time: 3 $\to$ 1)'), 
            TeX(r'($\hat{\zeta}_{1,4}:$ Time: 3 $\to$ 2)'),
            TeX(r'(P(obs. S2 | true S1))'), TeX(r'(P(obs. S3 | true S1))'), TeX(r'(P(obs. S1 | true S2))'),
            TeX(r'(P(obs. S3 | true S2))'), TeX(r'(P(obs. S1 | true S3))'), TeX(r'(P(obs. S2 | true S3))'),
            TeX(r'($\delta_1 = \mu$)'), TeX(r'($\delta_2 = \alpha$)'), TeX(r'($\delta_3 = \beta$)'),
            TeX(r'($\tau^2$)'), 
            TeX(r'($\Upsilon_{1,1}$)'), TeX(r'($\Upsilon_{2,1}$)'), TeX(r'($\Upsilon_{3,1}$)'), 
            TeX(r'($\Upsilon_{1,2}$)'), TeX(r'($\Upsilon_{2,2}$)'), TeX(r'($\Upsilon_{3,2}$)'),
            TeX(r'($\Upsilon_{1,3}$)'), TeX(r'($\Upsilon_{2,3}$)'), TeX(r'($\Upsilon_{3,3}$)'))

# -----------------------------------------------------------------------------
# Create mcmc trace plots and histograms
# -----------------------------------------------------------------------------
chain_list = vector(mode = "list", length = length(index_seeds))
post_means = matrix(nrow = length(index_seeds), ncol = length(labels))

ind = 0

for(seed in index_seeds){

    file_name = paste0(dir,'mcmc_out_',toString(seed), '_', trial_num, '.rda')
    
    if (file.exists(file_name)) {
        load(file_name)
        ind = ind + 1

        print(mcmc_out$accept)

        # Thinning the chain
        main_chain = mcmc_out$chain[index_post,]
        ind_keep = seq(1, nrow(main_chain), by=10)

      	chain_list[[ind]] = main_chain[ind_keep, ]
    	post_means[ind,] <- colMeans(main_chain[ind_keep, ])
    }
}

# Plot and save the mcmc trace plots and histograms.

pdf(paste0('Plots/mcmc_out_', trial_num, '.pdf'))
par(mfrow=c(4, 2))

stacked_chains = do.call( rbind, chain_list)
par_mean = par_median = upper = lower = rep( NA, length(labels))

labels_sub <- 1:length(labels)

load('Data/Simulation/true_par.rda')

for(r in 1:length(labels_sub)){

    if (simulation == 1) {
        plot( NULL, xlab=paste0("true val: ", round(true_par[r], 3)), ylab=NA, main=labels[r], xlim=c(1,nrow(chain_list[[1]])),
              ylim=range(stacked_chains[,r]) )
    } else {
        plot( NULL, xlab=NA, ylab=NA, main=labels[r], xlim=c(1,nrow(chain_list[[1]])),
              ylim=range(stacked_chains[,r]) )
    }

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
    
    if(simulation == 1) {
        abline( v=true_par[r], col='green', lwd=2, lty=2)
    }

}

dev.off()

print("zeta")
print(matrix(par_mean[par_index$zeta],ncol=2))

print("misclass mean")
resp_fnc = matrix(c(1, exp(par_mean[par_index$misclass][1]), exp(par_mean[par_index$misclass][2]),
                      exp(par_mean[par_index$misclass][3]), 1, exp(par_mean[par_index$misclass][4]),
                      exp(par_mean[par_index$misclass][5]), exp(par_mean[par_index$misclass][6]),1), 
                      ncol=3, byrow=TRUE)
resp_fnc = resp_fnc / rowSums(resp_fnc)
print(resp_fnc)

print("delta")
print(par_mean[par_index$delta])

print('tau2')
print(par_mean[par_index$tau2])

# print('upsilon')
# print(matrix(par_mean[par_index$upsilon], ncol = 3))

# print("zeta mean")
# print(matrix( par_mean[par_index$zeta], ncol=2))

# print("a mean")
# print(exp(par_mean[par_index$log_a_vec]))

# print("c mean")
# print(exp(par_mean[par_index$log_c]))
