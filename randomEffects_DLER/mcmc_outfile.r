# This script file produces trace plots and histograms of the mcmc output files
library(latex2exp)
library(tidyverse)
library(gridExtra)
library(egg)

dir = 'Model_out/' 

# Model type -------------------------------------------------------------------
# 1: baseline model (age, sex, pEdu)
# 2: DLER
# 3: all covariates

covariate_struct = 3
# ------------------------------------------------------------------------------

# Information defining which approach to take ----------------------------------
simulation = T
case_b = T

if(simulation) {
    trial_num = covariate_struct
    index_seeds = c(1:100)
} else {
    trial_num = covariate_struct+3
    index_seeds = c(1:5)
}

# ------------------------------------------------------------------------------

# Size of posterior sample from mcmc chains
n_post = 19500; burnin = 500; steps = 20000
# n_post = 20000; burnin = 500; steps = 100000

# Matrix row indices for the posterior sample to use
index_post = (steps - burnin - n_post + 1):(steps - burnin)

if(covariate_struct == 1) {
    par_index = list(zeta=1:24, misclass=38:41, delta = 25:27, tau2 = 28, 
                     sigma2 = 29:31, gamma = 32:34, delta_new = 35:37)

    labels <- c(TeX(r'($\hat{\zeta}_{0,1}:$ baseline: 1 $\to$ 2)'), 
                TeX(r'($\hat{\zeta}_{0,2}:$ baseline: 1 $\to$ 3)'),
                TeX(r'($\hat{\zeta}_{0,3}:$ baseline: 2 $\to$ 1)'),
                TeX(r'($\hat{\zeta}_{0,4}:$ baseline: 2 $\to$ 3)'),
                TeX(r'($\hat{\zeta}_{0,5}:$ baseline: 3 $\to$ 1)'),
                TeX(r'($\hat{\zeta}_{0,6}:$ baseline: 3 $\to$ 2)'),
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
                TeX(r'($\mu$)'), TeX(r'($\alpha$)'), TeX(r'($\beta$)'),
                TeX(r'($\log(\tau^2)$)'), TeX(r'($\log(\sigma_1^2)$)'), 
                TeX(r'($\log(\sigma_2^2)$)'), TeX(r'($\log(\sigma_3^2)$)'),
                TeX(r'($\hat{\gamma}_1:$ age)'), TeX(r'($\hat{\gamma}_2:$ sex1)'), 
                TeX(r'($\hat{\gamma}_3:$ yes edu)'),
                TeX(r'($\hat{\delta}_1:$ DLER baseline)'), 
                TeX(r'($\hat{\delta}_2:$ DLER state 2)'), 
                TeX(r'($\hat{\delta}_3:$ DLER state 3)'),
                TeX(r'(logit P(obs S2 | true S1))'), TeX(r'(logit P(obs S3 | true S1))'),
                TeX(r'(logit P(obs S3 | true S2))'), TeX(r'(logit P(obs S2 | true S3))'),
                TeX(r'($\mu + \alpha$)'), TeX(r'($\mu + \beta$)'),
                TeX(r'($\tau^2 + \sigma_1^2$)'), TeX(r'($\tau^2 + \sigma_2^2$)'), 
                TeX(r'($\tau^2 + \sigma_3^2$)'), TeX(r'($\tau^2$)'),
                TeX(r'($\sigma_1^2$)'), TeX(r'($\sigma_2^2$)'), 
                TeX(r'($\sigma_3^2$)'))
} else if(covariate_struct == 2) {
    par_index = list(zeta=1:12, misclass=26:29, delta = 13:15, tau2 = 16, 
                     sigma2 = 17:19, gamma = 20:22, delta_new = 23:25)

    labels <- c(TeX(r'($\hat{\zeta}_{0,1}:$ baseline: 1 $\to$ 2)'), 
                TeX(r'($\hat{\zeta}_{0,2}:$ baseline: 1 $\to$ 3)'),
                TeX(r'($\hat{\zeta}_{0,3}:$ baseline: 2 $\to$ 1)'),
                TeX(r'($\hat{\zeta}_{0,4}:$ baseline: 2 $\to$ 3)'),
                TeX(r'($\hat{\zeta}_{0,5}:$ baseline: 3 $\to$ 1)'),
                TeX(r'($\hat{\zeta}_{0,6}:$ baseline: 3 $\to$ 2)'),
                TeX(r'($\hat{\zeta}_{1,1}:$ DLER: 1 $\to$ 2)'), 
                TeX(r'($\hat{\zeta}_{1,2}:$ DLER: 1 $\to$ 3)'),
                TeX(r'($\hat{\zeta}_{1,3}:$ DLER: 2 $\to$ 1)'),
                TeX(r'($\hat{\zeta}_{1,4}:$ DLER: 2 $\to$ 3)'),
                TeX(r'($\hat{\zeta}_{1,5}:$ DLER: 3 $\to$ 1)'),
                TeX(r'($\hat{\zeta}_{1,6}:$ DLER: 3 $\to$ 2)'),
                TeX(r'($\mu$)'), TeX(r'($\alpha$)'), TeX(r'($\beta$)'),
                TeX(r'($\log(\tau^2)$)'), TeX(r'($\log(\sigma_1^2)$)'), 
                TeX(r'($\log(\sigma_2^2)$)'), TeX(r'($\log(\sigma_3^2)$)'),
                TeX(r'($\hat{\gamma}_1:$ age)'), TeX(r'($\hat{\gamma}_2:$ sex1)'), 
                TeX(r'($\hat{\gamma}_3:$ yes edu)'),
                TeX(r'($\hat{\delta}_1:$ DLER baseline)'), 
                TeX(r'($\hat{\delta}_2:$ DLER state 2)'), 
                TeX(r'($\hat{\delta}_3:$ DLER state 3)'),
                TeX(r'(logit P(obs S2 | true S1))'), TeX(r'(logit P(obs S3 | true S1))'),
                TeX(r'(logit P(obs S3 | true S2))'), TeX(r'(logit P(obs S2 | true S3))'),
                TeX(r'($\mu + \alpha$)'), TeX(r'($\mu + \beta$)'),
                TeX(r'($\tau^2 + \sigma_1^2$)'), TeX(r'($\tau^2 + \sigma_2^2$)'), 
                TeX(r'($\tau^2 + \sigma_3^2$)'), TeX(r'($\tau^2$)'),
                TeX(r'($\sigma_1^2$)'), TeX(r'($\sigma_2^2$)'), 
                TeX(r'($\sigma_3^2$)'))
} else {
    par_index = list(zeta=1:30, misclass=44:47, delta = 31:33, tau2 = 34, 
                     sigma2 = 35:37, gamma = 38:40, delta_new = 41:43)

    labels <- c(TeX(r'($\hat{\zeta}_{0,1}:$ baseline: 1 $\to$ 2)'), 
                TeX(r'($\hat{\zeta}_{0,2}:$ baseline: 1 $\to$ 3)'),
                TeX(r'($\hat{\zeta}_{0,3}:$ baseline: 2 $\to$ 1)'),
                TeX(r'($\hat{\zeta}_{0,4}:$ baseline: 2 $\to$ 3)'),
                TeX(r'($\hat{\zeta}_{0,5}:$ baseline: 3 $\to$ 1)'),
                TeX(r'($\hat{\zeta}_{0,6}:$ baseline: 3 $\to$ 2)'),
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
                TeX(r'($\mu$)'), TeX(r'($\alpha$)'), TeX(r'($\beta$)'),
                TeX(r'($\log(\tau^2)$)'), TeX(r'($\log(\sigma_1^2)$)'), 
                TeX(r'($\log(\sigma_2^2)$)'), TeX(r'($\log(\sigma_3^2)$)'),
                TeX(r'($\hat{\gamma}_1:$ age)'), TeX(r'($\hat{\gamma}_2:$ sex1)'), 
                TeX(r'($\hat{\gamma}_3:$ yes edu)'),
                TeX(r'($\hat{\delta}_1:$ DLER baseline)'), 
                TeX(r'($\hat{\delta}_2:$ DLER state 2)'), 
                TeX(r'($\hat{\delta}_3:$ DLER state 3)'),
                TeX(r'(logit P(obs S2 | true S1))'), TeX(r'(logit P(obs S3 | true S1))'),
                TeX(r'(logit P(obs S3 | true S2))'), TeX(r'(logit P(obs S2 | true S3))'),
                TeX(r'($\mu + \alpha$)'), TeX(r'($\mu + \beta$)'),
                TeX(r'($\tau^2 + \sigma_1^2$)'), TeX(r'($\tau^2 + \sigma_2^2$)'), 
                TeX(r'($\tau^2 + \sigma_3^2$)'), TeX(r'($\tau^2$)'),
                TeX(r'($\sigma_1^2$)'), TeX(r'($\sigma_2^2$)'), 
                TeX(r'($\sigma_3^2$)'))
}


# -----------------------------------------------------------------------------
# Create mcmc trace plots and histograms
# -----------------------------------------------------------------------------
chain_list = vector(mode = "list", length = length(index_seeds))
post_means = matrix(nrow = length(index_seeds), ncol = length(labels))

ind = 0

cred_set = vector(mode = 'list', length = length(labels))
for(i in 1:length(cred_set)) cred_set[[i]] = matrix(nrow = 100, ncol = 2)

for(seed in index_seeds){
    
    file_name = NULL
    
    if(simulation) {
        if(case_b) {
            file_name = paste0(dir,'mcmc_out_',toString(seed), '_', trial_num, '_sim_30b_check.rda')  
        } else {
            file_name = paste0(dir,'mcmc_out_',toString(seed), '_', trial_num, '_sim_30.rda')
        }
    } else {
        if(case_b) {
            file_name = paste0(dir,'mcmc_out_',toString(seed), '_', trial_num, '_30b_check.rda')
        } else {
            file_name = paste0(dir,'mcmc_out_',toString(seed), '_', trial_num, '_30.rda')
        }
    }
    
    if (file.exists(file_name)) {
        load(file_name)
        ind = ind + 1

        print(file_name)

        # Thinning the chain
        main_chain = mcmc_out$chain[index_post,]
        ind_keep = seq(1, nrow(main_chain), by=10)
        # ind_keep = seq(1, nrow(main_chain), by=1)
        
        mu_alpha_sum = main_chain[,par_index$delta[1]] + main_chain[,par_index$delta[2]]
        mu_beta_sum = main_chain[,par_index$delta[1]] + main_chain[,par_index$delta[3]]
        
        tau_sig1_sum = exp(main_chain[,par_index$tau2]) + exp(main_chain[,par_index$sigma2[1]])
        tau_sig2_sum = exp(main_chain[,par_index$tau2]) + exp(main_chain[,par_index$sigma2[2]])
        tau_sig3_sum = exp(main_chain[,par_index$tau2]) + exp(main_chain[,par_index$sigma2[3]])
        
        tau2 = exp(main_chain[,par_index$tau2])
        sig1 = exp(main_chain[,par_index$sigma2[1]])
        sig2 = exp(main_chain[,par_index$sigma2[2]])
        sig3 = exp(main_chain[,par_index$sigma2[3]])

        main_chain = cbind(main_chain, mu_alpha_sum, mu_beta_sum, 
                           tau_sig1_sum, tau_sig2_sum, tau_sig3_sum,tau2,
                           sig1, sig2, sig3)

        chain_list[[ind]] = main_chain[ind_keep, ]
    	post_means[ind,] <- colMeans(main_chain[ind_keep, ])
        
        if(simulation) {
            for(j in 1:ncol(main_chain)) {
                cred_set[[j]][seed, 1] = round(quantile( main_chain[ind_keep,j],
                                                    prob=.025), 4)
                cred_set[[j]][seed, 2] = round(quantile( main_chain[ind_keep,j],
                                                    prob=.975), 4)
            }
        } 
    }
}

# True parameter values if the simulation
if(simulation) {
    load(paste0('Data/true_par_', covariate_struct, '_30.rda'))
    mu_alpha_sum = true_par[par_index$delta[1]] + true_par[par_index$delta[2]]
    mu_beta_sum = true_par[par_index$delta[1]] + true_par[par_index$delta[3]]
    
    tau_sig1_sum = exp(true_par[par_index$tau2]) + exp(true_par[par_index$sigma2[1]])
    tau_sig2_sum = exp(true_par[par_index$tau2]) + exp(true_par[par_index$sigma2[2]])
    tau_sig3_sum = exp(true_par[par_index$tau2]) + exp(true_par[par_index$sigma2[3]])
    
    tau2 = exp(true_par[par_index$tau2])
    sig1 = exp(true_par[par_index$sigma2[1]])
    sig2 = exp(true_par[par_index$sigma2[2]])
    sig3 = exp(true_par[par_index$sigma2[3]])
    true_par = c(true_par, mu_alpha_sum, mu_beta_sum, 
                 tau_sig1_sum, tau_sig2_sum, tau_sig3_sum,tau2,
                 sig1, sig2, sig3)
} 


# Calculate the coverage in the simulation
if(simulation) {
    cov_df = rep(NA, length(labels))
    for(i in 1:length(labels)) {
        val = true_par[i]
        cov_df[i] = mean(cred_set[[i]][,1] <= val & val <= cred_set[[i]][,2], na.rm=T)
    }
}

# Plot and save the mcmc trace plots and histograms.
pdf_title = NULL
if(simulation) {
    if(case_b) {
        pdf_title = paste0('Plots/mcmc_out_', trial_num, '_sim_30b_check.pdf')
    } else {
        pdf_title = paste0('Plots/mcmc_out_', trial_num, '_sim_30.pdf')   
    }
} else {
    if(case_b) {
        pdf_title = paste0('Plots/mcmc_out_', trial_num, '_30b_check.pdf')
    } else {
        pdf_title = paste0('Plots/mcmc_out_', trial_num, '_30.pdf')   
    }
}
print(pdf_title)

pdf(pdf_title)
par(mfrow=c(4, 2))

stacked_chains = do.call( rbind, chain_list)
par_mean = par_median = upper = lower = max_confidence_set = rep( NA, length(labels))

mle_ind = 1

for(r in 1:length(labels)){
    
    par_mean[r] = round( mean(stacked_chains[,r]), 4)
    par_median[r] = round( median(stacked_chains[,r]), 4)
    upper[r] = round( quantile( stacked_chains[,r], prob=.975), 4)
    lower[r] = round( quantile( stacked_chains[,r], prob=.025), 4)

    # Determining the max confidence set
    for(ci in seq(0.001, 0.499, by = 0.001)) {
        upp_ci = quantile( stacked_chains[,r], prob=1-ci)
        low_ci = quantile( stacked_chains[,r], prob=ci)
        if(0 < low_ci | 0 > upp_ci) {
            max_confidence_set[r] = round(1 - ci - ci, 3)
            break
        }
    }

    if(simulation) {
        plot( NULL, xlab=paste0("true val: ", round(true_par[r], 3)), ylab=NA, 
              main=labels[r], xlim=c(1,nrow(chain_list[[1]])),
              ylim=range(stacked_chains[,r]) )
    } else {
        
        plot( NULL, xlab=paste0("[", lower[r], ", ", upper[r], 
                                "], max conf. = ", 100*max_confidence_set[r], "%"), 
              ylab=NA, main=labels[r], xlim=c(1,nrow(chain_list[[1]])),
              ylim=range(stacked_chains[,r]) )
    }

    for(seed in 1:length(chain_list)) lines( chain_list[[seed]][,r], type='l', col=seed)

    hist( stacked_chains[,r], breaks=sqrt(nrow(stacked_chains)), ylab=NA, main=NA,
            freq=F, xlab=paste0('Mean = ',toString(par_mean[r]),
                                ' Median = ',toString(par_median[r])))
    abline( v=upper[r], col='red', lwd=2, lty=2)
    abline( v=lower[r], col='purple', lwd=2, lty=2)
    
    if(!simulation) {
        if(r == par_index$delta[1]) abline( v=  6.464088, col='blue', lwd=2, lty=2)
        if(r == par_index$delta[2]) abline( v= -0.2681087, col='blue', lwd=2, lty=2)
        if(r == par_index$delta[3]) abline( v= -0.1132974, col='blue', lwd=2, lty=2)   
    }
    
    if(simulation) {
        abline( v=true_par[r], col='green', lwd=2, lty=2)
    } else {
        mle_val = c(6.1959793814433, 6.35079064587973,
                    1.55661033186978, 1.32762953909378, 1.57489792721762)
        if(r > max(par_index$misclass)) {
            abline( v=mle_val[mle_ind], col='blue', lwd=2, lty=2)
            mle_ind = mle_ind + 1
        }
    }
}

VP = vector(mode = 'list', length = length(labels))
if(simulation) {
    for(r in 1:length(labels)) {
        # Adding the boxplots
        yVar = post_means[,r]
        x_label = paste0("Coverage is: ", round(cov_df[r], digits=3))

        plot_df = data.frame(yVar = yVar, disc_type = covariate_struct)
        VP[[r]] = ggplot(plot_df, aes(x=disc_type, y = yVar)) +
        geom_violin(trim=FALSE) +
        geom_boxplot(width=0.1) +
        ggtitle(labels[r]) +
        ylab(paste0("Parameter Value: ", round(true_par[r], 3))) +
        xlab(x_label) +
        geom_hline(yintercept=true_par[r], linetype="dashed", color = "red") +
        theme(text = element_text(size = 7))
    }

    grid.arrange(VP[[1]], VP[[2]], VP[[3]], VP[[4]], VP[[5]],
                VP[[6]], VP[[7]], VP[[8]], VP[[9]], ncol=3, nrow =3)
    grid.arrange(VP[[10]], VP[[11]], VP[[12]], VP[[13]], VP[[14]],
                VP[[15]], VP[[16]], VP[[17]], VP[[18]], ncol=3, nrow =3)
    if(covariate_struct > 1) {
        grid.arrange(VP[[19]], VP[[20]], VP[[21]], VP[[22]], VP[[23]],
                     VP[[24]], VP[[25]], VP[[26]], VP[[27]], ncol=3, nrow =3)
        grid.arrange(VP[[28]], VP[[29]], VP[[30]], VP[[31]], VP[[32]],
                     VP[[33]], VP[[34]], VP[[35]], VP[[36]], ncol=3, nrow =3)
        grid.arrange(VP[[37]], VP[[38]], VP[[39]], VP[[40]], VP[[41]], 
                     VP[[42]], VP[[43]], VP[[44]], VP[[45]], ncol=3, nrow =3) 
        grid.arrange(VP[[46]], VP[[47]], VP[[48]], VP[[49]], VP[[50]], 
                     VP[[51]], VP[[52]], VP[[53]], VP[[54]], ncol=3, nrow =3)    
    }
}
dev.off()


# Calculate the 95% credible sets centered at the post. median of log(sigma_123)
log_sigmas = stacked_chains[,par_index$sigma2]
interval_width = matrix(nrow = 3, ncol = 2)
for(i in 1:ncol(log_sigmas)) {
    print(paste0("sigma", i))
    med_i = median(log_sigmas[,i])
    
    for(j in seq(0, 10, by = 0.00001)) {
        sum_j = sum((log_sigmas[,i] <= med_i + j) & (log_sigmas[,i] >= med_i - j))
        if(sum_j >= .95 * nrow(log_sigmas)) {
            interval_width[i, 1] = j
            interval_width[i, 2] = sum_j
            break
        }
    }
}

print("95% Credible Sets for log(sigma)")
print(cbind(apply(log_sigmas, 2, median) - interval_width[,1],
            apply(log_sigmas, 2, median) + interval_width[,1]))
print("95% Credible width")
print(interval_width)

if(!simulation) {
    # Probability of transitioning in 30s with certain covariate combinations -----
    load('../Data/data_format_30.rda')
    data_format = data_format_30
    
    # Removing the participants with missing labels
    miss_info = c(26296, 29698, 30625, 401, 423, 419, 457)
    data_format = data_format[!(data_format[,"ID.."] %in% miss_info), ]
    
    ages = NULL
    dler_val = NULL
    for(a in unique(data_format[,"ID.."])) {
        ages = c(ages, unique(data_format[data_format[,"ID.."] == a, "Age"]))
        dler_val = c(dler_val, unique(data_format[data_format[,"ID.."] == a, "DLER_avg"]))
    }
    s_dler = sd(dler_val)
    print(paste0("Mean DLER: ", mean(dler_val)))
    print(paste0("Standard dev. DLER: ", s_dler))
    
    sex = c(0,1)
    pEdu = c(0,1)
    dler = c(-s_dler,0,s_dler)
    
    if(covariate_struct == 1) {
        zeta_est = matrix(par_median[par_index$zeta], nrow=6, ncol=4)
    } else if(covariate_struct == 2) {
        zeta_est = matrix(par_median[par_index$zeta], nrow=6, ncol=2)
    } else {
        zeta_est = matrix(par_median[par_index$zeta], nrow=6, ncol=5)   
    }
    
    # Probability evolution curves ------------------------------------------------
    prob_evo <- function( zeta_est, t_i, cov_val){
        
        # Compute prob evolutions of transitioning 1->2, 1->3, 2->3
        prob_evo_mat = matrix( 0, 3, length(t_i))
        P_t = matrix(c(1,0,0), ncol = 3)
        prob_evo_mat[,1] = c(P_t)
        
        for(k in 2:length(t_i)){
            
            qs = zeta_est %*% cov_val
            Q = matrix(c( 1,          exp(qs[1]), exp(qs[2]),
                          exp(qs[3]),          1, exp(qs[4]),
                          exp(qs[5]), exp(qs[6]),          1), ncol = 3, byrow=T)
            P = Q / rowSums(Q)
            
            P_t = P_t %*% P
            prob_evo_mat[,k] = c(P_t)
        }
        
        return(prob_evo_mat)
    }
    
    # baseline, age, sex, pEdu, dler
    dler_vals = seq(2.6, -1.4, by = -0.2)
    t_i = 1:200
    
    # Sex = 0
    case_0 = vector(mode = 'list', length = length(dler_vals))
    # Sex = 1
    case_1 = vector(mode = 'list', length = length(dler_vals))
    for(i in 1:length(dler_vals)) {
        if(covariate_struct == 1) {
            cov_val_0 = matrix(c(1, 0, 0, 0), ncol = 1)
            cov_val_1 = matrix(c(1, 0, 1, 0), ncol = 1)
        } else if(covariate_struct == 2) {
            cov_val_0 = matrix(c(1, dler_vals[i]), ncol = 1)
            cov_val_1 = matrix(c(1, dler_vals[i]), ncol = 1)
        } else {
            cov_val_0 = matrix(c(1, 0, 0, 0, dler_vals[i]), ncol = 1)
            cov_val_1 = matrix(c(1, 0, 1, 0, dler_vals[i]), ncol = 1)
        }
        
        case_0[[i]] = prob_evo(zeta_est, t_i, cov_val_0)
        case_1[[i]] = prob_evo(zeta_est, t_i, cov_val_1)
    }

    if(covariate_struct == 1) {
        qs = zeta_est %*% cov_val_0

        Q = matrix(c( 1,          exp(qs[1]), exp(qs[2]),
                        exp(qs[3]),          1, exp(qs[4]),
                        exp(qs[5]), exp(qs[6]),          1), ncol = 3, byrow=T)
        P = Q / rowSums(Q)
        print("Transition Probability matrix (sex = 0): ")
        print(P)

        qs = zeta_est %*% cov_val_1

        Q = matrix(c( 1,          exp(qs[1]), exp(qs[2]),
                        exp(qs[3]),          1, exp(qs[4]),
                        exp(qs[5]), exp(qs[6]),          1), ncol = 3, byrow=T)
        P = Q / rowSums(Q)
        print("Transition Probability matrix (sex = 1): ")
        print(P)
    }

    if(covariate_struct != 1) {
        
        prob0_1_to_1 = cbind(t_i, case_0[[1]][1,], dler_vals[1] + mean(dler_val))
        prob1_1_to_1 = cbind(t_i, case_1[[1]][1,], dler_vals[1] + mean(dler_val))
        prob0_1_to_2 = cbind(t_i, case_0[[1]][2,], dler_vals[1] + mean(dler_val))
        prob1_1_to_2 = cbind(t_i, case_1[[1]][2,], dler_vals[1] + mean(dler_val))
        prob0_1_to_3 = cbind(t_i, case_0[[1]][3,], dler_vals[1] + mean(dler_val))
        prob1_1_to_3 = cbind(t_i, case_1[[1]][3,], dler_vals[1] + mean(dler_val))
        for(c in 2:length(dler_vals)) {
            temp1 = cbind(t_i, case_0[[c]][1,], dler_vals[c] + mean(dler_val))
            prob0_1_to_1 = rbind(prob0_1_to_1, temp1)
            
            temp1 = cbind(t_i, case_1[[c]][1,], dler_vals[c] + mean(dler_val))
            prob1_1_to_1 = rbind(prob1_1_to_1, temp1)
            
            temp1 = cbind(t_i, case_0[[c]][2,], dler_vals[c] + mean(dler_val))
            prob0_1_to_2 = rbind(prob0_1_to_2, temp1)
            
            temp1 = cbind(t_i, case_1[[c]][2,], dler_vals[c] + mean(dler_val))
            prob1_1_to_2 = rbind(prob1_1_to_2, temp1)
            
            temp1 = cbind(t_i, case_0[[c]][3,], dler_vals[c] + mean(dler_val))
            prob0_1_to_3 = rbind(prob0_1_to_3, temp1)
            
            temp1 = cbind(t_i, case_1[[c]][3,], dler_vals[c] + mean(dler_val))
            prob1_1_to_3 = rbind(prob1_1_to_3, temp1)
        }
        
        colnames(prob0_1_to_1) = colnames(prob1_1_to_1) = colnames(prob0_1_to_2) = 
            colnames(prob1_1_to_2) = colnames(prob0_1_to_3) = colnames(prob1_1_to_3) = 
            c('t', 'prob', 'dler')
        
        prob0_1_to_1 = as.data.frame(prob0_1_to_1)
        prob1_1_to_1 = as.data.frame(prob1_1_to_1)
        prob0_1_to_2 = as.data.frame(prob0_1_to_2)
        prob1_1_to_2 = as.data.frame(prob1_1_to_2)
        prob0_1_to_3 = as.data.frame(prob0_1_to_3)
        prob1_1_to_3 = as.data.frame(prob1_1_to_3)
        
        # 1 -> 1 transition ----------------------------------------------------
        if(covariate_struct == 3) {
            sub_title_0 = TeX(r'((sex $=$ 0) )')
        } else {
            sub_title_0 = " "
        }
        plot0_1_to_1 = ggplot(data = prob0_1_to_1, aes(x = t, y = prob, 
                                                       color = as.integer(dler),
                                                       group = dler)) +
            geom_line() +
            ylim(c(0,1)) + 
            scale_colour_gradient(name = "DLER value", 
                                  low = "blue", high = "red") +
            labs(x = TeX(r'(Time post baseline (x30 sec) )'), 
                 y = TeX(r'(Probability of baseline state)'), 
                 title = TeX(r'(Probability evolution curve from baseline state)'),
                 subtitle = sub_title_0) +
            theme(legend.position = "bottom",
                  legend.key.height = unit(0.25, 'cm'), #change legend key height
                  legend.key.width = unit(1, 'cm'), #change legend key width
                  legend.title = element_text(size=10), #change legend title font size
                  legend.text = element_text(size=7)) +
            guides(colour = guide_colorbar(title.position = "left",title.vjust = 1))
        
        ggsave(filename = paste0("Plots/probEvo_trial", trial_num, "_1to1_0.png"), 
               plot = plot0_1_to_1, width = 1500, height = 1000, units = "px")

        if(covariate_struct == 3) {
            plot1_1_to_1 = ggplot(data = prob1_1_to_1, aes(x = t, y = prob, 
                                                           color = as.integer(dler), 
                                                           group = dler)) +
                geom_line() +
                ylim(c(0,1)) + 
                scale_colour_gradient(name = "DLER value", 
                                      low = "blue", high = "red") +
                labs(x = TeX(r'(Time post baseline (x30 sec) )'), 
                     y = TeX(r'(Probability of baseline state)'), 
                     title = TeX(r'(Probability evolution curve from baseline state)'),
                     subtitle = TeX(r'((sex $=$ 1) )')) +
                theme(legend.position = "bottom",
                      legend.key.height = unit(0.25, 'cm'), #change legend key height
                      legend.key.width = unit(1, 'cm'), #change legend key width
                      legend.title = element_text(size=10), #change legend title font size
                      legend.text = element_text(size=7)) +
                guides(colour = guide_colorbar(title.position = "left",title.vjust = 1))   
            
            ggsave(filename = paste0("Plots/probEvo_trial", trial_num, "_1to1_1.png"), 
                   plot = plot1_1_to_1, width = 1500, height = 1000, units = "px")
        }
        
        # 1 -> 2 transition ----------------------------------------------------
        plot0_1_to_2 = ggplot(data = prob0_1_to_2, aes(x = t, y = prob, 
                                                       color = as.integer(dler),
                                                       group = dler)) +
            geom_line() +
            ylim(c(0,1)) + 
            scale_colour_gradient(name = "DLER value", 
                                  low = "blue", high = "red") +
            labs(x = TeX(r'(Time post baseline (x30 sec) )'), 
                 y = TeX(r'(Transition probability)'), 
                 title = TeX(r'(Probability evolution curve (state 1 $\to$ 2) )'),
                 subtitle = sub_title_0) +
            theme(legend.position = "bottom",
                  legend.key.height = unit(0.25, 'cm'), #change legend key height
                  legend.key.width = unit(1, 'cm'), #change legend key width
                  legend.title = element_text(size=10), #change legend title font size
                  legend.text = element_text(size=7)) +
            guides(colour = guide_colorbar(title.position = "left",title.vjust = 1))
        
        # ggsave(filename = paste0("Plots/probEvo_trial", trial_num, "_1to2_0.png"), 
        #        plot = plot0_1_to_2, width = 1500, height = 1000, units = "px")
        
        if(covariate_struct == 3) {
            plot1_1_to_2 = ggplot(data = prob1_1_to_2, aes(x = t, y = prob, 
                                                           color = as.integer(dler), 
                                                           group = dler)) +
                geom_line() +
                ylim(c(0,1)) + 
                scale_colour_gradient(name = "DLER value", 
                                      low = "blue", high = "red") +
                labs(x = TeX(r'(Time post baseline (x30 sec) )'), 
                     y = TeX(r'(Transition probability)'), 
                     title = TeX(r'(Probability evolution curve (state 1 $\to$ 2) )'),
                     subtitle = TeX(r'((sex $=$ 1) )')) +
                theme(legend.position = "bottom",
                      legend.key.height = unit(0.25, 'cm'), #change legend key height
                      legend.key.width = unit(1, 'cm'), #change legend key width
                      legend.title = element_text(size=10), #change legend title font size
                      legend.text = element_text(size=7)) +
                guides(colour = guide_colorbar(title.position = "left",title.vjust = 1))
            
            # ggsave(filename = paste0("Plots/probEvo_trial", trial_num, "_1to2_1.png"), 
            #        plot = plot1_1_to_2, width = 1500, height = 1000, units = "px")
        }
        
        # 1 -> 3 transition ----------------------------------------------------
        plot0_1_to_3 = ggplot(data = prob0_1_to_3, aes(x = t, y = prob, 
                                                       color = as.integer(dler),
                                                       group = dler)) +
            geom_line() +
            ylim(c(0,1)) + 
            scale_colour_gradient(name = "DLER value", 
                                  low = "blue", high = "red") +
            labs(x = TeX(r'(Time post baseline (x30 sec) )'), 
                 y = TeX(r'(Transition probability)'), 
                 title = TeX(r'(Probability evolution curve (state 1 $\to$ 3) )'),
                 subtitle = sub_title_0) +
            theme(legend.position = "bottom",
                  legend.key.height = unit(0.25, 'cm'), #change legend key height
                  legend.key.width = unit(1, 'cm'), #change legend key width
                  legend.title = element_text(size=10), #change legend title font size
                  legend.text = element_text(size=7)) +
            guides(colour = guide_colorbar(title.position = "left",title.vjust = 1))
        
        # ggsave(filename = paste0("Plots/probEvo_trial", trial_num, "_1to3_0.png"), 
        #        plot = plot0_1_to_3, width = 1500, height = 1000, units = "px")
        
        if(covariate_struct == 3) {
            plot1_1_to_3 = ggplot(data = prob1_1_to_3, aes(x = t, y = prob, 
                                                           color = as.integer(dler), 
                                                           group = dler)) +
                geom_line() +
                ylim(c(0,1)) + 
                scale_colour_gradient(name = "DLER value", 
                                      low = "blue", high = "red") +
                labs(x = TeX(r'(Time post baseline (x30 sec) )'), 
                     y = TeX(r'(Transition probability)'), 
                     title = TeX(r'(Probability evolution curve (state 1 $\to$ 3) )'),
                     subtitle = TeX(r'((sex $=$ 1) )')) +
                theme(legend.position = "bottom",
                      legend.key.height = unit(0.25, 'cm'), #change legend key height
                      legend.key.width = unit(1, 'cm'), #change legend key width
                      legend.title = element_text(size=10), #change legend title font size
                      legend.text = element_text(size=7)) +
                guides(colour = guide_colorbar(title.position = "left",title.vjust = 1))   
            
            # ggsave(filename = paste0("Plots/probEvo_trial", trial_num, "_1to3_1.png"), 
            #        plot = plot1_1_to_3, width = 1500, height = 1000, units = "px")
        }
        
        app0 = ggarrange(plot0_1_to_2 + 
                             labs(title = TeX(r'(Probability evolution: baseline to state 2)')) +
                             theme(axis.title.x = element_blank(),
                                   plot.title = element_text(size = 8, face = "bold"),
                                   plot.subtitle = element_text(size = 7)), 
                         plot0_1_to_3 + 
                             labs(title = TeX(r'(Probability evolution: baseline to state 3)')) +
                             theme(axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   axis.title.y = element_blank(),
                                   axis.title.x = element_blank(),
                                   legend.position="none",
                                   plot.title = element_text(size = 8, face = "bold"),
                                   plot.subtitle = element_text(size = 7)),
                         nrow = 1, ncol = 2)
        
        ggsave(filename = paste0("Plots/probEvo_trial", trial_num, "_combo_0.png"), 
               plot = app0, width = 1500, height = 1000, units = "px")
        
        if(covariate_struct == 3) {
            app1 = ggarrange(plot1_1_to_2 + 
                                 labs(title = TeX(r'(Probability evolution: baseline to state 2)'))+
                                 theme(axis.title.x = element_blank(),
                                       plot.title = element_text(size = 8, face = "bold"),
                                       plot.subtitle = element_text(size = 7)), 
                             plot1_1_to_3 + 
                                 labs(title = TeX(r'(Probability evolution: baseline to state 3)')) +
                                 theme(axis.text.y = element_blank(),
                                       axis.ticks.y = element_blank(),
                                       axis.title.y = element_blank(),
                                       axis.title.x = element_blank(),
                                       legend.position="none",
                                       plot.title = element_text(size = 8, face = "bold"),
                                       plot.subtitle = element_text(size = 7)),
                             nrow = 1, ncol = 2)   
            
            ggsave(filename = paste0("Plots/probEvo_trial", trial_num, "_combo_1.png"), 
                   plot = app1, width = 1500, height = 1000, units = "px")
        }
          
    } else {
        # Sex = 0
        prob0_1_to_1 = data.frame('t' = t_i, 'prob' = case_0[[1]][1,])
        plot0_1_to_1 = ggplot(data = prob0_1_to_1, aes(x = t, y = prob)) +
            ylim(c(0,1)) + 
            geom_line(linewidth=1) +
            labs(x = TeX(r'(Time post baseline (x30 sec) )'), 
                 y = TeX(r'(Probability of baseline state)'), 
                 title = TeX(r'(Probability evolution curve from baseline state)'),
                 subtitle = TeX(r'((sex $=$ 0) )'))
        
        prob0_1_to_2 = data.frame('t' = t_i, 'prob' = case_0[[1]][2,])
        plot0_1_to_2 = ggplot(data = prob0_1_to_2, aes(x = t, y = prob)) +
            ylim(c(0,1)) + 
            geom_line(linewidth=1) +
            labs(x = TeX(r'(Time post baseline (x30 sec) )'), 
                 y = TeX(r'(Transition probability)'), 
                 title = TeX(r'(Probability evolution curve (state 1 $\to$ 2) )'),
                 subtitle = TeX(r'((sex $=$ 0) )'))
        
        prob0_1_to_3 = data.frame('t' = t_i, 'prob' = case_0[[1]][3,])
        plot0_1_to_3 = ggplot(data = prob0_1_to_3, aes(x = t, y = prob)) +
            ylim(c(0,1)) + 
            geom_line(linewidth=1) +
            labs(x = TeX(r'(Time post baseline (x30 sec) )'), 
                 y = TeX(r'(Transition probability)'), 
                 title = TeX(r'(Probability evolution curve (state 1 $\to$ 3) )'),
                 subtitle = TeX(r'((sex $=$ 0) )'))
        
        ggsave(filename = paste0("Plots/probEvo_trial", trial_num, "_1to1_0.png"), 
               plot = plot0_1_to_1, width = 1500, height = 1000, units = "px")
        # ggsave(filename = paste0("Plots/probEvo_trial", trial_num, "_1to2_0.png"), 
        #        plot = plot0_1_to_2, width = 1500, height = 1000, units = "px")
        # ggsave(filename = paste0("Plots/probEvo_trial", trial_num, "_1to3_0.png"), 
        #        plot = plot0_1_to_3, width = 1500, height = 1000, units = "px")
        
        app0 = ggarrange(plot0_1_to_2 + 
                             labs(title = TeX(r'(Probability evolution: baseline to state 2)')) +
                             theme(axis.title.x = element_blank(),
                                   plot.title = element_text(size = 8, face = "bold"),
                                   plot.subtitle = element_text(size = 7)), 
                         plot0_1_to_3 + 
                             labs(title = TeX(r'(Probability evolution: baseline to state 3)')) +
                             theme(axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   axis.title.y = element_blank(),
                                   axis.title.x = element_blank(),
                                   legend.position="none",
                                   plot.title = element_text(size = 8, face = "bold"),
                                   plot.subtitle = element_text(size = 7)),
                         nrow = 1, ncol = 2)
        ggsave(filename = paste0("Plots/probEvo_trial", trial_num, "_combo_0.png"), 
               plot = app0, width = 1500, height = 1000, units = "px")


        # Sex = 1
        prob0_1_to_1 = data.frame('t' = t_i, 'prob' = case_1[[1]][1,])
        plot0_1_to_1 = ggplot(data = prob0_1_to_1, aes(x = t, y = prob)) +
            ylim(c(0,1)) + 
            geom_line(linewidth=1) +
            labs(x = TeX(r'(Time post baseline (x30 sec) )'), 
                 y = TeX(r'(Probability of baseline state)'), 
                 title = TeX(r'(Probability evolution curve from baseline state)'),
                 subtitle = TeX(r'((sex $=$ 1) )'))
        
        prob0_1_to_2 = data.frame('t' = t_i, 'prob' = case_1[[1]][2,])
        plot0_1_to_2 = ggplot(data = prob0_1_to_2, aes(x = t, y = prob)) +
            ylim(c(0,1)) + 
            geom_line(linewidth=1) +
            labs(x = TeX(r'(Time post baseline (x30 sec) )'), 
                 y = TeX(r'(Transition probability)'), 
                 title = TeX(r'(Probability evolution curve (state 1 $\to$ 2) )'),
                 subtitle = TeX(r'((sex $=$ 1) )'))
        
        prob0_1_to_3 = data.frame('t' = t_i, 'prob' = case_1[[1]][3,])
        plot0_1_to_3 = ggplot(data = prob0_1_to_3, aes(x = t, y = prob)) +
            ylim(c(0,1)) + 
            geom_line(linewidth=1) +
            labs(x = TeX(r'(Time post baseline (x30 sec) )'), 
                 y = TeX(r'(Transition probability)'), 
                 title = TeX(r'(Probability evolution curve (state 1 $\to$ 3) )'),
                 subtitle = TeX(r'((sex $=$ 1) )'))
        
        ggsave(filename = paste0("Plots/probEvo_trial", trial_num, "_1to1_1.png"), 
               plot = plot0_1_to_1, width = 1500, height = 1000, units = "px")
        # ggsave(filename = paste0("Plots/probEvo_trial", trial_num, "_1to2_1.png"), 
        #        plot = plot0_1_to_2, width = 1500, height = 1000, units = "px")
        # ggsave(filename = paste0("Plots/probEvo_trial", trial_num, "_1to3_1.png"), 
        #        plot = plot0_1_to_3, width = 1500, height = 1000, units = "px")
        
        app0 = ggarrange(plot0_1_to_2 + 
                             labs(title = TeX(r'(Probability evolution: baseline to state 2)')) +
                             theme(axis.title.x = element_blank(),
                                   plot.title = element_text(size = 8, face = "bold"),
                                   plot.subtitle = element_text(size = 7)), 
                         plot0_1_to_3 + 
                             labs(title = TeX(r'(Probability evolution: baseline to state 3)')) +
                             theme(axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   axis.title.y = element_blank(),
                                   axis.title.x = element_blank(),
                                   legend.position="none",
                                   plot.title = element_text(size = 8, face = "bold"),
                                   plot.subtitle = element_text(size = 7)),
                         nrow = 1, ncol = 2)
        ggsave(filename = paste0("Plots/probEvo_trial", trial_num, "_combo_1.png"), 
               plot = app0, width = 1500, height = 1000, units = "px")
    }
    
}