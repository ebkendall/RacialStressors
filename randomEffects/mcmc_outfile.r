# This script file produces trace plots and histograms of the mcmc output files
library(latex2exp)

dir = 'Model_out/' 

# Model type -------------------------------------------------------------------
# 1: baseline only
# 2: baseline & DLER
# 3: all covariates

covariate_struct = 3
# ------------------------------------------------------------------------------

# Information defining which approach to take ----------------------------------
trial_num = 10
simulation = F
case_b = T
# ------------------------------------------------------------------------------

# Size of posterior sample from mcmc chains
n_post = 195000
# Step number at which the adaptive tuning scheme was frozen
burnin = 5000
# Total number of steps the mcmc algorithm is computed for
steps = 200000
# Matrix row indices for the posterior sample to use
index_post = (steps - burnin - n_post + 1):(steps - burnin)

index_seeds = c(1:5)

if(covariate_struct == 1) {
    par_index = list(zeta=1:6, misclass=14:17, delta = 7:9, tau2 = 10, 
                     sigma2 = 11:13, gamma = -1)
    labels <- c(TeX(r'($\hat{\zeta}_{0,1}:$ baseline: 1 $\to$ 2)'), 
                TeX(r'($\hat{\zeta}_{0,2}:$ baseline: 1 $\to$ 3)'),
                TeX(r'($\hat{\zeta}_{0,3}:$ baseline: 2 $\to$ 1)'),
                TeX(r'($\hat{\zeta}_{0,4}:$ baseline: 2 $\to$ 3)'),
                TeX(r'($\hat{\zeta}_{0,5}:$ baseline: 3 $\to$ 1)'),
                TeX(r'($\hat{\zeta}_{0,6}:$ baseline: 3 $\to$ 2)'),
                TeX(r'($\delta_1 = \mu$)'), TeX(r'($\delta_2 = \alpha$)'), TeX(r'($\delta_3 = \beta$)'),
                TeX(r'($\log(\tau^2)$)'), TeX(r'($\log(\sigma_1^2)$)'), TeX(r'($\log(\sigma_2^2)$)'), TeX(r'($\log(\sigma_3^2)$)'),
                TeX(r'(logit P(obs S2 | true S1))'), TeX(r'(logit P(obs S3 | true S1))'),
                TeX(r'(logit P(obs S3 | true S2))'), TeX(r'(logit P(obs S2 | true S3))'),
                TeX(r'($\mu + \alpha$)'), TeX(r'($\mu + \beta$)'),
                TeX(r'($\tau^2 + \sigma_1^2$)'), TeX(r'($\tau^2 + \sigma_2^2$)'), 
                TeX(r'($\tau^2 + \sigma_3^2$)'), TeX(r'($\tau^2$)'),
                TeX(r'($\sigma_1^2$)'), TeX(r'($\sigma_2^2$)'), 
                TeX(r'($\sigma_3^2$)'))
} else if(covariate_struct == 2) {
    par_index = list(zeta=1:12, misclass=21:24, delta = 13:15, tau2 = 16, 
                     sigma2 = 17:19, gamma = 20)

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
                TeX(r'($\delta_1 = \mu$)'), TeX(r'($\delta_2 = \alpha$)'), TeX(r'($\delta_3 = \beta$)'),
                TeX(r'($\log(\tau^2)$)'), TeX(r'($\log(\sigma_1^2)$)'), TeX(r'($\log(\sigma_2^2)$)'), TeX(r'($\log(\sigma_3^2)$)'),
                TeX(r'($\hat{\gamma}_1:$ DLER)'),
                TeX(r'(logit P(obs S2 | true S1))'), TeX(r'(logit P(obs S3 | true S1))'),
                TeX(r'(logit P(obs S3 | true S2))'), TeX(r'(logit P(obs S2 | true S3))'),
                TeX(r'($\mu + \alpha$)'), TeX(r'($\mu + \beta$)'),
                TeX(r'($\tau^2 + \sigma_1^2$)'), TeX(r'($\tau^2 + \sigma_2^2$)'), 
                TeX(r'($\tau^2 + \sigma_3^2$)'), TeX(r'($\tau^2$)'),
                TeX(r'($\sigma_1^2$)'), TeX(r'($\sigma_2^2$)'), 
                TeX(r'($\sigma_3^2$)'))
} else {
    par_index = list(zeta=1:30, misclass=42:45, delta = 31:33, tau2 = 34, sigma2 = 35:37,
                 gamma = 38:41)

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
                TeX(r'($\delta_1 = \mu$)'), TeX(r'($\delta_2 = \alpha$)'), TeX(r'($\delta_3 = \beta$)'),
                TeX(r'($\log(\tau^2)$)'), TeX(r'($\log(\sigma_1^2)$)'), TeX(r'($\log(\sigma_2^2)$)'), TeX(r'($\log(\sigma_3^2)$)'),
                TeX(r'($\hat{\gamma}_1:$ age)'), TeX(r'($\hat{\gamma}_2:$ sex1)'), 
                TeX(r'($\hat{\gamma}_3:$ yes edu)'), TeX(r'($\hat{\gamma}_4:$ DLER)'),
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

# load(paste0('Data/mean_age_', trial_num, '.rda'))

ind = 0

for(seed in index_seeds){
    
    file_name = NULL
    
    if(simulation) {
        if(case_b) {
            file_name = paste0(dir,'mcmc_out_',toString(seed), '_', trial_num, '_sim_30b.rda')   
        } else {
            file_name = paste0(dir,'mcmc_out_',toString(seed), '_', trial_num, '_sim_30.rda')      
        }
    } else {
        if(case_b) {
            file_name = paste0(dir,'mcmc_out_',toString(seed), '_', trial_num, '_30b.rda')   
        } else {
            file_name = paste0(dir,'mcmc_out_',toString(seed), '_', trial_num, '_30.rda')      
        }
    }
    
    if (file.exists(file_name)) {
        load(file_name)
        ind = ind + 1

        print(mcmc_out$accept)

        # Thinning the chain
        main_chain = mcmc_out$chain[index_post,]
        ind_keep = seq(1, nrow(main_chain), by=10)
        
        # # Un-doing the centering
        # main_chain[,par_index$zeta[1:6]] = main_chain[,par_index$zeta[1:6]] +
        #     mean_age * main_chain[,par_index$zeta[7:12]]
        # main_chain[,par_index$delta[1]] = main_chain[,par_index$delta[1]] +
        #     mean_age * main_chain[,par_index$gamma[1]]
        
        
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
    }
}

print(length(labels))
print(ncol(main_chain))

# Plot and save the mcmc trace plots and histograms.
pdf_title = NULL
if(simulation) {
    if(case_b) {
        pdf_title = paste0('Plots/mcmc_out_', trial_num, '_sim_30b.pdf')
    } else {
        pdf_title = paste0('Plots/mcmc_out_', trial_num, '_sim_30.pdf')   
    }
} else {
    if(case_b) {
        pdf_title = paste0('Plots/mcmc_out_', trial_num, '_30b.pdf')
    } else {
        pdf_title = paste0('Plots/mcmc_out_', trial_num, '_30.pdf')   
    }
}
print(pdf_title)

pdf(pdf_title)
par(mfrow=c(4, 2))

stacked_chains = do.call( rbind, chain_list)
par_mean = par_median = upper = lower = rep( NA, length(labels))

if(simulation) {
    load('Data/true_par_7_30.rda')
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
    
    if(!simulation) {
        if(r == par_index$delta[1]) abline( v=6.46408805031447, col='blue', lwd=2, lty=2)
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

print(par_median)
save(par_median, file = paste0('Model_out/par_median', trial_num, '.rda'))

dev.off()

# if(!simulation) {
#     # Probability of transitioning in 30s with certain covariate combinations -----
#     load('../Data/data_format_30.rda')
#     data_format = data_format_30
    
#     # Removing the participants with missing labels
#     miss_info = c(26296, 29698, 30625, 401, 423, 419, 457)
#     data_format = data_format[!(data_format[,"ID.."] %in% miss_info), ]
    
#     ages = NULL
#     dler_val = NULL
#     for(a in unique(data_format[,"ID.."])) {
#         ages = c(ages, unique(data_format[data_format[,"ID.."] == a, "Age"]))
#         dler_val = c(dler_val, unique(data_format[data_format[,"ID.."] == a, "DLER_avg"]))
#     }
#     s_dler = sd(dler_val)
#     print(paste0("Mean DLER: ", mean(dler_val)))
#     print(paste0("Standard dev. DLER: ", s_dler))
    
#     sex = c(0,1)
#     pEdu = c(0,1)
#     dler = c(-s_dler,0,s_dler)
#     zeta_est = matrix(par_median[par_index$zeta], nrow=6, ncol=5)
    
#     for(s in sex) {
#         for(p in pEdu) {
#             for(d in dler) {
#                 print(paste0("sex: ", s, ", pEdu: ", p, ", DLER: ", d))
#                 cov_val = matrix(c(1, 0, s, p, d), ncol = 1)
#                 qs = zeta_est %*% cov_val
#                 Q = matrix(c(         1, exp(qs[1]), exp(qs[2]),
#                                       exp(qs[3]),          1, exp(qs[4]),
#                                       exp(qs[5]), exp(qs[6]),          1), ncol = 3, byrow=T)
#                 P = Q / rowSums(Q)
#                 print(paste0('state 1 -> 2: ', round(P[1,2], digits=4)))
#                 print(paste0('state 1 -> 3: ', round(P[1,3], digits=4)))
#                 print(paste0('sum of probs: ', round(P[1,2] + P[1,3], digits=4)))
#             }
#         }
#     }
    
#     # Probability evolution curves ------------------------------------------------
#     prob_evo <- function( par, par_index, t_i, cov_val){
        
#         zeta_est = matrix(par[par_index$zeta], nrow=6, ncol=5)
        
#         # Compute prob evolutions of transitioning 1->2, 1->3, 2->3
#         prob_evo_mat = matrix( 0, 3, length(t_i))
#         P_t = matrix(c(1,0,0), ncol = 3)
#         prob_evo_mat[,1] = c(P_t)
        
#         for(k in 2:length(t_i)){
            
#             qs = zeta_est %*% cov_val
#             Q = matrix(c(         1, exp(qs[1]), exp(qs[2]),
#                                   exp(qs[3]),          1, exp(qs[4]),
#                                   exp(qs[5]), exp(qs[6]),          1), ncol = 3, byrow=T)
#             P = Q / rowSums(Q)
            
#             P_t = P_t %*% P
#             prob_evo_mat[,k] = c(P_t)
#         }
        
#         return(prob_evo_mat)
#     }
    
#     # DLER range: -1.31355932 -> 2.56879362
#     # Sex  = 0, 1
    
#     # baseline, age, sex, pEdu, dler
#     dler_vals = seq(2.6, -1.4, by = -0.2)
#     t_i = 1:200
    
#     colfunc <- colorRampPalette(c("red", "green"))
#     col_dler = colfunc(length(dler_vals))
    
#     # Sex = 0
#     case_0 = vector(mode = 'list', length = length(dler_vals))
#     # Sex = 1
#     case_1 = vector(mode = 'list', length = length(dler_vals))
#     for(i in 1:length(dler_vals)) {
#         cov_val_0 = matrix(c(1, 0, 0, 0, dler_vals[i]), ncol = 1)
#         cov_val_1 = matrix(c(1, 0, 1, 0, dler_vals[i]), ncol = 1)
        
#         case_0[[i]] = prob_evo(par_median, par_index, t_i, cov_val_0)
#         case_1[[i]] = prob_evo(par_median, par_index, t_i, cov_val_1)
#     }
    
#     png(filename = "Plots/probEvo1to2_0.png", width = 1500, height = 1000)
#     par(mar=c(10,10,10,10))
#     plot(x = t_i, y = case_0[[1]][2,], type = 'l', col = col_dler[1], lwd = 2,
#          xlab = TeX(r'(Time post baseline (x30 sec) )'), 
#          ylab = TeX(r'(Transition probability)'), ylim = c(0,1),
#          cex.lab=4, cex.axis=3, mgp=c(6,2,0))
#     mytitle = TeX(r'(Probability evolution curve (state 1 $\to$ 2) )')
#     mysubtitle = TeX(r'((sex $=$ 0) )')
#     mtext(side=3, line=3.5, at=-0.07, adj=0, cex=5, mytitle)
#     mtext(side=3, line=0.5, at=-0.07, adj=0, cex=4, mysubtitle)
#     for(i in 2:length(dler_vals)) {
#         lines(x = t_i, y = case_0[[i]][2,], col = col_dler[i], lwd = 2)
#     }
#     lines(x = t_i, y = case_0[[which(dler_vals == 0)]][2,], col = 'blue', lwd = 4)
#     legend( 'topleft', horiz=F, bty='n', x.intersp=.75,
#             legend=c( 'DLER = +2.6', 'DLER =   0.0', 'DLER =  -1.4'), pch=15, pt.cex=3, 
#             col=c(col_dler[1], 'blue', col_dler[length(col_dler)]), cex = 3)
#     dev.off()
    
#     png(filename = "Plots/probEvo1to2_1.png", width = 1500, height = 1000)
#     par(mar=c(10,10,10,10))
#     plot(x = t_i, y = case_1[[1]][2,], type = 'l', col = col_dler[1], lwd = 2,
#          xlab = TeX(r'(Time post baseline (x30 sec) )'), 
#          ylab = TeX(r'(Transition probability)'), ylim = c(0,1),
#          cex.lab=4, cex.axis=3, mgp=c(6,2,0))
#     mytitle = TeX(r'(Probability evolution curve (state 1 $\to$ 2) )')
#     mysubtitle = TeX(r'((sex $=$ 1) )')
#     mtext(side=3, line=3.5, at=-0.07, adj=0, cex=5, mytitle)
#     mtext(side=3, line=0.5, at=-0.07, adj=0, cex=4, mysubtitle)
#     for(i in 2:length(dler_vals)) {
#         lines(x = t_i, y = case_1[[i]][2,], col = col_dler[i], lwd = 2)
#     }
#     lines(x = t_i, y = case_1[[which(dler_vals == 0)]][2,], col = 'blue', lwd = 4)
#     legend( 'topleft', horiz=F, bty='n', x.intersp=.75,
#             legend=c( 'DLER = +2.6', 'DLER =   0.0', 'DLER =  -1.4'), pch=15, pt.cex=3, 
#             col=c(col_dler[1], 'blue', col_dler[length(col_dler)]), cex = 3)
#     dev.off()
    
#     png(filename = "Plots/probEvo1to3_0.png", width = 1500, height = 1000)
#     par(mar=c(10,10,10,10))
#     plot(x = t_i, y = case_0[[1]][3,], type = 'l', col = col_dler[1], lwd = 2,
#          xlab = TeX(r'(Time post baseline (x30 sec) )'), 
#          ylab = TeX(r'(Transition probability)'), ylim = c(0,1),
#          cex.lab=4, cex.axis=3, mgp=c(6,2,0))
#     mytitle = TeX(r'(Probability evolution curve (state 1 $\to$ 3) )')
#     mysubtitle = TeX(r'((sex $=$ 0) )')
#     mtext(side=3, line=3.5, at=-0.07, adj=0, cex=5, mytitle)
#     mtext(side=3, line=0.5, at=-0.07, adj=0, cex=4, mysubtitle)
#     for(i in 2:length(dler_vals)) {
#         lines(x = t_i, y = case_0[[i]][3,], col = col_dler[i], lwd = 2)
#     }
#     lines(x = t_i, y = case_0[[which(dler_vals == 0)]][3,], col = 'blue', lwd = 4)
#     legend( 'topleft', horiz=F, bty='n', x.intersp=.75,
#             legend=c( 'DLER = +2.6', 'DLER =   0.0', 'DLER =  -1.4'), pch=15, pt.cex=3, 
#             col=c(col_dler[1], 'blue', col_dler[length(col_dler)]), cex = 3)
#     dev.off()
    
#     png(filename = "Plots/probEvo1to3_1.png", width = 1500, height = 1000)
#     par(mar=c(10,10,10,10))
#     plot(x = t_i, y = case_1[[1]][3,], type = 'l', col = col_dler[1], lwd = 2,
#          xlab = TeX(r'(Time post baseline (x30 sec) )'), 
#          ylab = TeX(r'(Transition probability)'), ylim = c(0,1),
#          cex.lab=4, cex.axis=3, mgp=c(6,2,0))
#     mytitle = TeX(r'(Probability evolution curve (state 1 $\to$ 3) )')
#     mysubtitle = TeX(r'((sex $=$ 1) )')
#     mtext(side=3, line=3.5, at=-0.07, adj=0, cex=5, mytitle)
#     mtext(side=3, line=0.5, at=-0.07, adj=0, cex=4, mysubtitle)
#     for(i in 2:length(dler_vals)) {
#         lines(x = t_i, y = case_1[[i]][3,], col = col_dler[i], lwd = 2)
#     }
#     lines(x = t_i, y = case_1[[which(dler_vals == 0)]][3,], col = 'blue', lwd = 4)
#     legend( 'topleft', horiz=F, bty='n', x.intersp=.75,
#             legend=c( 'DLER = +2.6', 'DLER =   0.0', 'DLER =  -1.4'), pch=15, pt.cex=3, 
#             col=c(col_dler[1], 'blue', col_dler[length(col_dler)]), cex = 3)
#     dev.off()
    
# }


