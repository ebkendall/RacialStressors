# Rcpp packages
library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)
sourceCpp("mcmc_routine_c.cpp")

# Needed for OpenMP C++ parallel
Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
Sys.setenv("PKG_LIBS" = "-fopenmp")

# Initialization --------------------------------------------------------------
set.seed(2023)
dir = 'Model_out/'

# Information defining which approach to take ----------------------------------
trial_num = 9
simulation = F
case_b = T

index_seeds = c(1:5)
# ------------------------------------------------------------------------------


# Load the posterior samples of the HMM parameters ----------------------------
n_post = 495000; burnin = 5000; steps = 500000
index_post = (steps - burnin - n_post + 1):(steps - burnin)

par_chain = NULL

for (seed in index_seeds) {
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
    
    load(file_name)
    main_chain = mcmc_out$chain[index_post,]
    ind_keep = seq(1, nrow(main_chain), by=100)
    
    par_chain_i = main_chain[ind_keep, ]
    par_chain   = rbind(par_chain, par_chain_i)
}
# -----------------------------------------------------------------------------


# State space sampler ---------------------------------------------------------
if(simulation) {
    # Simulation
    load('Data/sim_data_1_30.rda')
    data_format = sim_data
} else {
    # Real data analysis
    load('Data/data_format_30.rda')
    data_format = data_format_30
    
    miss_info = c(26296, 29698, 30625, 401, 423, 419, 457)
    data_format = data_format[!(data_format[,"ID.."] %in% miss_info), ]
}

n_sub = length(unique(data_format[,'ID..']))

# par_index = list(zeta=1:30, misclass=42:45, delta = 31:33, tau2 = 34, sigma2 = 35:37,
#                  gamma = 38:41)
par_index = list(zeta=1:6, misclass=15:18, delta = 7:9, tau2 = 10, sigma2 = 11:13)

temp_data = as.matrix(data_format); rownames(temp_data) = NULL
id = as.numeric(temp_data[,"ID.."])
y_1 = as.numeric(temp_data[,"State"])
y_2 = as.numeric(temp_data[,"RSA"])
t = as.numeric(temp_data[,"Time"])
EIDs = unique(id)

if(simulation & case_b) {
    y_1 = as.numeric(temp_data[,"Alt_state"])
}

# cov_info = temp_data[,c("Age", "sex1", "edu_yes", "DLER_avg"), drop=F]
cov_info = matrix(0, nrow = nrow(temp_data), ncol = 1)

# Centering age
# if(!simulation) {
#     load(paste0('Data/mean_age_', trial_num, '.rda'))
#     load(paste0('Data/mean_dler_', trial_num, '.rda'))
#     cov_info[,'Age'] = cov_info[,'Age'] - mean_age
#     cov_info[,'DLER_avg'] = cov_info[,'DLER_avg'] - mean_dler   
# }

prior_mean = rep(0, ncol(par_chain))
prior_sd = rep(20, ncol(par_chain))

prior_par = list()
prior_par[[1]] = prior_mean
prior_par[[2]] = prior_sd

new_steps =  100000
new_burnin = 5000

B_chain_obs = vector(mode = 'list', length = length(EIDs))
for(i in 1:length(EIDs)) B_chain_obs[[i]] = y_1[id == EIDs[i]]

B_chain_MLE = viterbi_alg(EIDs, apply(par_chain, 2, median), par_index, id, y_2, y_1, cov_info, case_b)

B_chain = state_space_sampler_no_label(new_steps, new_burnin, EIDs, par_chain, 
                                       par_index, y_2, id, t, y_1, cov_info, case_b)

file_name = NULL
if(simulation) {
    if(case_b) {
        file_name = paste0("Model_out/B_chain_", trial_num, "_sim_30b.rda")
        file_name2 = paste0("Model_out/B_chain_", trial_num, "_sim_30b_MLE.rda")
    } else {
        file_name = paste0("Model_out/B_chain_", trial_num, "_sim_30.rda")   
        file_name2 = paste0("Model_out/B_chain_", trial_num, "_sim_30_MLE.rda")  
    }
} else {
    if(case_b) {
        file_name = paste0("Model_out/B_chain_", trial_num, "_30b_s1_up.rda")
        file_name2 = paste0("Model_out/B_chain_", trial_num, "_30b_s1_MLE_up.rda")
    } else {
        file_name = paste0("Model_out/B_chain_", trial_num, "_30.rda")   
        file_name2 = paste0("Model_out/B_chain_", trial_num, "_30_MLE.rda")   
    }
}
save(B_chain, file = file_name)
save(B_chain_MLE, file = file_name2)
# -----------------------------------------------------------------------------