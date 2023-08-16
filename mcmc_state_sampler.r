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
index_seeds = c(1,3:5)

# Information defining which approach to take ----------------------------------
trial_num = 1
simulation = F
thirty = T
use_labels = T
# ------------------------------------------------------------------------------


# Load the posterior samples of the HMM parameters ----------------------------
n_post = 5000; burnin = 5000; steps = 30000
index_post = (steps - burnin - n_post + 1):(steps - burnin)

par_chain = NULL

for (seed in index_seeds) {
    file_name = NULL
    if(simulation) {
        if(thirty) {
            file_name = paste0(dir,'mcmc_out_',toString(seed), '_', trial_num, '_sim_30.rda')   
        } else {
            file_name = paste0(dir,'mcmc_out_',toString(seed), '_', trial_num, '_sim_15.rda')
        }
    } else {
        if(thirty) {
            file_name = paste0(dir,'mcmc_out_',toString(seed), '_', trial_num, '_30.rda')   
        } else {
            file_name = paste0(dir,'mcmc_out_',toString(seed), '_', trial_num, '_15.rda')   
        }
    }
    
    load(file_name)
    
    main_chain = mcmc_out$chain[index_post,]
    ind_keep = seq(1, nrow(main_chain), by=10)
    
    par_chain_i = main_chain[ind_keep, ]
    par_chain   = rbind(par_chain, par_chain_i)
}
# -----------------------------------------------------------------------------


# State space sampler ---------------------------------------------------------
if(simulation) {
    # Simulation
    load('Data/sim_data_1_a.rda')
    load('Data/true_par_a.rda')
    data_format = sim_data
} else {
    # Real data analysis
    if(thirty) {
        load('Data/data_format_30.rda')
        data_format = data_format_30
    } else {
        load('Data/data_format_15.rda')
        data_format = data_format_15   
    } 
    
    miss_info = c(26296, 29698, 30625, 401, 423, 419, 457)
    data_format = data_format[!(data_format[,"ID.."] %in% miss_info), ]
}

n_sub = length(unique(data_format[,'ID..']))

par_index = list( zeta=1:5, misclass=6:9,
                  delta = 10:12, tau2 = 13, sigma2 = 14)

temp_data = as.matrix(data_format); rownames(temp_data) = NULL
id = as.numeric(temp_data[,"ID.."])
y_1 = as.numeric(temp_data[,"State"])
y_2 = as.numeric(temp_data[,"RSA"])
t = as.numeric(temp_data[,"Time"])
EIDs = unique(id)

new_steps =  50000
new_burnin = 10000

B_chain = state_space_sampler(new_steps, new_burnin, EIDs, colMeans(par_chain), 
                              par_index, y_1, y_2, id, t)

file_name = NULL
if(simulation) {
    file_name = paste0("Model_out/B_chain_", trial_num, "_sim.rda")
} else {
    if(thirty) {
        file_name = paste0("Model_out/B_chain_", trial_num, "_30.rda")
    } else {
        file_name = paste0("Model_out/B_chain_", trial_num, "_15.rda")
    }
}
save(B_chain, file = file_name)
# -----------------------------------------------------------------------------