# Rcpp packages
library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)
sourceCpp("mcmc_routine_c.cpp")

set.seed(2023)

dir = 'Model_out/' 

# covariate_struct value indicators --------------------------------------------
# 1: baseline model (age, sex, pEdu)
# 2: DLER only
# 3: all covariates

args = commandArgs(TRUE)
covariate_struct = as.numeric(args[1])
# ------------------------------------------------------------------------------

trial_num = covariate_struct
simulation = T

# ------------------------------------------------------------------------------
n_post = 20000; burnin = 500; steps = 500000

# Matrix row indices for the posterior sample to use
index_post = (steps - burnin - n_post + 1):(steps - burnin)
# ------------------------------------------------------------------------------

init_par = NULL

if(simulation) {
    # Simulation
    # 30s epochs
    load(paste0('Data/sim_data_', covariate_struct, '_30.rda'))
    load(paste0('Data/true_par_', covariate_struct, '_30.rda'))
    data_format = sim_data

    index_seeds = 1
} else {
    # ********************************************
    # **** For real data contact first author ****
    # ********************************************
    index_seeds = c(1:5)
}

EIDs = unique(data_format[,"ID.."])

# Making the data into a matrix
temp_data = as.matrix(data_format); rownames(temp_data) = NULL
id = as.numeric(temp_data[,"ID.."])
y_1 = as.numeric(temp_data[,"State"])
y_2 = as.numeric(temp_data[,"RSA"])
t = as.numeric(temp_data[,"Time"])

if(simulation) {
    y_1 = as.numeric(temp_data[,"Alt_state"])
} 

# Defining the parameter and covariance structure
if(covariate_struct == 1) {
    # Baseline model w/ no DLER
    init_par = c(c(matrix(c(-2, 0, 0, 0,
                            -2, 0, 0, 0,
                            -2, 0, 0, 0,
                            -2, 0, 0, 0,
                            -2, 0, 0, 0,
                            -2, 0, 0, 0), ncol=4, byrow = T)),
             c(6.411967, 0, 0), 
             log(0.51^2),  
             c(log(0.8^2), 0, 0),
             0, 0, 0,
             0, 0, 0,
             -1, -1, -1, -1)
    par_index = list(zeta=1:24, misclass=38:41, delta = 25:27, tau2 = 28, 
                     sigma2 = 29:31, gamma = 32:34, delta_new = 35:37)
    
} else if(covariate_struct == 2) {
    # DLER only model
    init_par = c(c(matrix(c(-2, 0,
                            -2, 0,
                            -2, 0,
                            -2, 0,
                            -2, 0,
                            -2, 0), ncol=2, byrow = T)),
             c(6.411967, 0, 0), 
             log(0.51^2),  
             c(log(0.8^2), 0, 0),
             0, 0, 0,
             0, 0, 0,
             -1, -1, -1, -1)
    par_index = list(zeta=1:12, misclass=26:29, delta = 13:15, tau2 = 16, 
                     sigma2 = 17:19, gamma = 20:22, delta_new = 23:25)
} else {
    #  All covariates
    init_par = c(c(matrix(c(-2, 0, 0, 0, 0,
                            -2, 0, 0, 0, 0,
                            -2, 0, 0, 0, 0,
                            -2, 0, 0, 0, 0,
                            -2, 0, 0, 0, 0,
                            -2, 0, 0, 0, 0), ncol=5, byrow = T)),
             c(6.411967, 0, 0), 
             log(0.51^2),  
             c(log(0.8^2), 0, 0),
             0, 0, 0,
             0, 0, 0,
             -1, -1, -1, -1)
    par_index = list(zeta=1:30, misclass=44:47, delta = 31:33, tau2 = 34, 
                     sigma2 = 35:37, gamma = 38:40, delta_new = 41:43)
}

cov_info = temp_data[,c("Age", "sex1", "edu_yes", "DLER_avg"), drop=F] 
    
if(!simulation) {
    # Centering Age & DLER
    ages = NULL
    dler_val = NULL
    for(a in EIDs) {
        ages = c(ages, unique(data_format[data_format[,"ID.."] == a, "Age"]))
        dler_val = c(dler_val, unique(data_format[data_format[,"ID.."] == a, "DLER_avg"]))
    }
    mean_age = mean(ages)
    mean_dler = mean(dler_val)
    cov_info[,'Age'] = cov_info[,'Age'] - mean_age
    cov_info[,'DLER_avg'] = cov_info[,'DLER_avg'] - mean_dler    
}


ind = 0
chain_list = vector(mode = "list", length = length(index_seeds))
for(seed in index_seeds){
    
    file_name = NULL
    
    if(simulation) {
        file_name = paste0(dir,'mcmc_out_',toString(seed), '_', trial_num, '_sim_30b.rda') 
    } else {
        file_name = paste0(dir,'mcmc_out_',toString(seed), '_', trial_num, '_30b.rda')
    }
    
    if (file.exists(file_name)) {
        load(file_name)
        ind = ind + 1

        print(file_name)

        # Thinning the chain
        main_chain = mcmc_out$chain[index_post,]
        ind_keep = seq(1, nrow(main_chain), by=10)
        
        chain_list[[ind]] = main_chain[ind_keep, ]
    }
}

stacked_chains = do.call( rbind, chain_list)
init_par = apply(stacked_chains, 2, median)

# Understanding the states of interest:
non_baseline = matrix(nrow = length(EIDs), ncol = 2)
colnames(non_baseline) = c("EID", "num_non_base")
non_baseline[,1] = EIDs
for(i in 1:length(EIDs)) {
    y_1_sub = data_format[data_format[,"ID.."] == EIDs[i], "State"]
    non_baseline[i,2] = sum(y_1_sub != 1)
}

print("number of states needing to sample")
print(sort(unique(non_baseline[,2])))
print("total number of combinations for each state sequence")
print(3^sort(unique(non_baseline[,2])))
print("table of states")
print(table(non_baseline[,2]))

state_combos = vector(mode = 'list', length = 7)
curr_b = rep(1, 7)
for(a in 1:3) {
    curr_b[1] = a
    state_combos[[1]] = c(state_combos[[1]], curr_b[1])
    for(b in 1:3) {
        curr_b[2] = b
        state_combos[[2]] = c(state_combos[[2]], curr_b[1:2])
        for(c in 1:3) {
            curr_b[3] = c
            state_combos[[3]] = c(state_combos[[3]], curr_b[1:3])
            for(d in 1:3) {
                curr_b[4] = d
                state_combos[[4]] = c(state_combos[[4]], curr_b[1:4])
                for(e in 1:3) {
                    curr_b[5] = e
                    state_combos[[5]] = c(state_combos[[5]], curr_b[1:5])
                    for(f in 1:3) {
                        curr_b[6] = f
                        state_combos[[6]] = c(state_combos[[6]], curr_b[1:6])
                        for (g in 1:3) {
                            curr_b[7] = g
                            state_combos[[7]] = c(state_combos[[7]], curr_b[1:7])
                        }
                    }
                }
            }
        }
    }
}

for(i in 1:7) {
    temp_mat = matrix(state_combos[[i]], nrow = 3^i, byrow = T)
    state_combos[[i]] = temp_mat
}

B_MLE = vector(mode = 'list', length = length(EIDs))

for(ind in 1:length(EIDs)) {
    
    y_1_sub = y_1[data_format[,"ID.."] == EIDs[ind]]
    t_pts = c(max(which(y_1_sub == 1)), which(y_1_sub != 1)) - 1
    
    print(paste0(ind, ", ", EIDs[ind], ", num non-zero: ", non_baseline[ind, 2]))
    optimal_state_seq = brute_force_ss(EIDs[ind], ind - 1, init_par, par_index,
                                       id, y_2, y_1, cov_info, covariate_struct,
                                       state_combos, t_pts, length(EIDs))
    B_MLE[[ind]] = optimal_state_seq
}

if(simulation) {
    save(B_MLE, file = paste0("Model_out/B_MLE_", trial_num, "_sim.rda"))
} else {
    save(B_MLE, file = paste0("Model_out/B_MLE_", trial_num, ".rda"))   
}
