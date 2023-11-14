source("mcmc_routine.r")

args = commandArgs(TRUE)

ind = as.numeric(args[1]) 
set.seed(ind)
print(ind)

# Information defining which approach to take ----------------------------------
trial_num = 10
simulation = F
thirty = T
case_b = T
# ------------------------------------------------------------------------------

init_par = NULL

if(simulation) {
    # Simulation
    if(thirty) {
        # 30s epochs
        load('Data/sim_data_1_30.rda')
        load('Data/true_par_30.rda')
        data_format = sim_data
    } else {
        # 15s epochs
        load('Data/sim_data_1_15.rda')
        load('Data/true_par_15.rda')
        data_format = sim_data
    }
    
    # Initialize the chains at the true values
    if(case_b) {
        # No misclassification because NO y_1
        init_par = c(true_par[1:5], true_par[10:14])
        
        # misclass is kept in par_index for book-keeping in C++
        par_index = list( zeta=1:5, misclass=0,
                          delta = 6:8, tau2 = 9, sigma2 = 10:12)
    } else {
        # Misclassification exists
        init_par = true_par
        
        par_index = list( zeta=1:5, misclass=6:9,
                          delta = 10:12, tau2 = 13, sigma2 = 14:16)
    }
    
} else {
    # Real data analysis
    if(thirty) {
        # 30s epochs
        load('Data/data_format_30.rda')
        data_format = data_format_30
    } else {
        # 15s epochs
        load('Data/data_format_15.rda')
        data_format = data_format_15   
    }
    
    # Removing the participants with missing labels
    miss_info = c(26296, 29698, 30625, 401, 423, 419, 457)
    data_format = data_format[!(data_format[,"ID.."] %in% miss_info), ]
    
    # Initialize the chains are arbitrary starting values
    if(case_b) {
        # No misclassification because NO y_1
        init_par = c(c(matrix(c(-4, 0, 0, 0, 0,
                                -4, 0, 0, 0, 0,
                                -4, 0, 0, 0, 0,
                                -4, 0, 0, 0, 0,
                                -4, 0, 0, 0, 0), ncol=5, byrow = T)),
                     c(6.411967, 0, 0), 
                     log(0.51^2),  
                     c(log(0.8^2), 0, 0),
                     0, 0, 0, 0)
        
        # misclass is kept in par_index for book-keeping in C++
        par_index = list( zeta=1:25, misclass=0,
                          delta = 26:28, tau2 = 29, sigma2 = 30:32,
                          gamma = 33:36)
    } else {
        # Misclassification exists
        init_par = c(c(matrix(c(-4, 0, 0, 0, 0,
                                -4, 0, 0, 0, 0,
                                -4, 0, 0, 0, 0,
                                -4, 0, 0, 0, 0,
                                -4, 0, 0, 0, 0), ncol=5, byrow = T)),
                     c(-4, -4, -4, -4),
                     c(6.411967, 0, 0), 
                     log(0.51^2), 
                     c(log(0.8^2), 0, 0),
                     0, 0, 0, 0)   
        
        par_index = list( zeta=1:25, misclass=26:29,
                          delta = 30:32, tau2 = 33, sigma2 = 34:36,
                          gamma = 37:40)
    }
}

n_sub = length(unique(data_format[,'ID..']))

# Initializing the variance and mean terms ------------------------------------
s1_group = data_format[data_format$State == 1, ]
s2_group = data_format[data_format$State == 2, ]
s3_group = data_format[data_format$State == 3, ]

variance_calc <- function(s_group) {
    y_bar = mean(s_group$RSA)
    tau2_hat = 0
    tau_sigma_hat = 0
    for(i in unique(s_group$ID..)) {
        
        sub_i = s_group[s_group$ID.. == i, ]
        y_bar_i = mean(sub_i$RSA)
        for(j in 1:nrow(sub_i)) {
            tau2_hat = tau2_hat + (sub_i$RSA[j] - y_bar_i)^2
        }
        
        tau_sigma_hat = tau_sigma_hat + nrow(sub_i) * (y_bar_i - y_bar)^2
    }
    big_N = nrow(s_group)
    little_t = length(unique(s_group$ID..))
    n_i_vec = c(table(s_group$ID..))
    names(n_i_vec) = NULL
    n_0 = (1/(little_t - 1)) * (big_N - sum(n_i_vec^2) / big_N)
    
    tau2_hat = (1 / (big_N - little_t)) * tau2_hat
    sigma2_hat = (1 / n_0) * ((1/(little_t - 1)) * tau_sigma_hat - tau2_hat)
    return(c(tau2_hat, sigma2_hat, y_bar))
}

s1_vars = variance_calc(s1_group)
s2_vars = variance_calc(s2_group)
s3_vars = variance_calc(s3_group)
all_vars = variance_calc(data_format)

init_par[par_index$tau2] = log(s1_vars[1])
init_par[par_index$sigma2[1]] = log(s1_vars[2])
init_par[par_index$sigma2[2]] = log(s2_vars[2])
init_par[par_index$sigma2[3]] = log(s3_vars[2])
init_par[par_index$delta[1]] = s1_vars[3]
init_par[par_index$delta[2]] = s2_vars[3] - s1_vars[3]
init_par[par_index$delta[3]] = s3_vars[3] - s1_vars[3]
# ------------------------------------------------------------------------------

# Specifying the priors --------------------------------------------------------
prior_mean = rep(0, length(init_par))
prior_sd = rep(20, length(init_par))
if(case_b) {
    prior_mean[26:28] = c(s1_vars[3], 
                          s2_vars[3] - s1_vars[3], 
                          s3_vars[3] - s1_vars[3])
    
    prior_mean[29] = log(s1_vars[1])
    prior_sd[29] = 1
    prior_mean[30] = log(s1_vars[2])
    prior_sd[30] = 1
    prior_mean[31] = log(s2_vars[2])
    prior_sd[31] = 1
    prior_mean[32] = log(s3_vars[2])
    prior_sd[32] = 1
}

prior_par = list()
prior_par[[1]] = prior_mean
prior_par[[2]] = prior_sd
# ------------------------------------------------------------------------------

temp_data = as.matrix(data_format); rownames(temp_data) = NULL
id = as.numeric(temp_data[,"ID.."])
y_1 = as.numeric(temp_data[,"State"])
y_2 = as.numeric(temp_data[,"RSA"])
t = as.numeric(temp_data[,"Time"])
if(simulation) {
    cov_info = matrix(0, nrow = nrow(temp_data), ncol = 4)
} else {
    cov_info = temp_data[,c("Age", "sex1", "edu_yes", "DLER_avg"), drop=F]
}

steps = 50000
burnin = 5000

s_time = Sys.time()

print(init_par)
mcmc_out = mcmc_routine(y_1, y_2, t, id, init_par, prior_par, par_index,
             steps, burnin, n_sub, case_b, cov_info, simulation)

e_time = Sys.time() - s_time; print(e_time)

if(simulation) {
    if(thirty) {
        if(case_b) {
            save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, "_", trial_num, "_sim_30b.rda"))  
        } else {
            save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, "_", trial_num, "_sim_30.rda"))     
        }
    } else {
        if(case_b) {
            save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, "_", trial_num, "_sim_15b.rda"))     
        } else {
            save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, "_", trial_num, "_sim_15.rda"))        
        }
    }
} else {
    if(thirty) {
        if(case_b) {
            save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, "_", trial_num, "_30b.rda"))
        } else {
            save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, "_", trial_num, "_30.rda"))   
        }
    } else {
        if(case_b) {
            save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, "_", trial_num, "_15b.rda"))
        } else {
            save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, "_", trial_num, "_15.rda"))   
        }
    }
}
