source("mcmc_routine.r")

args = commandArgs(TRUE)

ind = as.numeric(args[1]) 
set.seed(ind)
print(ind)

# Information defining which approach to take ----------------------------------
trial_num = 2
simulation = T
case_b = F
# ------------------------------------------------------------------------------

init_par = NULL

if(simulation) {
    # Simulation
    # 30s epochs
    load('Data/sim_data_3_30.rda')
    load('Data/true_par_3_30.rda')
    data_format = sim_data
} else {
    # Real data analysis
    # 30s epochs
    load('Data/data_format_30.rda')
    data_format = data_format_30
    
    # Removing the participants with missing labels
    miss_info = c(26296, 29698, 30625, 401, 423, 419, 457)
    data_format = data_format[!(data_format[,"ID.."] %in% miss_info), ]
}

init_par = c(c(matrix(c(-4, 0, 0, 0, 0,
                        -4, 0, 0, 0, 0,
                        -4, 0, 0, 0, 0,
                        -4, 0, 0, 0, 0,
                        -4, 0, 0, 0, 0,
                        -4, 0, 0, 0, 0), ncol=5, byrow = T)),
                c(0, 0, 0), 
                log(0.51^2),  
                c(log(0.8^2), 0, 0),
                6.411967, 0, 0, 0, 0)
    
# misclass is kept in par_index for book-keeping in C++
par_index = list(zeta=1:30, misclass=0,delta = 31:33, tau2 = 34, sigma2 = 35:37,
                 gamma = 38:42)

n_sub = length(unique(data_format[,'ID..']))

# Initializing the variance and mean terms ------------------------------------
if(simulation) {
    init_par = true_par
} else {
    s1_group = data_format[data_format[,'State'] == 1, ]
    s2_group = data_format[data_format[,'State'] == 2, ]
    s3_group = data_format[data_format[,'State'] == 3, ]
    
    alt_var_calc <- function(s_group) {
        x_bar = mean(s_group[,'RSA'])
        ss = rep(0, nrow(s_group))
        for(i in 1:nrow(s_group)) {
            ss[i] = (s_group[i,'RSA'] - x_bar)^2
        }
        ss = mean(ss)
        cat("mean: ", x_bar, '\n', 
            "var: ", ss, '\n')
    }
    print("Empirical MLE Estimates")
    print("State 1:"); alt_var_calc(s1_group)
    print("State 2:"); alt_var_calc(s2_group)
    print("State 3:"); alt_var_calc(s3_group)
    
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
    
    print("Empirical Sums of Squares Estimates")
    print("State 1:")
    cat("mean: ", s1_vars[3], '\n', 
        "var: ", s1_vars[1] + s1_vars[2], '\n')
    print("State 2:")
    cat("mean: ", s2_vars[3], '\n', 
        "var: ", s2_vars[1] + s2_vars[2], '\n')
    print("State 3:")
    cat("mean: ", s3_vars[3], '\n', 
        "var: ", s3_vars[1] + s3_vars[2], '\n')
    
    init_par[par_index$tau2] = log(s1_vars[1])
    init_par[par_index$sigma2[1]] = log(s1_vars[2])
    init_par[par_index$sigma2[2]] = log(s2_vars[2])
    init_par[par_index$sigma2[3]] = log(s3_vars[2])
    init_par[par_index$delta[1]] = 0
    init_par[par_index$delta[2]] = s2_vars[3] - s1_vars[3]
    init_par[par_index$delta[3]] = s3_vars[3] - s1_vars[3]   
    init_par[par_index$gamma[1]] = s1_vars[3]
}
# ------------------------------------------------------------------------------

# Specifying the priors --------------------------------------------------------
prior_mean = rep(0, length(init_par))
prior_sd = rep(20, length(init_par))

prior_par = list()
prior_par[[1]] = prior_mean
prior_par[[2]] = prior_sd
# ------------------------------------------------------------------------------

temp_data = as.matrix(data_format); rownames(temp_data) = NULL
id = as.numeric(temp_data[,"ID.."])
y_1 = as.numeric(temp_data[,"State"])
y_2 = as.numeric(temp_data[,"RSA"])
t = as.numeric(temp_data[,"Time"])

if(simulation & case_b) {
    y_1 = as.numeric(temp_data[,"Alt_state"])
}

cov_info = temp_data[,c("Age", "sex1", "edu_yes", "DLER_avg"), drop=F]

# Centering age
if(!simulation) {
    mean_age = mean(cov_info[,'Age'])
    cov_info[,'Age'] = cov_info[,'Age'] - mean_age   
}

steps = 100000
burnin = 5000

s_time = Sys.time()

print(init_par)
mcmc_out = mcmc_routine(y_1, y_2, t, id, init_par, prior_par, par_index,
             steps, burnin, n_sub, case_b, cov_info, simulation)

e_time = Sys.time() - s_time; print(e_time)

if(simulation) {
    if(case_b) {
        save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, "_", trial_num, "_sim_30b.rda"))  
    } else {
        save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, "_", trial_num, "_sim_30.rda"))     
    }
} else {
    if(case_b) {
        save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, "_", trial_num, "_30b.rda"))
    } else {
        save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, "_", trial_num, "_30.rda"))   
    }
}
