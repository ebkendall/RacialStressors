source("mcmc_routine.r")

args = commandArgs(TRUE)

ind = as.numeric(args[1]) 
set.seed(ind)
print(ind)

# Information defining which approach to take ----------------------------------
trial_num = 2
simulation = F
thirty = T
case_b = F
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
                          delta = 6:8, tau2 = 9, sigma2 = 10)
    } else {
        # Misclassification exists
        init_par = true_par
        
        par_index = list( zeta=1:5, misclass=6:9,
                          delta = 10:12, tau2 = 13, sigma2 = 14)
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
        init_par = c(c(matrix(c(-4,
                                -4,
                                -4,
                                -4,
                                -4), ncol=1, byrow = T)),
                     c(6.411967, 0, 0), 
                     log(0.51^2),  log(0.8^2))
        
        # misclass is kept in par_index for book-keeping in C++
        par_index = list( zeta=1:5, misclass=0,
                          delta = 6:8, tau2 = 9, sigma2 = 10)
    } else {
        # Misclassification exists
        init_par = c(c(matrix(c(-4, 0, 0, 0, 0, 0, 0,
                                -4, 0, 0, 0, 0, 0, 0,
                                -4, 0, 0, 0, 0, 0, 0,
                                -4, 0, 0, 0, 0, 0, 0,
                                -4, 0, 0, 0, 0, 0, 0), ncol=7, byrow = T)),
                     c(-4, -4, -4, -4),
                     c(6.411967, 0, 0), 
                     log(0.51^2),  log(0.8^2))   
        
        par_index = list( zeta=1:35, misclass=36:39,
                          delta = 40:42, tau2 = 43, sigma2 = 44)
    }
}

n_sub = length(unique(data_format[,'ID..']))

# Uninformed priors
prior_mean = rep(0, length(init_par))
prior_sd = rep(20, length(init_par))

prior_par = list()
prior_par[[1]] = prior_mean
prior_par[[2]] = prior_sd

temp_data = as.matrix(data_format); rownames(temp_data) = NULL
id = as.numeric(temp_data[,"ID.."])
y_1 = as.numeric(temp_data[,"State"])
y_2 = as.numeric(temp_data[,"RSA"])
t = as.numeric(temp_data[,"Time"])
cov_info = temp_data[,c("Age", "sex1", "sex2", "edu_yes", "edu_no", "sum_DLER"), drop=F]

steps = 500000
burnin = 5000

s_time = Sys.time()

mcmc_out = mcmc_routine(y_1, y_2, t, id, init_par, prior_par, par_index,
             steps, burnin, n_sub, case_b, cov_info)

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
