source("mcmc_routine.r")

args = commandArgs(TRUE)

ind = as.numeric(args[1]) 
set.seed(ind)
print(ind)

trial_num = 3

simulation = F
thirty = T
init_par = NULL

if(simulation) {
    # Simulation
    load('Data/sim_data_1_a.rda')
    load('Data/true_par_a.rda')
    data_format = sim_data
    
    init_par = true_par
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
    
    init_par = c(c(matrix(c(-4,
                            -4,
                            -4,
                            -4), ncol=1, byrow = T)),
                 c(-4, -1, -1, -4, -4, -1),
                 c(6.411967, 0, 0), 
                 log(0.51^2),  log(0.8^2),
                 -4, -4)
}

n_sub = length(unique(data_format[,'ID..']))

par_index = list( zeta=1:4, misclass=5:10,
                  delta = 11:13, tau2 = 14, sigma2 = 15,
                  init = 16:17)

prior_mean = c(-4, -4, -4, -4,
                0,  0,  0,  0,  0,  0,
               6.5, -2, -0.5,
               -1.386, -0.45,
               -4, -4) 
prior_sd   = c(20, 20, 20, 20,
               1, .1, .1, 1, 1, .1,
               1, 0.5^2, 0.1^1,
               0.1^2, 0.1^2,
               5, 5)

prior_par = list()
prior_par[[1]] = prior_mean
prior_par[[2]] = prior_sd

temp_data = as.matrix(data_format); rownames(temp_data) = NULL
id = as.numeric(temp_data[,"ID.."])
y_1 = as.numeric(temp_data[,"State"])
y_2 = as.numeric(temp_data[,"RSA"])
t = as.numeric(temp_data[,"Time"])

steps = 20000
burnin = 5000

s_time = Sys.time()

mcmc_out = mcmc_routine(y_1, y_2, t, id, init_par, prior_par, par_index,
             steps, burnin, n_sub)

e_time = Sys.time() - s_time; print(e_time)

if(simulation) {
    save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, "_", trial_num, "_sim.rda"))    
} else {
    if(thirty) {
        save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, "_", trial_num, "_30.rda"))
    } else {
        save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, "_", trial_num, "_15.rda"))
    }
}
