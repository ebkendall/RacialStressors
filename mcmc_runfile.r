source("mcmc_routine.r")

args = commandArgs(TRUE)

ind = as.numeric(args[1]) 
set.seed(ind)
print(ind)

trial_num = 1

simulation = T
init_par = NULL

if(simulation) {
    # Simulation
    load('Data/sim_data_1_a.rda')
    load('Data/true_par_a.rda')
    data_format = sim_data
    
    init_par = true_par
} else {
    # Real data analysis
    # load('Data/data_format_30.rda')
    # data_format = data_format_30
    load('Data/data_format_15.rda')
    data_format = data_format_15   
    
    init_par = c(c(matrix(c(-4,
                            -4,
                            -4,
                            -4), ncol=1, byrow = T)),
                 c(-4, -4),
                 c(6.411967, 0, 0), 
                 0.1,  0.1)
}

n_sub = length(unique(data_format[,'ID..']))

par_index = list( zeta=1:4, misclass=5:6,
                  delta = 7:9, tau2 = 10, sigma2 = 11)

prior_mean = rep(0 , length(init_par)) 
prior_sd   = rep(20, length(init_par))

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
    save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, "_", trial_num, ".rda"))
}
