source("mcmc_routine.r")

args = commandArgs(TRUE)

ind = as.numeric(args[1]) 
set.seed(ind)
print(ind)

trial_num = 12

# Real data analysis
# load('Data/data_format_30.rda')
# data_format = data_format_30
load('Data/data_format_15.rda')
data_format = data_format_15

# Simulation
# load('Data/Simulation/sim_data_1_c.rda')
# load('Data/Simulation/true_par_b.rda')
# data_format = sim_data

n_sub = length(unique(data_format[,'ID..']))


init_par = c(c(matrix(c(-4,
                        -4,
                        -4,
                        -4,
                        -4), ncol=1, byrow = T)),
            c(-4, -4, -4, -4, -4, -4),
            c(6.411967, 0, 0), 
            1,  1)

par_index = list( zeta=1:5, misclass=6:11,
                  delta = 12:14, tau2 = 15, sigma2 = 16)

# Initializing using the most recent MCMC -------------------------------------
# load('Model_out/mcmc_out_3_7.rda')
# init_par[par_index$delta_i] = c(mcmc_out$big_delta_i[[20]])
# init_par[-par_index$delta_i] = colMeans(mcmc_out$chain[20000:25000, ])
# rm(mcmc_out)
# -----------------------------------------------------------------------------

# prior_mean = c(c(matrix(c(-4,
#                           -10,
#                           -4,
#                           -10,
#                           -10), ncol=1, byrow = T)),
#                c(-5, 0, 0, -5, -5, 0))
# 
# prior_sd = c(c(matrix(c(2,
#                         2,
#                         2,
#                         2,
#                         2), ncol=1, byrow = T)),
#              c(5, 5, 5, 5, 5, 5))

prior_mean = rep(0 ,16) # 11
prior_sd   = rep(20,16) # 11

prior_par = list()
prior_par[[1]] = prior_mean
prior_par[[2]] = prior_sd

temp_data = as.matrix(data_format); rownames(temp_data) = NULL
id = as.numeric(temp_data[,"ID.."])
y_1 = as.numeric(temp_data[,"State"])
y_2 = as.numeric(temp_data[,"RSA"])
t = as.numeric(temp_data[,"Time"])

steps = 30000
burnin = 5000
n_cores = 20

s_time = Sys.time()

mcmc_out = mcmc_routine(y_1, y_2, t, id, init_par, prior_par, par_index,
             steps, burnin, n_cores, n_sub)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, "_", trial_num, ".rda"))
