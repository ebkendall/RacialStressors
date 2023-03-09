source("mcmc_routine.r")

args = commandArgs(TRUE)

ind = as.numeric(args[1]) 
set.seed(ind)
print(ind)

trial_num = 28

# Real data analysis
# load('Data/data_format.rda')
load('Data/data_format_5.rda')
data_format = as.matrix(data_format_5)


# Simulation
# load('Data/Simulation/sim_data_1_c.rda')

n_sub = length(unique(data_format[,'id']))

# load('Data/new_delta_est.rda')
# load('Data/Simulation/true_par_b.rda')

init_par = c(c(matrix(c(-4,0,
                        -4,0,
                        -4,0,
                        -4,0,
                        -4,0), ncol=2, byrow = T)),
            c(-4, -4, -4, -4, -4, -4),
            c(6.411967, 0, 0), 
            1, 
            c(diag(3)),
            c(rep(6.411967, n_sub), rep(0, n_sub), rep(0,n_sub)))

par_index = list( zeta=1:10, misclass=11:16,
                  delta = 17:19, tau2 = 20, upsilon = 21:29,
                  delta_i = 30:413)
# par_index = list( zeta=1:5, misclass=6:11,
#                   delta = 12:14, tau2 = 15, upsilon = 16:24,
#                   delta_i = 25:297)
# par_index = list( zeta=1:8, misclass=9:14,
#                   delta = 15:17, tau2 = 18, upsilon = 19:27,
#                   delta_i = 28:300)

# Initializing using the most recent MCMC -------------------------------------
# load('Model_out/mcmc_out_1_26.rda')
# init_par[par_index$delta_i] = c(mcmc_out$big_delta_i[[20]])
# # init_par[-par_index$delta_i] = colMeans(mcmc_out$chain)
# init_par[-par_index$delta_i] = true_par
# rm(mcmc_out)
# -----------------------------------------------------------------------------

# prior_mean = c(c(matrix(c(-4,  2.1,
#                           -8, -1.7,
#                           -4,  1.8,
#                           -6, -1.7,
#                           -6, -1.7), ncol=2, byrow = T)),
#                c(-5, -5, -5, -5, -5, 0))

# prior_sd = c(c(matrix(c(2,2,
#                         2,2,
#                         2,2,
#                         2,2,
#                         2,2), ncol=2, byrow = T)),
#              c(5, 5, 5, 5, 5, 5))
prior_mean = rep(0,16)
prior_sd = rep(20, 16)

prior_par = list()
prior_par[[1]] = prior_mean
prior_par[[2]] = prior_sd

temp_data = as.matrix(data_format); rownames(temp_data) = NULL
id = as.numeric(temp_data[,"id"])
y_1 = as.numeric(temp_data[,"state"])
y_2 = as.numeric(temp_data[,"rsa"])
t = as.numeric(temp_data[,"time"])

steps = 20000
burnin = 5000
n_cores = 20

s_time = Sys.time()

mcmc_out = mcmc_routine(y_1, y_2, t, id, init_par, prior_par, par_index,
             steps, burnin, n_cores, n_sub)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, "_", trial_num, ".rda"))
