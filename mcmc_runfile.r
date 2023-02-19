source("mcmc_routine.r")

args = commandArgs(TRUE)

ind = as.numeric(args[1]) 
set.seed(ind)
print(ind)

trial_num = 19

load('Data/data_format.rda')
n_sub = length(unique(data_format$ID..))
load('Data/new_delta_est.rda')

init_par = c(c(matrix(c(-4,0,
                        -4,0,
                        -4,0,
                        -4,0), ncol=2, byrow = T)),
            c(-4, -4, -4, -4, -4, -4),
            c(6.411967, 0, 0), 
            1, 
            c(diag(3)),
            c(rep(6.411967, 91), rep(0, 91), rep(0,91)))

# par_index = list( zeta=1:10, misclass=11:16,
#                   delta = 17:19, tau2 = 20, upsilon = 21:29,
#                   delta_i = 30:302)
# par_index = list( zeta=1:5, misclass=6:11,
#                   delta = 12:14, tau2 = 15, upsilon = 16:24,
#                   delta_i = 25:297)
par_index = list( zeta=1:8, misclass=9:14,
                  delta = 15:17, tau2 = 18, upsilon = 19:27,
                  delta_i = 28:300)

# Initializing using the most recent MCMC -------------------------------------
# load(paste0('Model_out/mcmc_out_2_1.rda'))
# init_par[par_index$delta_i] = c(mcmc_out$big_mu_i[[10]])
# init_par[-par_index$delta_i] = colMeans(mcmc_out$chain)
# rm(mcmc_out)
# -----------------------------------------------------------------------------

prior_mean = rep(0, 8)
prior_sd = rep(20,8)

prior_par = list()
prior_par[[1]] = prior_mean
prior_par[[2]] = prior_sd

temp_data = as.matrix(data_format); rownames(temp_data) = NULL
id = temp_data[,"ID.."]
y_1 = temp_data[,"State"]
y_2 = temp_data[,"RSA"]
t = temp_data[,"Time"]

steps = 20000
burnin = 5000
n_cores = 20

s_time = Sys.time()

mcmc_out = mcmc_routine(y_1, y_2, t, id, init_par, prior_par, par_index,
             steps, burnin, n_cores, n_sub)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, "_", trial_num, ".rda"))
