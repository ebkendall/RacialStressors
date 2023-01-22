source("mcmc_routine.r")

args = commandArgs(TRUE)

ind = as.numeric(args[1]) 
set.seed(ind)
print(ind)

trial_num = 1

load('Data/data_format.rda')
n_sub = length(unique(data_format$ID..))

init_par = c(c(matrix(c( -11,  2.5,
                          -6, -1.7,
                          -7, 0.84,
                        -5.6, -1.7,
                        -5.2, -1.8), ncol=2, byrow = T)),
            c(-4, -4),
            c(6.411967, 6.481880, 6.335972), 
            1, 
            c(diag(3)),
            rep(1,273))

par_index = list( beta=1:10, misclass=11:12,
                  mu_tilde = 13:15, tau2 = 16, upsilon = 17:25,
                  mu_i = 26:298)
# par_index = list( beta=1:12, pi_logit=13:14,
#                   mu_tilde = 15:17, tau2 = 18, upsilon = 19:27,
#                   mu_i = 28:300)

# Initializing using the most recent MCMC -------------------------------------
load(paste0('Model_out/mcmc_out_2_10.rda'))
init_par[par_index$mu_i] = c(mcmc_out$M[[10]])
rm(mcmc_out)
# -----------------------------------------------------------------------------

prior_mean = c(c(matrix(c(-5, 0,
                          -5, 0,
                          -5, 0,
                          -5, 0,
                          -5, 0), ncol=2, byrow = T)))
            #    c(-5, -5))  
prior_sd = c(c(matrix(c(10, 5,
                        10, 5,
                        10, 5,
                        10, 5,
                        10, 5), ncol=2, byrow = T))) 
            #  c(5, 5))
prior_par = list()
prior_par[[1]] = prior_mean
prior_par[[2]] = prior_sd

temp_data = as.matrix(data_format); rownames(temp_data) = NULL
id = temp_data[,"ID.."]
y_1 = temp_data[,"State"]
y_2 = temp_data[,"RSA"]
t = temp_data[,"Time"]

steps = 10000
burnin = 5000
n_cores = 20

s_time = Sys.time()

mcmc_out = mcmc_routine(y_1, y_2, t, id, init_par, prior_par, par_index,
             steps, burnin, n_cores, n_sub)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, "_", trial_num, ".rda"))
