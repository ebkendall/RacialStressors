source("mcmc_routine.r")

args = commandArgs(TRUE)

ind = as.numeric(args[1]) 
set.seed(ind)
print(ind)

trial_num = 10

load('Data/data_format.rda')
n_sub = length(unique(data_format$ID..))

init_par = c(c(matrix(c(-4, 0,
                        -4, 0,
                        -4, 0,
                        -4, 0,
                        -4, 0,
                        -4, 0), ncol=2, byrow = T)),
            c(-4, -4, -4, -4, -4, -4),
            c(-4, -4),
            c(6.411967, 6.481880, 6.335972), 
            1, 
            c(diag(3)))

par_index = list( beta=1:12, misclass = 13:18, pi_logit=19:20,
                  mu_tilde = 21:23, tau2 = 24, upsilon = 25:33)

# Initializing using the most recent MCMC -------------------------------------
load(paste0('Model_out/mcmc_out_2_9.rda'))
par_means = colMeans(mcmc_out$chain[4000:5000, ])
init_par[par_index$mu_tilde] = par_means[par_index$mu_tilde]
init_par[par_index$tau2] = par_means[par_index$tau2]
init_par[par_index$upsilon] = par_means[par_index$upsilon]
# init_par[par_index$tau2] = exp(init_par[par_index$tau2])
# init_par[par_index$mu_tilde] = c(6.411967, 6.481880, 6.335972)
rm(mcmc_out)
# -----------------------------------------------------------------------------

prior_mean = c(c(matrix(c(-5, 0,
                          -5, 0,
                          -5, 0,
                          -5, 0,
                          -5, 0,
                          -5, 0), ncol=2, byrow = T)),
               c(-5, -5, -5, -5, -5, -5),
               c(-5, -5))  
prior_sd = c(c(matrix(c(10, 5,
                        10, 5,
                        10, 5,
                        10, 5,
                        10, 5,
                        10, 5), ncol=2, byrow = T)),
             c(5, 5, 5, 5, 5, 5),
             c(5, 5))
prior_par = data.frame( prior_mean= prior_mean,
                        prior_sd= prior_sd)

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
