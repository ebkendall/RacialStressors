library(mvtnorm)

# Model type -------------------------------------------------------------------
# 1: baseline only
# 2: baseline & DLER
# 3: all covariates

covariate_struct = 3
# ------------------------------------------------------------------------------

# Load the current data ------------------------------------------------------
load('../Data/data_format_30.rda')
data_format = data_format_30
miss_info = c(26296, 29698, 30625, 401, 423, 419, 457)
data_format = data_format[!(data_format[,"ID.."] %in% miss_info), ]

N = 3 * length(unique(data_format$ID..))

EIDs = unique(data_format$ID..)

# True parameter values ------------------------------------------------------
if(covariate_struct == 1) {
    # Baseline only model
    par_index = list(zeta=1:24, misclass=38:41, delta = 25:27, tau2 = 28, 
                     sigma2 = 29:31, gamma = 32:34, delta_new = 35:37)
    
    load(paste0("Model_out/mcmc_out_", 5, "_", 7, "_30b.rda"))
    true_par = apply(mcmc_out$chain[190000:200000, ], 2, median)
    save(true_par, file = paste0('Data/true_par_', covariate_struct, '_30.rda'))
    
} else if(covariate_struct == 2) {
    # Baseline & DLER model
    par_index = list(zeta=1:12, misclass=26:29, delta = 13:15, tau2 = 16, 
                     sigma2 = 17:19, gamma = 20:22, delta_new = 23:25)
    
    load(paste0("Model_out/mcmc_out_", 5, "_", 8, "_30b.rda"))
    true_par = apply(mcmc_out$chain[190000:200000, ], 2, median)
    save(true_par, file = paste0('Data/true_par_', covariate_struct, '_30.rda'))
    
} else {
    #  All covariates
    par_index = list(zeta=1:30, misclass=44:47, delta = 31:33, tau2 = 34, 
                     sigma2 = 35:37, gamma = 38:40, delta_new = 41:43)
    
    load(paste0("Model_out/mcmc_out_", 5, "_", 9, "_30b.rda"))
    true_par = apply(mcmc_out$chain[190000:200000, ], 2, median)
    print(true_par[par_index$zeta])
    true_par[par_index$zeta][1:2] = -2
    true_par[par_index$zeta][c(7,13,19,25)] = true_par[par_index$zeta][c(8,14,20,26)] 
    true_par[par_index$zeta][3:4] = -2
    true_par[par_index$zeta][c(10,16,22,28)] = true_par[par_index$zeta][c(9,15,21,27)] 
    true_par[par_index$zeta][5:6] = c(-4,-4)
    print(matrix(true_par[par_index$zeta], nrow = 6))
    save(true_par, file = paste0('Data/true_par_', covariate_struct, '_30.rda'))
}

# Covariate information -------------------------------------------------------
cov_mat_big = matrix(0, nrow=length(unique(data_format[,"ID.."])), ncol = 4)
for(i in 1:length(unique(data_format[,"ID.."]))) {
    sub_dat = as.matrix(data_format[data_format[,"ID.."] == unique(data_format[,"ID.."])[i], ])
    cov_mat_big[i, ] = c(sub_dat[1, c("Age", "sex1", "edu_yes", "DLER_avg")])
}
colnames(cov_mat_big) = c("Age", "sex1", "edu_yes", "DLER_avg")

# Centering Age & DLER
ages = mean(cov_mat_big[,"Age"])
dler_val = mean(cov_mat_big[,"DLER_avg"])

cov_mat_big[,'Age'] = cov_mat_big[,'Age'] - ages
cov_mat_big[,'DLER_avg'] = cov_mat_big[,'DLER_avg'] - dler_val    

# Simulate the data -----------------------------------------------------------
tau2 = exp(true_par[par_index$tau2])
sigma2_vec = exp(true_par[par_index$sigma2])
zeta = matrix(true_par[par_index$zeta], nrow = 6)
mu_alpha_beta = true_par[par_index$delta]
delta = true_par[par_index$delta_new]
gamma = matrix(true_par[par_index$gamma], ncol = 1)

print("Gamma")
print(c(gamma))
print("Delta")
print(delta)

for(ind in 1:100) {
    print(ind)
    set.seed(ind)
    sim_data = NULL

    for(i in 1:N) {
        id_info = sample(x = EIDs, size = 1, replace = T)
        n_i = sum(data_format[,"ID.."] == id_info)
        b_i = NULL
        s_i = NULL
        t_pts = time = data_format[data_format[,"ID.."] == id_info, "Time"]
        
        sub_dat = data_format[data_format[,"ID.."] == id_info, ,drop=F]
        y_1_sub = c(sub_dat[,"State"])
        
        cov_i_small = cov_mat_big[sample(x = 1:nrow(cov_mat_big), size = 1, replace = T), ,drop=F]
        cov_i = matrix(rep(cov_i_small, n_i), nrow = n_i, byrow = T)
        colnames(cov_i) = c("Age", "sex1", "edu_yes", "DLER_avg")
        cov_gamma_i = cov_i[, c("Age", "sex1", "edu_yes"), drop = F]
        cov_dler_i = c(cov_i[, "DLER_avg"])
        
        if(covariate_struct == 1) {
            z_i = matrix(c(1, cov_i[1,c("Age", "sex1", "edu_yes")]), nrow=1)   
        } else if(covariate_struct == 2) {
            z_i = matrix(c(1, cov_i[1,"DLER_avg"]), nrow=1)
        } else {
            z_i = matrix(c(1, cov_i[1,]), nrow=1)
        }

        start_run = TRUE
        b_i = 1
        s_i = rep(1, n_i)
        for(k in 2:n_i) {
            # if(y_1_sub[k] != 1 && y_1_sub[k-1] == 1) start_run = TRUE
            if(start_run) {
                q1   = exp(z_i %*% t(zeta[1, , drop=F])) 
                q2   = exp(z_i %*% t(zeta[2, , drop=F]))
                q3   = exp(z_i %*% t(zeta[3, , drop=F]))
                q4   = exp(z_i %*% t(zeta[4, , drop=F]))
                q5   = exp(z_i %*% t(zeta[5, , drop=F]))
                q6   = exp(z_i %*% t(zeta[6, , drop=F]))
                
                Q = matrix(c( 1,  q1, q2,
                             q3,   1, q4,
                             q5,  q6,  1), ncol=3, byrow=T)
                P_i = Q / rowSums(Q)
                
                # Sample the true, latent state sequence
                b_i = c( b_i, sample(1:3, size=1, prob=P_i[tail(b_i,1),]))
                s_i[k] = 99
            } else {
                b_i = c(b_i, 1)
            }
        }
        
        rsa_i = rep(0, n_i)
        
        n_i_s1 = sum(b_i == 1)
        n_i_s2 = sum(b_i == 2)
        n_i_s3 = sum(b_i == 3)

        if(n_i_s1 > 0) {
            mean_s1 = rep(mu_alpha_beta[1] 
                          + cov_gamma_i[1,,drop = F] %*% gamma
                          + cov_dler_i[1] * delta[1], n_i_s1)
            var_s1 = matrix(sigma2_vec[1], nrow = n_i_s1, ncol = n_i_s1)
            diag(var_s1) = diag(var_s1) + tau2
            rsa_i_s1 = rmvnorm(n = 1, mean = mean_s1, sigma = var_s1)
            
            rsa_i[b_i == 1] = rsa_i_s1
        }
        if(n_i_s2 > 0) {
            mean_s2 = rep(mu_alpha_beta[1] + mu_alpha_beta[2] 
                          + cov_gamma_i[1,,drop = F] %*% gamma
                          + cov_dler_i[1] * delta[2], n_i_s2)
            var_s2 = matrix(sigma2_vec[2], nrow = n_i_s2, ncol = n_i_s2)
            diag(var_s2) = diag(var_s2) + tau2
            rsa_i_s2 = rmvnorm(n = 1, mean = mean_s2, sigma = var_s2)
            
            rsa_i[b_i == 2] = rsa_i_s2
        }
        if(n_i_s3 > 0) {
            mean_s3 = rep(mu_alpha_beta[1] + mu_alpha_beta[3] 
                          + cov_gamma_i[1,,drop = F] %*% gamma
                          + cov_dler_i[1] * delta[3], n_i_s3)
            var_s3 = matrix(sigma2_vec[3], nrow = n_i_s3, ncol = n_i_s3)
            diag(var_s3) = diag(var_s3) + tau2
            rsa_i_s3 = rmvnorm(n = 1, mean = mean_s3, sigma = var_s3)
            
            rsa_i[b_i == 3] = rsa_i_s3
        }

        sim_data_sub = cbind(rep(i, n_i), time, b_i, c(rsa_i), s_i, cov_i)

        sim_data = rbind(sim_data, sim_data_sub)

    }

    colnames(sim_data) = c("ID..", "Time", "State", "RSA", "Alt_state", 
                           "Age", "sex1", "edu_yes", "DLER_avg")

    print('Proption of occurances in each state:')
    print(table(sim_data[,'State']))
    print(table(sim_data[,'State'])/dim(sim_data)[1])
    
    count_transitions = matrix(0, nrow=3, ncol=3)
    total_trans = 0
    for(i in unique(sim_data[,"ID.."])){
        
        b_i_mle = as.numeric(c(sim_data[sim_data[,"ID.."] == i, 'State']))
        
        for(t in 1:(length(b_i_mle) - 1)) {
            count_transitions[b_i_mle[t], b_i_mle[t+1]] = 
                count_transitions[b_i_mle[t], b_i_mle[t+1]] + 1
            total_trans = total_trans + 1
        }
    }
    
    print("Transition distribution")
    print(count_transitions)
    
    save(sim_data, file = paste0("Data/sim_data_", covariate_struct, "_", ind, ".rda"))
}