library(mvtnorm)

thirty = T

# Load the current data ------------------------------------------------------
load('Data/data_format_30.rda')
data_format = data_format_30
miss_info = c(26296, 29698, 30625, 401, 423, 419, 457)
data_format = data_format[!(data_format[,"ID.."] %in% miss_info), ]

# N = 500
N = length(unique(data_format$ID..))

n_sim = 7

# True parameter values ------------------------------------------------------
par_index = list(zeta=1:30, misclass=42:45, delta = 31:33, tau2 = 34, sigma2 = 35:37,
                 gamma = 38:41)

load('Model_out/mcmc_out_1_7_30b.rda')
true_par = colMeans(mcmc_out$chain[295000:495000, ])
# zeta = matrix(colMeans(mcmc_out$chain[,par_index$zeta]), ncol = 5)
# delta = c(6.46408805, -0.26810867, -0.11329740)
# log_tau2 = -1.08023658
# log_sigma2 = c(0.20535332, -0.05919292, 0.26003737)
# gamma = matrix(colMeans(mcmc_out$chain[,par_index$gamma]), ncol = 1)
# true_par = c(c(zeta), delta, log_tau2, log_sigma2, gamma)

save(true_par, file = paste0('Data/true_par_', n_sim, '_30.rda'))

# Simulate the data -----------------------------------------------------------

tau2 = exp(true_par[par_index$tau2])
Sigma2 = diag(exp(true_par[par_index$sigma2]))
zeta = matrix(true_par[par_index$zeta], nrow = 6)
delta = true_par[par_index$delta]
gamma = matrix(true_par[par_index$gamma], ncol = 1)

cov_info = data_format[,c("Age", "sex1", "edu_yes", "DLER_avg"), drop=F]
cov_info$Age = as.numeric(cov_info$Age)
cov_info$sex1 = as.numeric(cov_info$sex1)
cov_info$edu_yes = as.numeric(cov_info$edu_yes)
cov_info$DLER_avg = as.numeric(cov_info$DLER_avg)
cov_info = as.matrix(cov_info)

ages = NULL
dler_val = NULL
for(a in unique(data_format[,"ID.."])) {
    ages = c(ages, unique(data_format[data_format[,"ID.."] == a, "Age"]))
    dler_val = c(dler_val, unique(data_format[data_format[,"ID.."] == a, "DLER_avg"]))
}
mean_age = mean(ages)
mean_dler = mean(dler_val)
cov_info[,'Age'] = cov_info[,'Age'] - mean_age
cov_info[,'DLER_avg'] = cov_info[,'DLER_avg'] - mean_dler

ind = n_sim
set.seed(ind)
sim_data = NULL

EIDs = unique(data_format$ID..)
for(i in 1:N) {
    print(i)
    id  = i
    # id_info = sample(x = EIDs, size = 1, replace = T)
    id_info = EIDs[i]
    n_i = sum(data_format[,"ID.."] == id_info)
    b_i = NULL
    s_i = NULL
    t_pts = time = data_format[data_format[,"ID.."] == id_info, "Time"]
    cov_i = cov_info[data_format$ID.. ==id_info, ]
    sub_dat = data_format[data_format[,"ID.."] == id_info, ]
    y_1_sub = c(sub_dat[,"State"])

    start_run = FALSE
    b_i = 1
    s_i = 1
    for(k in 2:n_i) {
        if((y_1_sub[k-1] == 1) && (y_1_sub[k] != 1)) start_run = TRUE
        if(start_run) {
            z_i = matrix(c(1, cov_i[1,]), nrow=1)
            
            q1   = exp(z_i %*% t(zeta[1, , drop=F])) 
            q2   = exp(z_i %*% t(zeta[2, , drop=F]))
            q3   = exp(z_i %*% t(zeta[3, , drop=F]))
            q4   = exp(z_i %*% t(zeta[4, , drop=F]))
            q5   = exp(z_i %*% t(zeta[5, , drop=F]))
            q6   = exp(z_i %*% t(zeta[6, , drop=F]))
            
            # transitions: 1->2, 2->1, 2->3, 3->1, 3->2
            Q = matrix(c(  1,  q1, q2,
                           q3,   1, q4,
                           q5,  q6,  1), ncol=3, byrow=T)
            P_i = Q / rowSums(Q)
            
            # Sample the true, latent state sequence
            b_i = c( b_i, sample(1:3, size=1, prob=P_i[tail(b_i,1),]))
            s_i = c( s_i, 99)
        } else {
            b_i = c(b_i, 1)
            s_i = c(s_i, 1)
        }
    }

    rsa_i = NULL
    for(k in 1:n_i) {
        if(b_i[k] == 1) {
            rsa_i[k] = rnorm(n = 1, 
                             mean = delta[1] + cov_i[k,,drop = F] %*% gamma, 
                             sd = sqrt(tau2 + Sigma2[1,1]))
        } else if(b_i[k] == 2) {
            rsa_i[k] = rnorm(n = 1, 
                             mean = delta[1] + delta[2] + cov_i[k,,drop = F] %*% gamma, 
                             sd = sqrt(tau2 + Sigma2[2,2]))
        } else {
            rsa_i[k] = rnorm(n = 1, 
                             mean = delta[1] + delta[3] + cov_i[k,,drop = F] %*% gamma, 
                             sd = sqrt(tau2 + Sigma2[3,3]))
        }
    }

    sim_data_sub = cbind(rep(id, n_i), time, b_i, c(rsa_i), s_i, cov_i)

    sim_data = rbind(sim_data, sim_data_sub)

}

colnames(sim_data) = c("ID..", "Time", "State", "RSA", "Alt_state", 
                       "Age", "sex1", "edu_yes", "DLER_avg")

print(head(sim_data))

if(thirty) {
    save(sim_data, file = paste0("Data/sim_data_", ind, "_30.rda"))
} else {
    save(sim_data, file = paste0("Data/sim_data_", ind, "_15.rda"))   
}

cat('\n','Proption of occurances in each state:','\n')
print(table(sim_data[,'State'])/dim(sim_data)[1])
cat('\n')

# print(sum(sim_data[,'State'] != sim_data[,'True_state']))
# diff_ind = which(sim_data[,'State'] != sim_data[,'True_state'])
# if(length(diff_ind) > 0) {
#     diff_ind = c(diff_ind, diff_ind+1, diff_ind - 1)
#     diff_ind = sort(diff_ind)
#     print(diff_ind)
#     print(sim_data[diff_ind, c('State', 'True_state')])   
# }


# # Forcing some mislabels
# load('Data/Simulation/sim_data_1_b.rda')
# id = unique(sim_data[,"ID.."])
# sim_data = cbind(sim_data, 0)
# colnames(sim_data)[6] = "changed"
# 
# for(i in id) {
#     sub_dat = sim_data[sim_data[,"ID.."] == i, ]
#     
#     if(sum(sub_dat[,"State"]==3) > 8) {
#         sub_dat[,"changed"] = 1
#         start_ind = max(which(diff(sub_dat[,"State"]) == 1)) + 1
#         sub_dat[start_ind:(start_ind + 3), "State"] = 2
#     }
#     
#     sim_data[sim_data[,"ID.."] == i, ] = sub_dat
# }
# 
# save(sim_data, file = 'Data/Simulation/sim_data_1_c.rda')
