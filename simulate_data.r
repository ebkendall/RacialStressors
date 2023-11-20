library(mvtnorm)

thirty = T

# Load the current data ------------------------------------------------------
load('Data/data_format_30.rda')
N = 500
# N = length(unique(data_format_30$ID..))
data_format = data_format_30

n_sim = 2

# True parameter values ------------------------------------------------------
par_index = list( zeta=1:25, misclass=0,
                  delta = 26:28, tau2 = 29, sigma2 = 30:32,
                  gamma = 33:36)

load('Model_out/mcmc_out_1_10_30b.rda')
zeta = matrix(colMeans(mcmc_out$chain[,par_index$zeta]), ncol = 5)
delta = c(6.46408805, -0.26810867, -0.11329740)
log_tau2 = -1.08023658
log_sigma2 = c(0.20535332, -0.05919292, 0.26003737)
gamma = matrix(colMeans(mcmc_out$chain[,par_index$gamma]), ncol = 1)

true_par = c(c(zeta), delta, log_tau2, log_sigma2, gamma)
save(true_par, file = paste0('Data/true_par_', n_sim, '_30.rda'))

# Simulate the data -----------------------------------------------------------

tau2 = exp(log_tau2)
Sigma2 = diag(exp(log_sigma2))

cov_info = data_format[,c("Age", "sex1", "edu_yes", "DLER_avg"), drop=F]
cov_info$Age = as.numeric(cov_info$Age)
cov_info$sex1 = as.numeric(cov_info$sex1)
cov_info$edu_yes = as.numeric(cov_info$edu_yes)
cov_info$DLER_avg = as.numeric(cov_info$DLER_avg)
cov_info = as.matrix(cov_info)

for(ind in 2:n_sim) {
    set.seed(ind)
    sim_data = NULL

    EIDs = unique(data_format$ID..)
    for(i in 1:N) {
        print(i)
        id  = i
        id_info = sample(x = EIDs, size = 1, replace = T)
        # id_info = EIDs[i]
        n_i = sum(data_format[,"ID.."] == id_info)
        b_i = NULL
        s_i = NULL
        t_pts = time = data_format[data_format[,"ID.."] == id_info, "Time"]
        cov_i = cov_info[data_format$ID.. ==id_info, ]
    
        for(k in 1:n_i) {
              if(k == 1) {
                b_i = 1
                s_i = 1
              } else {
                z_i = matrix(c(1, cov_i[1,]), nrow=1)
    
                q1   = exp(z_i %*% t(zeta[1, , drop=F])) 
                q2   = exp(z_i %*% t(zeta[2, , drop=F]))
                q3   = exp(z_i %*% t(zeta[3, , drop=F]))
                q4   = exp(z_i %*% t(zeta[4, , drop=F]))
                q5   = exp(z_i %*% t(zeta[5, , drop=F]))
          
                # transitions: 1->2, 2->1, 2->3, 3->1, 3->2
                Q = matrix(c(  1,  q1,  0,
                              q2,   1, q3,
                              q4,  q5,  1), ncol=3, byrow=T)
                P_i = Q / rowSums(Q)
          
                # Sample the true, latent state sequence
                b_i = c( b_i, sample(1:3, size=1, prob=P_i[tail(b_i,1),]))
              }
        }
        
        s_i = rep(99, length(b_i))
        jj = 1
        while(b_i[jj] == 1) {
            s_i[jj] = 1
            jj = jj+1
            if(jj > length(b_i)) break
        }
    
        V_i = cbind(1, cbind(as.numeric(b_i == 2), as.numeric(b_i == 3)))
        X_i = matrix(1, nrow = nrow(V_i), ncol=1) %x% cov_i[1,,drop = F]
    
        delta_i = rmvnorm( n=1, mean= delta, sigma=Sigma2)
    
        mean_Y_2 = V_i %*% matrix(delta_i, ncol = 1) + X_i %*% gamma
    
        rsa_i = rmvnorm(n=1, mean = mean_Y_2, sigma = diag(rep(tau2,n_i)) )
    
        sim_data_sub = cbind(rep(id, n_i), time, s_i, c(rsa_i), b_i)
    
        sim_data = rbind(sim_data, sim_data_sub)
    
    }

    colnames(sim_data) = c("ID..", "Time", "State", "RSA", "True_state")
    if(thirty) {
        save(sim_data, file = paste0("Data/sim_data_", ind, "_30.rda"))
    } else {
        save(sim_data, file = paste0("Data/sim_data_", ind, "_15.rda"))   
    }
    
    cat('\n','Proption of occurances in each state:','\n')
    print(table(sim_data[,'True_state'])/dim(sim_data)[1])
    cat('\n')
    
    # print(sum(sim_data[,'State'] != sim_data[,'True_state']))
    # diff_ind = which(sim_data[,'State'] != sim_data[,'True_state'])
    # if(length(diff_ind) > 0) {
    #     diff_ind = c(diff_ind, diff_ind+1, diff_ind - 1)
    #     diff_ind = sort(diff_ind)
    #     print(diff_ind)
    #     print(sim_data[diff_ind, c('State', 'True_state')])   
    # }

}


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
