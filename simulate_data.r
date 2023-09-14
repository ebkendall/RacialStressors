library(mvtnorm)

thirty = T

# Load the current data ------------------------------------------------------
if(thirty) {
  load('Data/data_format_30.rda')
  N = 500
  # N = length(unique(data_format_30$ID..))
  data_format = data_format_30
} else {
  load('Data/data_format_15.rda')
  N = length(unique(data_format_15$ID..))
  data_format = data_format_15
}

n_sim = 1

# True parameter values ------------------------------------------------------
par_index = list( zeta=1:5, misclass=6:9, delta = 10:12, tau2 = 13, sigma2 = 14)

if(thirty) {
    zeta = matrix(c(-1.98, -7, -0.37, -8.59, -8.72), ncol = 1)
    misclass = c(-8.4, -8.65, -7.23, -8.37)
    delta = c(6.49, -0.281, -0.114)
    log_tau2 = 0.422
    log_sigma2 = -4.8
} else {
    zeta = matrix(c(-2.78, -7.97, -1.484, -9.03, -9.16), ncol = 1)
    misclass = c(-9.17, -9.25, -7.99, -9.07)
    delta = c(6.43, -0.329, -0.115)
    log_tau2 = 0.534
    log_sigma2 = -4.89
}

true_par = c(c(zeta), misclass, delta, log_tau2, log_sigma2)
if(thirty) {
    save(true_par, file = 'Data/true_par_30.rda')   
} else {
    save(true_par, file = 'Data/true_par_15.rda')
}

# Simulate the data -----------------------------------------------------------

M = matrix(c(1, exp(misclass[1]), exp(misclass[2]),
             0,                1, exp(misclass[3]),
             0, exp(misclass[4]),                1), nrow = 3, byrow = T)
M = M / rowSums(M)

tau2 = exp(log_tau2)
sigma2 = exp(log_sigma2)

for(ind in 1:n_sim) {
    set.seed(ind)
    sim_data = NULL

    for(i in 1:N) {
        id  = i
        id_info = sample(x = unique(data_format[,"ID.."]), size = 1, replace = T)
        n_i = sum(data_format[,"ID.."] == id_info)
        b_i = NULL
        s_i = NULL
        t_pts = time = data_format[data_format[,"ID.."] == id_info, "Time"]
    
        for(k in 1:n_i) {
              if(k == 1) {
                b_i = 1
                s_i = 1
              } else {
                z_i = matrix(c(1), nrow=1)
    
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
    
                # Sample the observed state w/ misclassification probability
                s_i = c( s_i, sample(1:3, size=1, prob=M[tail(b_i,1),]))
              }
        }
    
        V_i = cbind(1, cbind(as.numeric(b_i == 2), as.numeric(b_i == 3)))
        # I have sampled both the true state sequence as well as the 
        # observed labels. Now we simulate the RSA sequence
    
        upsilon = diag(rep(sigma2, 3))
        delta_i = rmvnorm( n=1, mean= delta, sigma=upsilon)
    
        mean_Y_2 = V_i %*% matrix(delta_i, ncol = 1)
    
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
    
    print(sum(sim_data[,'State'] != sim_data[,'True_state']))

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
