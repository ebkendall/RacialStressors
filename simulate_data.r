library(mvtnorm)

thirty = T

# Load the current data ------------------------------------------------------
if(thirty) {
  load('Data/data_format_30.rda')
  N = length(unique(data_format_30$ID..))
  data_format = data_format_30
} else {
  load('Data/data_format_15.rda')
  N = length(unique(data_format_15$ID..))
  data_format = data_format_15
}

n_sim = 1

# True parameter values ------------------------------------------------------
# Can only mis-classify states 2 and 3
par_index = list( zeta=1:4, misclass=5:6,
                  delta = 7:9, tau2 = 10, sigma2 = 11)

# Baseline only transition matrix
# Logit transition probability parameters
# 1->2, 2->3, 3->1, 3->2
zeta = matrix(c(   -2,#  2.1,
                   -1,#  1.8,
                   -8.5,# -1.7,
                   -8.6), ncol = 1, byrow = T)# -1.7

# Logit misclassification probability parameters
misclass = c(-4, -4)

M = matrix(c(1,                0,                0,
             0,                1, exp(misclass[1]),
             0, exp(misclass[2]),                1), nrow = 3, byrow = T)
M = M / rowSums(M)

# Parent means for random effect
delta = c(6.6, -2, -0.5)

# Log of variance for parent parameters
tau2 = 0.5

# Covariance term for the random effect
sigma2 = 0.01

true_par = c(c(zeta), misclass, delta, log(tau2), log(sigma2))
save(true_par, file = 'Data/true_par_a.rda')

# Simulate the data -----------------------------------------------------------

for(ind in 1:n_sim) {
    set.seed(ind)
    sim_data = NULL

    for(i in 1:N) {
        n_i = sum(data_format[,"ID.."] == unique(data_format[,"ID.."])[i])
        id  = unique(data_format[,"ID.."])[i]
        b_i = NULL
        s_i = NULL
        t_pts = time = data_format[data_format[,"ID.."] == id, "Time"]
    
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
          
                      # transitions: 1->2, 2->3, 3->1, 3->2
                Q = matrix(c(  1,  q1,  0,
                              0,   1, q2,
                              q3,  q4,  1), ncol=3, byrow=T)
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
    save(sim_data, file = paste0("Data/sim_data_", ind, "_a.rda"))
    
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
