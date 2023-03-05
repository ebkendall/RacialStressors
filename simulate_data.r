library(mvtnorm)

# Load the current data ------------------------------------------------------
load('Data/Old_data/data_format.rda')
N = length(unique(data_format$ID..))
n_sim = 1

# True parameter values -------------------------------------------------------
par_index = list( zeta=1:8, misclass=9:14,
                  delta = 15:17, tau2 = 18, upsilon = 19:27,
                  delta_i = 28:300)
# par_index = list( zeta=1:4, misclass=5:10,
#                   delta = 11:13, tau2 = 14, upsilon = 15:23,
#                   delta_i = 24:296)

# Baseline only transition matrix
# Logit transition probability parameters
# 1->2, 2->3, 3->1, 3->2
zeta = matrix(c(   -3,  2.1,
                   -4,  2.1,
                   -6, -1.7,
                   -6, -1.7), ncol = 2, byrow = T)

# Logit misclassification probability parameters
misclass = c(   -6-4, -5-5,
             -2-8,    -2-8,
             -6-4, -6-4   )

M = matrix(c(               1, exp(misclass[1]), exp(misclass[2]),
             exp(misclass[3]),                1, exp(misclass[4]),
             exp(misclass[5]), exp(misclass[6]),               1), nrow = 3, byrow = T)
M = M / rowSums(M)

# Parent means for random effect
delta = c(6.6, -2, -0.5)

# Log of variance for parent parameters
tau2 = 0.1

# Covariance term for the random effect
upsilon = diag(c(1, 0.1^2, 0.1^2))

true_par = c(c(zeta), misclass, delta, tau2, c(upsilon))
save(true_par, file = 'Data/Simulation/true_par.rda')

# Simulate the data -----------------------------------------------------------

for(ind in 1:n_sim) {
      set.seed(ind)
      sim_data = NULL

      for(i in 1:N) {
            n_i = sum(data_format[,"ID.."] == unique(data_format[,"ID.."])[i])
            id  = unique(data_format[,"ID.."])[i]
            b_i = NULL
            s_i = NULL
            t_pts = 0:(n_i - 1)
            time = data_format[data_format[,"ID.."] == id, "Time"]

            for(k in 1:n_i) {
                  if(k == 1) {
                        b_i = 1
                        s_i = 1
                  }
                  else {
                        z_i = c(1, t_pts[k] / 10)

                        q1   = exp(z_i %*% t(zeta[1, , drop=F])) 
				q2   = exp(z_i %*% t(zeta[2, , drop=F]))
				q3   = exp(z_i %*% t(zeta[3, , drop=F]))
				q4   = exp(z_i %*% t(zeta[4, , drop=F]))

	                  # transitions: 1->2, 2->3, 3->1, 3->2
				Q = matrix(c(  1,  q1,  0,
				               0,   1, q2,
				              q3,  q4,  1), ncol=3, byrow=T)
				P_i = Q / rowSums(Q)

				# Sample the latent state sequence
				b_i = c( b_i, sample(1:3, size=1, prob=P_i[tail(b_i,1),]))

                        # Sample the misclassification probability
                        s_i = c( s_i, sample(1:3, size=1, prob=M[tail(b_i,1),]))
                  }
            }

            V_i = cbind(1, cbind(as.numeric(b_i == 2), as.numeric(b_i == 3)))
            # I have sampled both the true state sequence as well as the 
            # observed labels. Now we simulate the RSA sequence

            delta_i = rmvnorm( n=1, mean= delta, sigma=upsilon)

            mean_Y_2 = V_i %*% matrix(delta_i, ncol = 1)

            rsa_i = rmvnorm(n=1, mean = mean_Y_2, sigma = tau2 * diag(n_i))

            sim_data_sub = cbind(rep(id, n_i), time, s_i, c(rsa_i), b_i)

            sim_data = rbind(sim_data, sim_data_sub)

      }

      colnames(sim_data) = c("ID..", "Time", "State", "RSA", "True_state")
      save(sim_data, file = paste0("Data/Simulation/sim_data_", ind, ".rda"))

      cat('\n','Proption of occurances in each state:','\n')
	print(table(sim_data[,'True_state'])/dim(sim_data)[1])
	cat('\n')
	
	print(sum(sim_data[,'State'] != sim_data[,'True_state']))

}



