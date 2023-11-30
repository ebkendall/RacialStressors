load('Data/data_format_30.rda')
data_format = data_format_30
miss_info = c(26296, 29698, 30625, 401, 423, 419, 457)
data_format = data_format[!(data_format[,"ID.."] %in% miss_info), ]

temp_data = as.matrix(data_format); rownames(temp_data) = NULL

# Centering age
mean_age = mean(temp_data[,'Age'])
temp_data[,'Age'] = temp_data[,'Age'] - mean_age

par_index = list(delta = 1:3, tau2 = 4, sigma2 = 5:7, mu = 8)
par_index_og = list(delta = 1:3, tau2 = 4, sigma2 = 5:7)

likelihood_fnc <- function(par, data, par_index) {
      log_like = 0

      for(ii in 1:nrow(data)) {
            V_k_t = matrix(0, nrow = 3, ncol = 1)
            if(data[ii, "State"] == 1) {
                  V_k_t[1,1] = 1
            } else if(data[ii, "State"] == 2) {
                  V_k_t[1:2,1] = 1
            } else {
                  V_k_t[c(1,3),1] = 1
            }

            delta = matrix(par[par_index$delta],ncol = 1)
            sigma2_vec = exp(par[par_index$sigma2])
            tau2 = exp(par[par_index$tau2])

            big_Sigma_inv = diag(c(1/sigma2_vec))
            tau2_inv = 1/tau2

            W = tau2_inv * (V_k_t %*% t(V_k_t)) + big_Sigma_inv

            inv_W = solve(W)
            det_inv_W = det(inv_W)

            hold1 = (1/sqrt(2*pi*tau2)) * 
                       (1/sqrt(prod(sigma2_vec))) * sqrt(det_inv_W)

            x_ii = as.matrix(data[ii, c("Age", "sex1", "edu_yes", "DLER_avg"),drop=F])
            d_i = data[ii, "RSA"]

            temp1 = (d_i / tau2) * V_k_t + big_Sigma_inv %*% delta
            temp2 = t(temp1) %*% inv_W %*% temp1
            temp3 = t(delta) %*% big_Sigma_inv %*% delta

            exp_pow = (-1/(2*tau2)) * (d_i^2) - 0.5 * temp3 + 0.5 * temp2
            hold2 = exp(exp_pow)
            d_2_final = hold1 * hold2

            log_like = log_like + log(d_2_final)
      }
      return(log_like)
}

likelihood_fnc2 <- function(par, data, par_index) {
    log_like = 0
    
    delta = matrix(par[par_index$delta], ncol = 1)
    sigma2_vec = exp(par[par_index$sigma2])
    tau2 = exp(par[par_index$tau2])
    mu = par[par_index$mu]
    
    for(ii in 1:nrow(data)) {
        like_val = 0
        if(data[ii, "State"] == 1) {
            mean_ii = mu + delta[1]
            var_ii = tau2 + sigma2_vec[1]
            log_like = log_like + dnorm(x = data[ii, "RSA"], 
                                        mean = mean_ii,
                                        sd = sqrt(var_ii),
                                        log = T)
        } else if(data[ii, "State"] == 2) {
            mean_ii = mu + delta[2]
            var_ii = tau2 + sigma2_vec[2]
            log_like = log_like + dnorm(x = data[ii, "RSA"], 
                                        mean = mean_ii,
                                        sd = sqrt(var_ii),
                                        log = T)
        } else {
            mean_ii = mu + delta[3]
            var_ii = tau2 + sigma2_vec[3]
            log_like = log_like + dnorm(x = data[ii, "RSA"], 
                                        mean = mean_ii,
                                        sd = sqrt(var_ii),
                                        log = T)
        }
    }
    return(log_like)
}

likelihood_fnc_og <- function(par, data, par_index) {
    log_like = 0
    
    delta = matrix(par[par_index$delta], ncol = 1)
    sigma2_vec = exp(par[par_index$sigma2])
    tau2 = exp(par[par_index$tau2])
    
    for(ii in 1:nrow(data)) {
        like_val = 0
        if(data[ii, "State"] == 1) {
            mean_ii = delta[1]
            var_ii = tau2 + sigma2_vec[1]
            log_like = log_like + dnorm(x = data[ii, "RSA"], 
                                        mean = mean_ii,
                                        sd = sqrt(var_ii),
                                        log = T)
        } else if(data[ii, "State"] == 2) {
            mean_ii = delta[1] + delta[2]
            var_ii = tau2 + sigma2_vec[1] + sigma2_vec[2]
            log_like = log_like + dnorm(x = data[ii, "RSA"], 
                                        mean = mean_ii,
                                        sd = sqrt(var_ii),
                                        log = T)
        } else {
            mean_ii = delta[1] + delta[3]
            var_ii = tau2 + sigma2_vec[1] + sigma2_vec[3]
            log_like = log_like + dnorm(x = data[ii, "RSA"], 
                                        mean = mean_ii,
                                        sd = sqrt(var_ii),
                                        log = T)
        }
    }
    return(log_like)
}


init_par = c(         0, -0.2582171, -0.1166778, 
             -1.0842574,  
              0.1973717, -0.0714276,  0.2777276, 
              6.4557765)

numML2 <- optim(init_par, likelihood_fnc2, data=temp_data, par_index = par_index,
               control = list(fnscale=-1, trace = TRUE, maxit = 10000))
numML2$par

cat("State 1\n", "mean: ", numML2$par[1] + numML2$par[8], '\n', 
    "var: ", exp(numML2$par[4]) + exp(numML2$par[5]))

cat("State 2\n", "mean: ", numML2$par[2] + numML2$par[8], '\n', 
    "var: ", exp(numML2$par[4]) + exp(numML2$par[6]))

cat("State 3\n", "mean: ", numML2$par[3] + numML2$par[8], '\n', 
    "var: ", exp(numML2$par[4]) + exp(numML2$par[7]))

init_par2 = c( 6.4557765, -0.2582171, -0.1166778, 
              -1.0842574,  
               0.1973717, -0.0714276,  0.2777276)

numML_diff2 <- optim(init_par2, likelihood_fnc_og, data=temp_data, par_index = par_index_og,
                control = list(fnscale=-1, trace = TRUE, maxit = 10000))
numML_diff2$par

cat("State 1\n", "mean: ", numML_diff2$par[1], '\n', 
    "var: ", exp(numML_diff2$par[4]) + exp(numML_diff2$par[5]))

cat("State 2\n", "mean: ", numML_diff2$par[2], '\n', 
    "var: ", exp(numML_diff2$par[4]) + exp(numML_diff2$par[5]) + exp(numML_diff2$par[6]))

cat("State 3\n", "mean: ", numML_diff2$par[3], '\n', 
    "var: ", exp(numML_diff2$par[4]) + exp(numML_diff2$par[5]) + exp(numML_diff2$par[7]))
