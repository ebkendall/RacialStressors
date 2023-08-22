library(mvtnorm, quietly=T)
library(foreach, quietly=T)
library(doParallel, quietly=T)
library(deSolve, quietly=T)
library(LaplacesDemon, quietly=T)

# Rcpp packages
library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)
sourceCpp("mcmc_routine_c.cpp")

# Needed for OpenMP C++ parallel
Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
Sys.setenv("PKG_LIBS" = "-fopenmp")


init_par = c(c(matrix(c(-4, 0,
                        -4, 0,
                        -4, 0,
                        -4, 0,
                        -4, 0), ncol=2, byrow = T)),
             c(-4, -4),
             c(6.411967, 6.481880, 6.335972), 
             1, 
             c(diag(3)),
             rep(1,273))

par_index = list( beta=1:10, misclass=11:12,
                  mu_tilde = 13:15, tau2 = 16, upsilon = 17:25,
                  mu_i = 26:298)

prior_mean = c(c(matrix(c(-5, 0,
                          -5, 0,
                          -5, 0,
                          -5, 0,
                          -5, 0), ncol=2, byrow = T)),
               c(-5, -5))  
prior_sd = c(c(matrix(c(10, 5,
                        10, 5,
                        10, 5,
                        10, 5,
                        10, 5), ncol=2, byrow = T)),
             c(5, 5))
prior_par = list()
prior_par[[1]] = prior_mean
prior_par[[2]] = prior_sd


test = test_functions(init_par, prior_par, par_index)

print(is.finite(test))

print((test < -100))


