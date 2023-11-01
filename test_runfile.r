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


init_par = c(c(matrix(c(-4, 0, 0, 0, 0,
                        -4, 0, 0, 0, 0,
                        -4, 0, 0, 0, 0,
                        -4, 0, 0, 0, 0,
                        -4, 0, 0, 0, 0), ncol=5, byrow = T)),
             c(-4, -4, -4, -4),
             c(6.411967, 0, 0), 
             log(0.51^2), 
             c(log(0.8^2), log(0.3^2), log(0.3^2)),
             0, 0, 0, 0)   

par_index = list( zeta=1:25, misclass=26:29,
                  delta = 30:32, tau2 = 33, sigma2 = 34:36,
                  beta = 37:40)

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


