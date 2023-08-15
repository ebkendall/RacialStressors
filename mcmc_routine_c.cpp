#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

const double pi = 3.14159265358979323846;

const arma::mat adj_mat = { {1, 1, 0},
                            {0, 1, 1},
                            {1, 1, 1}};

//  FUNCTIONS: ---------------------------------------------------------------

double D_2_calc(const int state, const double y_2_k, const double tau2, 
                const double sigma2, const arma::vec delta) {
    arma::vec V_k_t;
    if(state == 1) {
        V_k_t = {1,0,0};
    } else if(state == 2) {
        V_k_t = {1,1,0};
    } else {
        V_k_t = {1,0,1};
    }
    
    arma::mat I = arma::eye(3,3);
    arma::mat W = (1/tau2) * (V_k_t * V_k_t.t()) + (1/sigma2) * I;
    arma::mat inv_W = arma::inv(W);
    double det_inv_W = arma::det(inv_W);
    
    double hold1 = (1/sqrt(2*pi*tau2)) * (1/sqrt(sigma2 * sigma2 * sigma2)) * sqrt(det_inv_W);
    
    arma::vec temp1 = (y_2_k / tau2) * V_k_t + (1/sigma2) * delta;
    double temp2 = arma::as_scalar(temp1.t() * inv_W * temp1);
    double temp3 = arma::as_scalar(delta.t() * delta);
    
    double exp_pow = (-1/(2*tau2)) * (y_2_k * y_2_k) - (1/(2*sigma2)) * temp3 + 0.5 * temp2;
    double hold2 = exp(exp_pow);
    
    double d_2_final = hold1 * hold2;
    return d_2_final;
}

// [[Rcpp::export]]
double fn_log_post_continuous(const arma::vec &EIDs, const arma::vec &pars,  
                              const arma::field<arma::vec> &prior_par, const arma::field<arma::uvec> &par_index,
                              const arma::vec &y_1, const arma::vec &id, const arma::vec &y_2) {
    
    // par_index KEY: (0) zeta, (1) misclass, (2) delta, (3) tau2, (4) sigma2, (5) init
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    arma::vec in_vals(EIDs.n_elem, arma::fill::zeros);
    
    // Initial state probabilities
    arma::vec vec_init_content = pars.elem(par_index(5) - 1);
    arma::vec init_logit = {1, exp(vec_init_content(0)), exp(vec_init_content(1))};
    arma::vec init = init_logit / arma::accu(init_logit);
    
    // Manually populate the misclassification probabilities
    arma::vec vec_misclass_content = pars.elem(par_index(1) - 1);
    arma::mat M = { {1, exp(vec_misclass_content(0)), exp(vec_misclass_content(1))},
                    {exp(vec_misclass_content(2)), 1, exp(vec_misclass_content(3))},
                    {exp(vec_misclass_content(4)), exp(vec_misclass_content(5)), 1}};
    arma::vec m_row_sums = arma::sum(M, 1);
    M = M.each_col() / m_row_sums;
    
    // Populate the Transition probability matrix (independent of time)
    arma::vec vec_zeta_content = pars.elem(par_index(0) - 1);
    arma::mat zeta = arma::reshape(vec_zeta_content, 4, 1); 
    // arma::colvec z_i = {1, k_scale}; // using the current time point
    arma::colvec z_i = {1};
    
    double q1_sub = arma::as_scalar(zeta.row(0) * z_i);
    double q1 = exp(q1_sub);
    double q2_sub = arma::as_scalar(zeta.row(1) * z_i);
    double q2 = exp(q2_sub);
    double q3_sub = arma::as_scalar(zeta.row(2) * z_i);
    double q3 = exp(q3_sub);
    double q4_sub = arma::as_scalar(zeta.row(3) * z_i);
    double q4 = exp(q4_sub);
    
    arma::mat Q = { {  1,  q1,  0},
                    {  0,   1,  q2},
                    { q3,  q4,   1}};
    arma::vec q_row_sums = arma::sum(Q, 1);
    arma::mat P = Q.each_col() / q_row_sums;
    
    arma::vec delta = pars.elem(par_index(2) - 1);
    
    double log_tau2 = arma::as_scalar(pars.elem(par_index(3) - 1));
    double tau2 = exp(log_tau2);

    double log_sigma2 = arma::as_scalar(pars.elem(par_index(4) - 1));
    double sigma2 = exp(log_sigma2);
    
    omp_set_num_threads(6);
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        int i = EIDs(ii);
        
        // Subsetting the data
        arma::uvec sub_ind = arma::find(id == i);
        arma::vec y_1_i = y_1.elem(sub_ind);
        arma::vec y_2_i = y_2.elem(sub_ind);
        
        // Likelihood component from y_1
        arma::vec misclass_fill = M.col(y_1_i(0) - 1);
        arma::mat D_i_1 = arma::diagmat(misclass_fill);
        
        // Likelihood component from y_2
        double d_1 = D_2_calc(1, y_2_i(0), tau2, sigma2, delta);
        double d_2 = D_2_calc(2, y_2_i(0), tau2, sigma2, delta);
        double d_3 = D_2_calc(3, y_2_i(0), tau2, sigma2, delta);
        
        arma::vec d_fill = {d_1, d_2, d_3};
        arma::mat D_i_2 = arma::diagmat(d_fill);
        
        arma::mat init_transpose = init.t();
        
        arma::mat f_i = init_transpose * D_i_1 * D_i_2;
        
        for(int k = 1; k < y_1_i.n_elem; k++) {
            
            // Likelihood component from y_1
            arma::vec misclass_fill = M.col(y_1_i(k) - 1);
            arma::mat D_i_1 = arma::diagmat(misclass_fill);
            
            // Likelihood component from y_2
            double d_1 = D_2_calc(1, y_2_i(k), tau2, sigma2, delta);
            double d_2 = D_2_calc(2, y_2_i(k), tau2, sigma2, delta);
            double d_3 = D_2_calc(3, y_2_i(k), tau2, sigma2, delta);
            
            arma::vec d_fill = {d_1, d_2, d_3};
            arma::mat D_i_2 = arma::diagmat(d_fill);
            
            arma::mat temp = f_i * P * D_i_1 * D_i_2;
            f_i = temp;
        }
        
        in_vals(ii) = log(arma::accu(f_i));
    }
    
    double in_value = arma::accu(in_vals);
    
    // Likelihood components from the Metropolis priors
    arma::vec p_mean = prior_par(0);
    arma::mat p_sd = arma::diagmat(prior_par(1));

    arma::mat x = pars;
    double log_prior_dens = arma::as_scalar(dmvnorm(x.t(), p_mean, p_sd, true));
    
    in_value = in_value + log_prior_dens;
    
    return in_value;
}


// [[Rcpp::export]]
void test_functions(const arma::vec &pars, const arma::field<arma::vec> &prior_par, 
                    const arma::field<arma::uvec> &par_index) {

    arma::vec delta = {1,2,3};
    double temp1 = D_2_calc(1, 1, 1, 1, delta);
    double temp2 = D_2_calc(2, 1, 1, 1, delta);
    double temp3 = D_2_calc(3, 1, 1, 1, delta);
    
    arma::field <arma::vec> temp(3);
    temp(0) = {1,2,3};
    temp(1) = {4,5};
    temp(2) = {6,7,8,9};
    
}

// STATE SPACE SAMPLER: ------------------------------------------------------
arma::field<arma::field<arma::mat>> Omega_set(const arma::mat &G) {
    int N = G.n_cols; // dimension of adj matrix
    
    arma::field<arma::mat> c(N);
    arma::field<arma::mat> b(N, N);
    arma::field<arma::mat> a(N);
    
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            b(i, j) = arma::mat(1, 2, arma::fill::zeros);
        }
    }
    
    // a -------------------------------------------------------
    for(int i = 0; i < N; i++) {
        arma::mat a_i(1, 2, arma::fill::zeros);
        arma::uvec sec_elem = arma::find(G.row(i) == 1);
        
        for(int j = 0; j < sec_elem.n_elem; j++) {
            int sec_ind = sec_elem(j);
            arma::uvec third_elem = arma::find(G.row(sec_ind) == 1);
            
            for(int k = 0; k < third_elem.n_elem; k++) {
                int third_ind = third_elem(k);
                arma::mat temp(1,2);
                temp(0,0) = sec_ind; temp(0,1) = third_ind;
                a_i = arma::join_vert(a_i, temp+1);
            }
        }
        
        a_i = a_i.rows(1, a_i.n_rows-1);
        a(i) = a_i;
    }
    
    // b -------------------------------------------------------
    for(int i = 0; i < N; i++) {
        arma::uvec sec_elem = arma::find(G.row(i) == 1);
        
        for(int j = 0; j < sec_elem.n_elem; j++) {
            int sec_ind = sec_elem(j);
            arma::uvec third_elem = arma::find(G.row(sec_ind) == 1);
            
            for(int k = 0; k < third_elem.n_elem; k++) {
                int third_ind = third_elem(k);
                arma::uvec fourth_elem = arma::find(G.row(third_ind) == 1);
                
                for(int l = 0; l < fourth_elem.n_elem; l++) {
                    int fourth_ind = fourth_elem(l);
                    arma::mat temp(1,2);
                    temp(0,0) = sec_ind; temp(0,1) = third_ind;
                    b(i, fourth_ind) = arma::join_vert(b(i, fourth_ind), temp+1);
                }
            }
        }
    }
    
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            if (b(i,j).n_rows > 1) {
                b(i, j) = b(i, j).rows(1, b(i,j).n_rows - 1);
            } else{
                arma::mat temp(1,2); temp(0,0) = -1; temp(0,1) = -1;
                b(i,j) = temp;
            }
        }
    }
    
    // c -------------------------------------------------------
    for(int i = 0; i < N; i++) {
        arma::mat c_i(1, 2, arma::fill::zeros);
        arma::uvec sec_elem = arma::find(G.col(i) == 1);
        
        for(int j = 0; j < sec_elem.n_elem; j++) {
            int sec_ind = sec_elem(j);
            arma::uvec third_elem = arma::find(G.col(sec_ind) == 1);
            
            for(int k = 0; k < third_elem.n_elem; k++) {
                int third_ind = third_elem(k);
                arma::mat temp(1,2);
                temp(0,0) = third_ind; temp(0,1) = sec_ind;
                c_i = arma::join_vert(c_i, temp+1);
            }
        }
        
        c_i = c_i.rows(1, c_i.n_rows-1);
        c(i) = c_i;
    }
    
    arma::field<arma::field<arma::mat>> Omega_List(3);
    Omega_List(0) = c; Omega_List(1) = b; Omega_List(2) = a;
    
    return Omega_List;
}

const arma::field<arma::field<arma::mat>> Omega_List_GLOBAL = Omega_set(adj_mat);

arma::mat Omega_fun_cpp_new(const int k, const int n_i, const arma::vec &b_i) {
    
    arma::mat Omega_set;
    
    // b(k) is either 1, 2, 3 therefore subtract 1 for the index
    if (k == 1) {
        // () -> () -> 1-3
        Omega_set = Omega_List_GLOBAL(0)(b_i(2) - 1);
    } else if (k <= n_i - 2) {
        // 1-3 -> () -> () -> 1-3
        Omega_set = Omega_List_GLOBAL(1)(b_i(k - 2) - 1, b_i(k + 1) - 1);
    } else if (k == n_i - 1) {
        // 1-3 -> () -> ()
        Omega_set = Omega_List_GLOBAL(2)(b_i(n_i - 3) - 1);
    }
    
    return Omega_set;
    
}

double log_f_i_cpp(const int i, const int ii, const arma::vec &pars, 
                   const arma::field<arma::uvec> &par_index,
                   const arma::vec &y_1, arma::vec t_pts, const arma::vec &id, 
                   const arma::vec &B, const arma::vec &y_2, const int n_sub) {
    // par_index KEY: (0) zeta, (1) misclass, (2) delta, (3) tau2, (4) sigma2, (5) init
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    double in_value = 0;
    
    arma::vec eids = id;
    arma::uvec sub_ind = arma::find(eids == i);
    int n_i = sub_ind.max() - sub_ind.min() + 1;
    
    // Subsetting the data to relate only to this participant
    arma::mat b_i = B;
    arma::vec y_1_sub = y_1.elem(sub_ind);
    arma::mat y_2_sub = y_2.elem(sub_ind);
    
    arma::vec vec_zeta_content = pars.elem(par_index(0) - 1);
    arma::mat zeta = arma::reshape(vec_zeta_content, 4, 1); 
    
    arma::vec vec_init_content = pars.elem(par_index(5) - 1);
    arma::vec init_logit = {1, exp(vec_init_content(0)), exp(vec_init_content(1))};
    arma::vec P_init = init_logit / arma::accu(init_logit);
    
    arma::vec vec_misclass_content = pars.elem(par_index(1) - 1);
    arma::mat M = { {1, exp(vec_misclass_content(0)), exp(vec_misclass_content(1))},
                    {exp(vec_misclass_content(2)), 1, exp(vec_misclass_content(3))},
                    {exp(vec_misclass_content(4)), exp(vec_misclass_content(5)), 1}};
    arma::vec m_row_sums = arma::sum(M, 1);
    M = M.each_col() / m_row_sums;
    
    arma::vec delta = pars.elem(par_index(2) - 1);
    
    double log_tau2 = arma::as_scalar(pars.elem(par_index(3) - 1));
    double tau2 = exp(log_tau2);
    
    double log_sigma2 = arma::as_scalar(pars.elem(par_index(4) - 1));
    double sigma2 = exp(log_sigma2);
    
    // Full likelihood evaluation is not needed for updating pairs of b_i components
    for(int w=0; w < t_pts.n_elem; ++w){
        int k = t_pts(w);
        if(k==0){
            int b_k = b_i(k);
            int y_1_k = y_1_sub(k);
            double d_0 = D_2_calc(b_k, y_2_sub(k), tau2, sigma2, delta);
            in_value = in_value + log(P_init[b_k - 1]) + log(M(b_k - 1, y_1_k-1)) + log(d_0);
        } else{
            arma::colvec z_i = {1};
            
            double q1_sub = arma::as_scalar(zeta.row(0) * z_i);
            double q1 = exp(q1_sub);
            double q2_sub = arma::as_scalar(zeta.row(1) * z_i);
            double q2 = exp(q2_sub);
            double q3_sub = arma::as_scalar(zeta.row(2) * z_i);
            double q3 = exp(q3_sub);
            double q4_sub = arma::as_scalar(zeta.row(3) * z_i);
            double q4 = exp(q4_sub);
            
            arma::mat Q = { {  1,  q1,  0},
                            {  0,   1,  q2},
                            { q3,  q4,   1}};
            arma::vec q_row_sums = arma::sum(Q, 1);
            arma::mat P_i = Q.each_col() / q_row_sums;
            
            int b_k_1 = b_i(k-1);
            int b_k = b_i(k);
            int y_1_k = y_1_sub(k);
            
            double d_k = D_2_calc(b_k, y_2_sub(k), tau2, sigma2, delta);
            in_value = in_value + log(P_i( b_k_1 - 1, b_k - 1)) 
                        + log(M(b_k - 1, y_1_k-1)) + log(d_k);
        }
    }
    
    return in_value;
}

arma::vec update_b_i_cpp(const arma::vec &EIDs, const arma::vec &pars,
                          const arma::field<arma::uvec> &par_index,
                          const arma::vec &y_1, const arma::vec &id,
                          arma::vec &b_curr, const arma::vec &y_2) {
    
    // par_index KEY: (0) zeta, (1) misclass, (2) delta, (3) tau2, (4) sigma2
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    arma::vec B_return(b_curr.n_elem, arma::fill::zeros);
    
    omp_set_num_threads(10);
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(id == i);
        
        arma::vec b_i = b_curr.elem(sub_ind);
        
        int n_i = sub_ind.n_elem; 
        
        for (int k = 0; k < n_i - 1; k++) {
            
            arma::vec t_pts;
            if (k == n_i - 2) {
                t_pts = arma::linspace(k, k+1, 2);
            } else {
                t_pts = arma::linspace(k, k+2, 3);
            }
            
            arma::vec pr_B = b_i;
            
            // Sample and update the two neighboring states
            arma::mat Omega_set = Omega_fun_cpp_new(k + 1, n_i, b_i);
            
            int sampled_index = arma::randi(arma::distr_param(1, Omega_set.n_rows));
            
            pr_B.rows(k, k+1) = Omega_set.row(sampled_index-1).t();
            
            double log_target_prev = log_f_i_cpp(i, ii, pars, par_index, y_1, 
                                                 t_pts, id, b_i, y_2, 
                                                 EIDs.n_elem);
            
            double log_target = log_f_i_cpp(i, ii, pars, par_index, y_1, 
                                            t_pts, id, pr_B, y_2, 
                                            EIDs.n_elem);
            
            // Note that the proposal probs cancel in the MH ratio
            double diff_check = log_target - log_target_prev;
            double min_log = log(arma::randu(arma::distr_param(0,1)));
            if(diff_check > min_log){
                b_i = pr_B;
            }
        }
        B_return.elem(sub_ind) = b_i;
    }
    
    return B_return;
}

// [[Rcpp::export]]
arma::mat state_space_sampler(const int steps, const int burnin, 
                              const arma::vec &EIDs, const arma::vec &pars,  
                              const arma::field<arma::uvec> &par_index,
                              const arma::vec &y_1, const arma::vec &y_2,
                              const arma::vec &id, const arma::vec &t) {
    
    
    arma::mat B_master(steps - burnin, y_1.n_elem, arma::fill::zeros);
    arma::vec prev_B = y_1;
    
    for(int ttt = 0; ttt < steps; ttt++) {
        
        Rcpp::Rcout << ttt << std::endl;
        arma::vec curr_B = update_b_i_cpp(EIDs, pars, par_index, y_1, id, prev_B, y_2);
        if(ttt >= burnin)  B_master.row(ttt - burnin) = curr_B.t();
        prev_B = curr_B;
        
    }
    return B_master;
}



