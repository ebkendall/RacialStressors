#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

using namespace Rcpp;

// The adjacency matrix defining the possible state transitions --------------
const arma::mat adj_mat = { {1, 1, 1},
                            {1, 1, 1},
                            {1, 1, 1}};

//  FUNCTIONS: ---------------------------------------------------------------
double D_2_calc_fix(const int state_num, const arma::vec y_2_i_state, 
                    const double tau2, const arma::vec sigma_2_vec, 
                    const arma::vec delta,const arma::vec x_i, 
                    const arma::vec gamma, arma::vec delta_new, double dler_i) {
    // par_index KEY: (0) zeta, (1) misclass, (2) delta (mu, alpha, beta), 
    //                (3) tau2, (4) sigma2, (5) gamma, (6) delta_new

    double d_2_val = 0;
    int n_i = y_2_i_state.n_elem;
    double mu = delta(0);
    double alpha = delta(1);
    double beta = delta(2);
    
    if(state_num == 1) {

        double scalar_mean = arma::dot(x_i, gamma) + (dler_i * delta_new(0)) 
                                + mu;

        arma::vec mean_ii = arma::vec(n_i, arma::fill::ones);
        mean_ii = scalar_mean * mean_ii;
        
        arma::mat var_ii = arma::mat(n_i, n_i, arma::fill::ones);
        var_ii = sigma_2_vec(0) * var_ii;
        var_ii.diag() += tau2;

        // Checking whether variance is invertible based on the values given
        arma::mat C;
        bool success = arma::inv(C,var_ii);
        if(!success) {
            Rcpp::Rcout << "bad variance terms" << std::endl;
            d_2_val = -1 * arma::datum::inf;
        } else {
            arma::vec log_d_2 = dmvnorm(y_2_i_state.t(), mean_ii, var_ii, true);
            d_2_val = arma::as_scalar(log_d_2);
        }
        
    } else if(state_num == 2) {

        double scalar_mean = arma::dot(x_i, gamma) + (dler_i * delta_new(1)) 
                                + mu + alpha;
        
        arma::vec mean_ii = arma::vec(n_i, arma::fill::ones);
        mean_ii = scalar_mean * mean_ii;
        
        arma::mat var_ii = arma::mat(n_i, n_i, arma::fill::ones);
        var_ii = sigma_2_vec(1) * var_ii;
        var_ii.diag() += tau2;
        
        // Checking whether variance is invertible based on the values given
        arma::mat C;
        bool success = arma::inv(C,var_ii);
        if(!success) {
            Rcpp::Rcout << "bad variance terms" << std::endl;
            d_2_val = -1 * arma::datum::inf;
        } else {
            arma::vec log_d_2 = dmvnorm(y_2_i_state.t(), mean_ii, var_ii, true);
            d_2_val = arma::as_scalar(log_d_2);
        }
    } else {

        double scalar_mean = arma::dot(x_i, gamma) + (dler_i * delta_new(2)) 
                                + mu + beta;
        
        arma::vec mean_ii = arma::vec(n_i, arma::fill::ones);
        mean_ii = scalar_mean * mean_ii;
        
        arma::mat var_ii = arma::mat(n_i, n_i, arma::fill::ones);
        var_ii = sigma_2_vec(2) * var_ii;
        var_ii.diag() += tau2;
        
        // Checking whether variance is invertible based on the values given
        arma::mat C;
        bool success = arma::inv(C,var_ii);
        if(!success) {
            Rcpp::Rcout << "bad variance terms" << std::endl;
            d_2_val = -1 * arma::datum::inf;
        } else {
            arma::vec log_d_2 = dmvnorm(y_2_i_state.t(), mean_ii, var_ii, true);
            d_2_val = arma::as_scalar(log_d_2);
        }
    }
    
    return d_2_val;
}

// [[Rcpp::export]]
double fn_log_post_continuous(const arma::vec &EIDs, const arma::vec &pars,  
                              const arma::field<arma::vec> &prior_par, 
                              const arma::field<arma::uvec> &par_index,
                              const arma::vec &y_1, const arma::vec &id, 
                              const arma::vec &y_2, const arma::mat &cov_info,
                              arma::field<arma::vec> B, const int covariate_struct) {
    
    // par_index KEY: (0) zeta, (1) misclass, (2) delta (mu, alpha, beta), 
    //                (3) tau2, (4) sigma2, (5) gamma, (6) delta_new
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    arma::vec in_vals(EIDs.n_elem, arma::fill::zeros);
    
    // Initial state probabilities
    arma::vec init = {1, 0, 0};
    
    // Populate the parameters dependent on the covariate structure: zeta & gamma
    arma::vec vec_zeta_content = pars.elem(par_index(0) - 1);
    arma::mat zeta;
    arma::vec gamma = pars.elem(par_index(5) - 1);      // covariate effects
    arma::vec delta_new = pars.elem(par_index(6) - 1);  // state DLER effect
    
    if(covariate_struct == 1) {
        zeta = arma::reshape(vec_zeta_content, 6, 4); 
    } else if(covariate_struct == 2) {
        zeta = arma::reshape(vec_zeta_content, 6, 2); 
    } else {
        zeta = arma::reshape(vec_zeta_content, 6, 5); 
    }
    
    arma::vec delta = pars.elem(par_index(2) - 1); // c(mu, alpha, beta)
    
    double log_tau2 = arma::as_scalar(pars.elem(par_index(3) - 1));
    double tau2 = exp(log_tau2);
    
    arma::vec log_sigma2 = pars.elem(par_index(4) - 1);
    arma::vec sigma_2_vec = {exp(log_sigma2(0)), exp(log_sigma2(1)), exp(log_sigma2(2))};
    
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        int i = EIDs(ii);
        
        arma::vec b_i = B(ii);
        arma::uvec s1 = arma::find(b_i == 1);
        arma::uvec s2 = arma::find(b_i == 2);
        arma::uvec s3 = arma::find(b_i == 3);
        
        // Sub-setting the data
        arma::uvec sub_ind = arma::find(id == i);
        arma::vec y_1_i = y_1.elem(sub_ind);
        arma::vec y_2_i = y_2.elem(sub_ind);
        
        // Structuring the covariates for ** age, sex, pEdu **
        arma::mat cov_info_i = cov_info.rows(sub_ind);
        arma::colvec x_i = {cov_info_i(0,0), cov_info_i(0,1), cov_info_i(0,2)};
        arma::colvec z_i;

        if(covariate_struct == 1) {
            // age, sex, pEdu
            z_i = {1, cov_info_i(0,0), cov_info_i(0,1), cov_info_i(0,2)};
        } else if(covariate_struct == 2) {
            // DLER
            z_i = {1, cov_info_i(0,3)};
        } else {
            // age, sex, pEdu, DLER
            z_i = {1, cov_info_i(0,0), cov_info_i(0,1), 
                    cov_info_i(0,2), cov_info_i(0,3)};
        }

        // DLER information
        double dler_i = cov_info_i(0,3);
        
        // Evaluating the probability transition matrix
        double q1_sub = arma::as_scalar(zeta.row(0) * z_i);
        double q1 = exp(q1_sub);
        double q2_sub = arma::as_scalar(zeta.row(1) * z_i);
        double q2 = exp(q2_sub);
        double q3_sub = arma::as_scalar(zeta.row(2) * z_i);
        double q3 = exp(q3_sub);
        double q4_sub = arma::as_scalar(zeta.row(3) * z_i);
        double q4 = exp(q4_sub);
        double q5_sub = arma::as_scalar(zeta.row(4) * z_i);
        double q5 = exp(q5_sub);
        double q6_sub = arma::as_scalar(zeta.row(5) * z_i);
        double q6 = exp(q6_sub);
        
        arma::mat Q = { {  1,  q1,  q2},
                        { q3,   1,  q4},
                        { q5,  q6,   1}};
        arma::vec q_row_sums = arma::sum(Q, 1);
        arma::mat P = Q.each_col() / q_row_sums;
        
        // Likelihood contribution from state space ---------------------------
        // The initial probability for state 1 is always 1
        double log_state_val = 0;
        
        for(int k = 1; k < y_2_i.n_elem; k++) {
            if(y_1_i(k) > 1) {
                log_state_val = log_state_val + log(P(b_i(k-1) - 1, b_i(k) - 1));
            }
        }

        // Likelihood contribution from y_2 -----------------------------------
        double log_y_val = 0;
        
        if(s1.n_elem > 0) {
            int state_num = 1;
            arma::vec y_2_i_state = y_2_i.elem(s1);

            double d_num = D_2_calc_fix(state_num, y_2_i_state, tau2, 
                                        sigma_2_vec, delta, x_i, gamma,
                                        delta_new, dler_i);

            // Numerical stability check
            arma::vec d_check = {d_num};
            if(d_check.has_inf()) { 
                log_y_val = -1 * arma::datum::inf;
            } else {
                log_y_val = log_y_val + d_num;
            }
        }
        if(s2.n_elem > 0) {
            int state_num = 2;
            arma::vec y_2_i_state = y_2_i.elem(s2);

            double d_num = D_2_calc_fix(state_num, y_2_i_state, tau2, 
                                        sigma_2_vec, delta, x_i, gamma,
                                        delta_new, dler_i);

            // Numerical stability check
            arma::vec d_check = {d_num};
            if(d_check.has_inf()) { 
                log_y_val = -1 * arma::datum::inf;
            } else {
                log_y_val = log_y_val + d_num;
            }            
        }
        if(s3.n_elem > 0) {
            int state_num = 3;
            arma::vec y_2_i_state = y_2_i.elem(s3);

            double d_num = D_2_calc_fix(state_num, y_2_i_state, tau2, 
                                        sigma_2_vec, delta, x_i, gamma,
                                        delta_new, dler_i);

            // Numerical stability check
            arma::vec d_check = {d_num};
            if(d_check.has_inf()) { 
                log_y_val = -1 * arma::datum::inf;
            } else {
                log_y_val = log_y_val + d_num;
            }
        }

        // Check for numerical instability
        arma::vec log_y_check = {log_y_val};
        if(log_y_check.has_inf()) { 
            double inf_return = -1 * arma::datum::inf;
            in_vals(ii) = inf_return;
            break;
        }
        if(log_y_check.has_nan()) { 
            double inf_return = -1 * arma::datum::inf;
            in_vals(ii) = inf_return;
            break;
        }
        
        in_vals(ii) = log_state_val + log_y_val;
    }

    if(in_vals.has_inf()) { 
        Rcpp::Rcout << "in_vals has Inf" << std::endl;
        double inf_return = -1 * arma::datum::inf;
        return inf_return;
    }
    if(in_vals.has_nan()) { 
        Rcpp::Rcout << "in_vals has NaN" << std::endl;
        double inf_return = -1 * arma::datum::inf;
        return inf_return;
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

// STATE SPACE SAMPLER (no labels): -------------------------------------------
double log_f_i_cpp_no_label(const int i, const int ii, const arma::vec &pars, 
                            const arma::field<arma::uvec> &par_index,
                            arma::vec t_pts, const arma::vec &id, 
                            const arma::vec &B, const arma::vec &y_2, 
                            const int n_sub, const arma::mat &cov_info,
                            const arma::vec &y_1, const int covariate_struct) {
    // par_index KEY: (0) zeta, (1) misclass, (2) delta (mu, alpha, beta), 
    //                (3) tau2, (4) sigma2, (5) gamma, (6) delta_new
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    double in_value = 0;
    
    arma::vec eids = id;
    arma::uvec sub_ind = arma::find(eids == i);
    
    // Subsetting the data to relate only to this participant
    arma::mat b_i = B;
    arma::uvec s1 = arma::find(b_i == 1);
    arma::uvec s2 = arma::find(b_i == 2);
    arma::uvec s3 = arma::find(b_i == 3);
    
    arma::mat y_2_sub = y_2.elem(sub_ind);
    arma::vec y_1_sub = y_1.elem(sub_ind);
    
    // Populate the parameters dependent on the covariate structure: zeta & gamma
    arma::vec vec_zeta_content = pars.elem(par_index(0) - 1);
    arma::mat zeta;
    arma::vec gamma = pars.elem(par_index(5) - 1);      // covariate effects
    arma::vec delta_new = pars.elem(par_index(6) - 1);  // state DLER effect
    
    if(covariate_struct == 1) {
        zeta = arma::reshape(vec_zeta_content, 6, 4); 
    } else if(covariate_struct == 2) {
        zeta = arma::reshape(vec_zeta_content, 6, 2); 
    } else {
        zeta = arma::reshape(vec_zeta_content, 6, 5); 
    }

    // Structuring the covariates for ** age, sex, pEdu **
    arma::mat cov_info_i = cov_info.rows(sub_ind);
    arma::colvec x_i = {cov_info_i(0,0), cov_info_i(0,1), cov_info_i(0,2)};
    arma::colvec z_i;

    if(covariate_struct == 1) {
        // age, sex, pEdu
        z_i = {1, cov_info_i(0,0), cov_info_i(0,1), cov_info_i(0,2)};
    } else if(covariate_struct == 2) {
        // DLER
        z_i = {1, cov_info_i(0,3)};
    } else {
        // age, sex, pEdu, DLER
        z_i = {1, cov_info_i(0,0), cov_info_i(0,1), 
                cov_info_i(0,2), cov_info_i(0,3)};
    }

    // DLER information
    double dler_i = cov_info_i(0,3);   
    
    arma::vec P_init = {1, 0, 0};
    
    arma::vec delta = pars.elem(par_index(2) - 1); // c(mu, alpha, beta)
    
    double log_tau2 = arma::as_scalar(pars.elem(par_index(3) - 1));
    double tau2 = exp(log_tau2);

    arma::vec log_sigma2 = pars.elem(par_index(4) - 1);
    arma::vec sigma_2_vec = {exp(log_sigma2(0)), exp(log_sigma2(1)), exp(log_sigma2(2))};
    
    // Full likelihood evaluation is not needed for updating pairs of b_i components
    for(int w=0; w < t_pts.n_elem; ++w){
        int k = t_pts(w);
        if(k==0){
            // Currently NEVER have k==0 because initial state is set to 1
            int b_k = b_i(k);
            int y_1_k = y_1_sub(k);
            in_value = in_value + log(P_init[b_k - 1]);
        } else{
            // Evaluating the probability transition matrix
            double q1_sub = arma::as_scalar(zeta.row(0) * z_i);
            double q1 = exp(q1_sub);
            double q2_sub = arma::as_scalar(zeta.row(1) * z_i);
            double q2 = exp(q2_sub);
            double q3_sub = arma::as_scalar(zeta.row(2) * z_i);
            double q3 = exp(q3_sub);
            double q4_sub = arma::as_scalar(zeta.row(3) * z_i);
            double q4 = exp(q4_sub);
            double q5_sub = arma::as_scalar(zeta.row(4) * z_i);
            double q5 = exp(q5_sub);
            double q6_sub = arma::as_scalar(zeta.row(5) * z_i);
            double q6 = exp(q6_sub);
            
            arma::mat Q = { {  1,  q1,  q2},
                            { q3,   1,  q4},
                            { q5,  q6,   1}};
            arma::vec q_row_sums = arma::sum(Q, 1);
            arma::mat P_i = Q.each_col() / q_row_sums;
            
            int b_k_1 = b_i(k-1);
            int b_k = b_i(k);
            int y_1_k = y_1_sub(k);
            
            in_value = in_value + log(P_i( b_k_1 - 1, b_k - 1));
        }
    }
    
    double log_y_val = 0;
    
    if(s1.n_elem > 0) {
        int state_num = 1;
        arma::vec y_2_i_state = y_2_sub.elem(s1);
        log_y_val = log_y_val + D_2_calc_fix(state_num, y_2_i_state, tau2, 
                                             sigma_2_vec, delta, x_i, gamma,
                                             delta_new, dler_i);
    }
    if(s2.n_elem > 0) {
        int state_num = 2;
        arma::vec y_2_i_state = y_2_sub.elem(s2);
        log_y_val = log_y_val + D_2_calc_fix(state_num, y_2_i_state, tau2, 
                                             sigma_2_vec, delta, x_i, gamma,
                                             delta_new, dler_i);
    }
    if(s3.n_elem > 0) {
        int state_num = 3;
        arma::vec y_2_i_state = y_2_sub.elem(s3);
        log_y_val = log_y_val + D_2_calc_fix(state_num, y_2_i_state, tau2, 
                                             sigma_2_vec, delta, x_i, gamma,
                                             delta_new, dler_i);
    }
    
    in_value = in_value + log_y_val;
    
    return in_value;
}

// [[Rcpp::export]]
arma::field<arma::vec> update_b_i(const arma::vec &EIDs, const arma::vec &pars,
                                  const arma::field<arma::uvec> &par_index,
                                  const arma::vec &id,arma::field<arma::vec> B, 
                                  const arma::vec &y_2, const arma::vec &y_1,
                                  const arma::mat &cov_info, 
                                  const int covariate_struct){
    
    // par_index KEY: (0) zeta, (1) misclass, (2) delta (mu, alpha, beta), 
    //                (3) tau2, (4) sigma2, (5) gamma, (6) delta_new
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    arma::field<arma::vec> B_return(EIDs.n_elem);
    
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(id == i);
        
        arma::vec b_i = B(ii);
        
        int n_i = sub_ind.n_elem; 
        arma::vec y_1_sub = y_1.elem(sub_ind);

        bool start_run = false;
        
        for (int k = 0; k < n_i - 1; k++) {

            // state 1 observed without error
            if((y_1_sub(k) == 1) && (y_1_sub(k+1) != 1)) start_run = true;
            
            if(start_run) {
                arma::vec t_pts;
                if (k == n_i - 2) {
                    t_pts = arma::linspace(k, k+1, 2);
                } else {
                    t_pts = arma::linspace(k, k+2, 3);
                }
                
                arma::vec pr_B = b_i;
                
                // First time point always has state 1 as the first state
                if((y_1_sub(k) == 1) && (y_1_sub(k+1) != 1)) {
                    double sampled_index = arma::randi(arma::distr_param(1, 3));
                    arma::colvec os = {1, sampled_index};
                    pr_B.rows(k, k+1) = os;
                } else {
                    // Sample and update the two neighboring states
                    arma::mat Omega_set = Omega_fun_cpp_new(k + 1, n_i, b_i);
                    
                    int sampled_index = arma::randi(arma::distr_param(1, Omega_set.n_rows));
                    
                    pr_B.rows(k, k+1) = Omega_set.row(sampled_index-1).t();
                }
                
                double log_target_prev = log_f_i_cpp_no_label(i, ii, pars, par_index,
                                                              t_pts, id, b_i, y_2,
                                                              EIDs.n_elem, cov_info,
                                                              y_1, covariate_struct);
                
                double log_target = log_f_i_cpp_no_label(i, ii, pars, par_index,
                                                         t_pts, id, pr_B, y_2,
                                                         EIDs.n_elem, cov_info,
                                                         y_1, covariate_struct);
                
                // Note that the proposal probs cancel in the MH ratio
                double diff_check = log_target - log_target_prev;
                double min_log = log(arma::randu(arma::distr_param(0,1)));
                if(diff_check > min_log){
                    b_i = pr_B;
                }   
            }
        }
        
        B_return(ii) = b_i;
    }
    
    return B_return;
}

// [[Rcpp::export]]
arma::vec brute_force_ss(const int i, const int ii, const arma::vec &pars,
                         const arma::field<arma::uvec> &par_index,
                         const arma::vec &id, const arma::vec &y_2, 
                         const arma::vec &y_1, const arma::mat &cov_info, 
                         const int covariate_struct,
                         const arma::field<arma::mat> state_combos,
                         const arma::vec t_pts, int n_sub) {
    // par_index KEY: (0) zeta, (1) misclass, (2) delta (mu, alpha, beta), 
    //                (3) tau2, (4) sigma2, (5) gamma, (6) delta_new
    // "i" is the numeric EID value
    // "ii" is the index of the EID
    
    arma::vec eids = id;
    arma::uvec sub_ind = arma::find(eids == i);
    
    arma::vec y_1_sub = y_1.elem(sub_ind);
    arma::mat y_2_sub = y_2.elem(sub_ind);
    
    arma::uvec non_baseline = arma::find(y_1_sub != 1);
    
    arma::vec curr_b(sub_ind.n_elem, arma::fill::ones);
    arma::vec best_b(sub_ind.n_elem, arma::fill::ones);
    
    // Minimize the log-likelihood (i.e. maximize the negative log likelihood)
    double best_like = -1 * log_f_i_cpp_no_label(i, ii, pars, par_index, t_pts, id, 
                                            best_b, y_2, n_sub, cov_info, 
                                            y_1, covariate_struct);
    
    double start_like = best_like;
    
    Rcpp::Rcout << "num nonbaseline = " << non_baseline.n_elem << std::endl;
    
    // Minimum 7 depth
    for(int a = 1; a <= 3; a++) {
        Rcpp::Rcout << a << std::endl;
        curr_b(non_baseline(0)) = a;
        for(int b = 1; b <= 3; b++) {
            curr_b(non_baseline(1)) = b;
            for(int c = 1; c <= 3; c++) {
                curr_b(non_baseline(2)) = c;
                for(int d = 1; d <= 3; d++) {
                    curr_b(non_baseline(3)) = d;
                    for(int e = 1; e <= 3; e++) {
                        curr_b(non_baseline(4)) = e;
                        for(int f = 1; f <= 3; f++) {
                            curr_b(non_baseline(5)) = f;
                            for(int g = 1; g <= 3; g++) {
                                curr_b(non_baseline(6)) = g;
                                
                                if(non_baseline.n_elem == 8) {
                                    // Indices more = 1
                                    for(int h = 1; h <= 3; h++) {
                                        curr_b(non_baseline(7)) = h;
                                        
                                        // Evaluate the likelihood to test
                                        double curr_like = -1 * log_f_i_cpp_no_label(i, ii, pars, par_index, t_pts, id, 
                                                                                curr_b, y_2, n_sub, cov_info, 
                                                                                y_1, covariate_struct);
                                        if(curr_like < best_like) {
                                            best_b = curr_b;
                                            best_like = curr_like;
                                        }
                                    }
                                } else if(non_baseline.n_elem == 9){
                                    // Indices more = 2
                                    for(int h = 0; h < state_combos(1).n_rows; h++) {
                                        curr_b.subvec(non_baseline(7),non_baseline(8)) = state_combos(1).row(h).t();
                                        
                                        // Evaluate the likelihood to test
                                        double curr_like = -1 * log_f_i_cpp_no_label(i, ii, pars, par_index, t_pts, id, 
                                                                                curr_b, y_2, n_sub, cov_info, 
                                                                                y_1, covariate_struct);
                                        if(curr_like < best_like) {
                                            best_b = curr_b;
                                            best_like = curr_like;
                                        }
                                    }
                                    
                                } else if(non_baseline.n_elem == 10){
                                    // Indices more = 3
                                    for(int h = 0; h < state_combos(2).n_rows; h++) {
                                        curr_b.subvec(non_baseline(7),non_baseline(9)) = state_combos(2).row(h).t();
                                        
                                        // Evaluate the likelihood to test
                                        double curr_like = -1 * log_f_i_cpp_no_label(i, ii, pars, par_index, t_pts, id, 
                                                                                curr_b, y_2, n_sub, cov_info, 
                                                                                y_1, covariate_struct);
                                        if(curr_like < best_like) {
                                            best_b = curr_b;
                                            best_like = curr_like;
                                        }
                                    }
                                    
                                } else if(non_baseline.n_elem == 11){
                                    // Indices more = 4
                                    for(int h = 0; h < state_combos(3).n_rows; h++) {
                                        curr_b.subvec(non_baseline(7),non_baseline(10)) = state_combos(3).row(h).t();
                                        
                                        // Evaluate the likelihood to test
                                        double curr_like = -1 * log_f_i_cpp_no_label(i, ii, pars, par_index, t_pts, id, 
                                                                                curr_b, y_2, n_sub, cov_info, 
                                                                                y_1, covariate_struct);
                                        if(curr_like < best_like) {
                                            best_b = curr_b;
                                            best_like = curr_like;
                                        }
                                    }
                                    
                                } else if(non_baseline.n_elem == 12){
                                    // Indices more = 5
                                    for(int h = 0; h < state_combos(4).n_rows; h++) {
                                        curr_b.subvec(non_baseline(7),non_baseline(11)) = state_combos(4).row(h).t();
                                        
                                        // Evaluate the likelihood to test
                                        double curr_like = -1 * log_f_i_cpp_no_label(i, ii, pars, par_index, t_pts, id, 
                                                                                curr_b, y_2, n_sub, cov_info, 
                                                                                y_1, covariate_struct);
                                        if(curr_like < best_like) {
                                            best_b = curr_b;
                                            best_like = curr_like;
                                        }
                                    }
                                    
                                } else if(non_baseline.n_elem == 13){
                                    // Indices more = 6
                                    for(int h = 0; h < state_combos(5).n_rows; h++) {
                                        curr_b.subvec(non_baseline(7),non_baseline(12)) = state_combos(5).row(h).t();
                                        
                                        // Evaluate the likelihood to test
                                        double curr_like = -1 * log_f_i_cpp_no_label(i, ii, pars, par_index, t_pts, id, 
                                                                                curr_b, y_2, n_sub, cov_info, 
                                                                                y_1, covariate_struct);
                                        if(curr_like < best_like) {
                                            best_b = curr_b;
                                            best_like = curr_like;
                                        }
                                    }
                                    
                                } else if(non_baseline.n_elem == 14){
                                    // Indices more = 7
                                    for(int h = 0; h < state_combos(6).n_rows; h++) {
                                        curr_b.subvec(non_baseline(7),non_baseline(13)) = state_combos(6).row(h).t();
                                        
                                        // Evaluate the likelihood to test
                                        double curr_like = -1 * log_f_i_cpp_no_label(i, ii, pars, par_index, t_pts, id, 
                                                                                curr_b, y_2, n_sub, cov_info, 
                                                                                y_1, covariate_struct);
                                        if(curr_like < best_like) {
                                            best_b = curr_b;
                                            best_like = curr_like;
                                        }
                                    }
                                    
                                } else {
                                    // Evaluate the likelihood to test
                                    double curr_like = -1 * log_f_i_cpp_no_label(i, ii, pars, par_index, t_pts, id, 
                                                                            curr_b, y_2, n_sub, cov_info, 
                                                                            y_1, covariate_struct);
                                    if(curr_like < best_like) {
                                        best_b = curr_b;
                                        best_like = curr_like;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    Rcpp::Rcout << "start likelihood = " << start_like << std::endl;
    Rcpp::Rcout << "end likelihood = " << best_like << std::endl;
    return best_b;
    
}