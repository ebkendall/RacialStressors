#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

using namespace Rcpp;

const arma::mat adj_mat = { {1, 1, 1},
                            {1, 1, 1},
                            {1, 1, 1}};

//  FUNCTIONS: ---------------------------------------------------------------
double D_2_calc_fix(const int state, const double y_2_k, const double tau2, 
                    const arma::vec sigma_2_vec, const arma::vec delta,
                    const arma::vec x, const arma::vec gamma) {
    double d_2_val = 0;
    
    if(state == 1) {
        double mean_ii = arma::dot(x, gamma) + delta(0);
        double var_ii = tau2 + sigma_2_vec(0);
        double sd_ii = sqrt(var_ii);
        d_2_val = arma::normpdf(y_2_k, mean_ii, sd_ii);
    } else if(state == 2) {
        double mean_ii = arma::dot(x, gamma) + delta(0) + delta(1);
        double var_ii = tau2 + sigma_2_vec(1);
        double sd_ii = sqrt(var_ii);
        d_2_val = arma::normpdf(y_2_k, mean_ii, sd_ii);
    } else {
        double mean_ii = arma::dot(x, gamma) + delta(0) + delta(2);
        double var_ii = tau2 + sigma_2_vec(2);
        double sd_ii = sqrt(var_ii);
        d_2_val = arma::normpdf(y_2_k, mean_ii, sd_ii);
    }
    
    return d_2_val;
}

// [[Rcpp::export]]
double fn_log_post_continuous(const arma::vec &EIDs, const arma::vec &pars,  
                              const arma::field<arma::vec> &prior_par, 
                              const arma::field<arma::uvec> &par_index,
                              const arma::vec &y_1, const arma::vec &id, 
                              const arma::vec &y_2, const arma::mat &cov_info,
                              const bool case_b) {
    
    // par_index KEY: (0) zeta, (1) misclass, (2) delta, (3) tau2, (4) sigma2, 
    //                (5) gamma, (6) zeta_tilde, (7) sigma2_zeta
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    arma::vec in_vals(EIDs.n_elem, arma::fill::zeros);
    
    // Initial state probabilities
    arma::vec init = {1, 0, 0};
    
    // Populate the transition probability matrix (independent of time)
    arma::vec vec_zeta_content = pars.elem(par_index(0) - 1);
    arma::mat zeta = arma::reshape(vec_zeta_content, 6, 2); // ****************
    
    arma::vec delta = pars.elem(par_index(2) - 1);
    
    double log_tau2 = arma::as_scalar(pars.elem(par_index(3) - 1));
    double tau2 = exp(log_tau2);
    
    arma::vec log_sigma2 = pars.elem(par_index(4) - 1);
    arma::vec sigma_2_vec = {exp(log_sigma2(0)), exp(log_sigma2(1)), exp(log_sigma2(2))};
    
    arma::vec gamma = pars.elem(par_index(5) - 1);

    // Manually populate the misclassification probabilities
    arma::vec vec_misclass_content = pars.elem(par_index(1) - 1);
    arma::mat M = { {1, exp(vec_misclass_content(0)), exp(vec_misclass_content(1))},
                    {0, 1, exp(vec_misclass_content(2))},
                    {0, exp(vec_misclass_content(3)), 1}};
    arma::vec m_row_sums = arma::sum(M, 1);
    M = M.each_col() / m_row_sums;
    
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        int i = EIDs(ii);
        
        arma::rowvec val;
        double log_norm = 0;
        
        // Subsetting the data
        arma::uvec sub_ind = arma::find(id == i);
        arma::vec y_1_i = y_1.elem(sub_ind);
        arma::vec y_2_i = y_2.elem(sub_ind);
        
        // Evaluating the probability transition matrix
        arma::mat cov_info_i = cov_info.rows(sub_ind);
        // arma::colvec z_i = {1, cov_info_i(0,0), cov_info_i(0,1), cov_info_i(0,2), cov_info_i(0,3)};
        // arma::colvec x_i = {cov_info_i(0,0), cov_info_i(0,1), cov_info_i(0,2), cov_info_i(0,3)};
        arma::colvec z_i = {1, cov_info_i(0,0)}; // ****************************
        arma::colvec x_i = {cov_info_i(0,0)}; // *******************************
        
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
        
        // Likelihood component from y_2
        double d_1 = D_2_calc_fix(1, y_2_i(0), tau2, sigma_2_vec, delta, x_i, gamma);
        double d_2 = D_2_calc_fix(2, y_2_i(0), tau2, sigma_2_vec, delta, x_i, gamma);
        double d_3 = D_2_calc_fix(3, y_2_i(0), tau2, sigma_2_vec, delta, x_i, gamma);

        arma::vec d_fill = {d_1, d_2, d_3};
        if(d_fill.has_inf()) { 
            double inf_return = -1 * arma::datum::inf;
            in_vals(ii) = inf_return;
            continue;
        }
        if(d_fill.has_nan()) { 
            double inf_return = -1 * arma::datum::inf;
            in_vals(ii) = inf_return;
            continue;
        }
        arma::mat D_i_2 = arma::diagmat(d_fill);
        
        arma::mat init_transpose = init.t();
        
        arma::mat f_i = init_transpose * D_i_2;
        
        for(int k = 1; k < y_2_i.n_elem; k++) {
            
            // Likelihood component from y_2
            d_1 = D_2_calc_fix(1, y_2_i(k), tau2, sigma_2_vec, delta, x_i, gamma);
            d_2 = D_2_calc_fix(2, y_2_i(k), tau2, sigma_2_vec, delta, x_i, gamma);
            d_3 = D_2_calc_fix(3, y_2_i(k), tau2, sigma_2_vec, delta, x_i, gamma);

            d_fill = {d_1, d_2, d_3};
            if(d_fill.has_inf()) { 
                double inf_return = -1 * arma::datum::inf;
                in_vals(ii) = inf_return;
                break;
            }
            if(d_fill.has_nan()) { 
                double inf_return = -1 * arma::datum::inf;
                in_vals(ii) = inf_return;
                break;
            }
            
            D_i_2 = arma::diagmat(d_fill);
            
            arma::vec d_fill_1;
            if(case_b) {
                if(y_1_i(k) > 1) {
                    d_fill_1 = {1,1,1};
                } else {
                    d_fill_1 = {1,0,0};
                }
            } else {
                d_fill_1 = M.col(y_1_i(k) - 1);
            }

            arma::mat D_i_1 = arma::diagmat(d_fill_1);
            
            val = f_i * P * D_i_1 * D_i_2;
            
            arma::mat diag_val = arma::diagmat(val);
            arma::rowvec val_2 = val * diag_val;
            double norm_val = sqrt(arma::accu(val_2));
                        
            f_i = val / norm_val;
            log_norm = log_norm + log(norm_val);
        }

        if(in_vals.has_inf()) { 
            continue;
        } else if (in_vals.has_nan()) { 
            continue;
        } else {
            in_vals(ii) = log(arma::accu(f_i)) + log_norm;
        }
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
                            const bool case_b, const arma::vec &y_1) {
    // par_index KEY: (0) zeta, (1) misclass, (2) delta, (3) tau2, (4) sigma2, 
    //                (5) gamma, (6) zeta_tilde, (7) sigma2_zeta
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    double in_value = 0;
    
    arma::vec eids = id;
    arma::uvec sub_ind = arma::find(eids == i);
    
    // Subsetting the data to relate only to this participant
    arma::mat b_i = B;
    arma::mat y_2_sub = y_2.elem(sub_ind);
    arma::vec y_1_sub = y_1.elem(sub_ind);

    arma::mat cov_info_i = cov_info.rows(sub_ind);

    arma::vec vec_zeta_content = pars.elem(par_index(0) - 1);
    arma::mat zeta = arma::reshape(vec_zeta_content, 6, 5); 
    
    arma::colvec z_i = {1, cov_info_i(0, 0), cov_info_i(0, 1), cov_info_i(0, 2), cov_info_i(0, 3)};
    arma::colvec x_i = {cov_info_i(0, 0), cov_info_i(0, 1), cov_info_i(0, 2), cov_info_i(0, 3)};
    
    arma::vec P_init = {1, 0, 0};
    
    arma::vec delta = pars.elem(par_index(2) - 1);
    
    double log_tau2 = arma::as_scalar(pars.elem(par_index(3) - 1));
    double tau2 = exp(log_tau2);

    arma::vec log_sigma2 = pars.elem(par_index(4) - 1);
    arma::vec sigma_2_vec = {exp(log_sigma2(0)), exp(log_sigma2(1)), exp(log_sigma2(2))};
    
    arma::vec gamma = pars.elem(par_index(5) - 1);
    
    // Manually populate the misclassification probabilities
    arma::vec vec_misclass_content = pars.elem(par_index(1) - 1);
    arma::mat M = { {1, exp(vec_misclass_content(0)), exp(vec_misclass_content(1))},
                    {0, 1, exp(vec_misclass_content(2))},
                    {0, exp(vec_misclass_content(3)), 1}};
    arma::vec m_row_sums = arma::sum(M, 1);
    M = M.each_col() / m_row_sums;
    
    // Full likelihood evaluation is not needed for updating pairs of b_i components
    for(int w=0; w < t_pts.n_elem; ++w){
        int k = t_pts(w);
        if(k==0){
            // Currently NEVER have k==0 because initial state is set to 1
            int b_k = b_i(k);
            int y_1_k = y_1_sub(k);
            double d_0 = D_2_calc_fix(b_k, y_2_sub(k), tau2, sigma_2_vec, delta, x_i, gamma);
            in_value = in_value + log(P_init[b_k - 1]) + log(d_0);
            if(!case_b) {
                in_value = in_value + log(M(b_k - 1, y_1_k-1));
            }
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

            double d_k = D_2_calc_fix(b_k, y_2_sub(k), tau2, sigma_2_vec, delta, x_i, gamma);
            in_value = in_value + log(P_i( b_k_1 - 1, b_k - 1)) + log(d_k);
            if(!case_b) {
                in_value = in_value + log(M(b_k - 1, y_1_k-1));
            }
        }
    }
    
    return in_value;
}

arma::vec update_b_i_cpp_no_label( const arma::vec &EIDs, const arma::vec &pars,
                                   const arma::field<arma::uvec> &par_index,
                                   const arma::vec &id, arma::vec &b_curr, 
                                   const arma::vec &y_2, const arma::vec &y_1,
                                   const arma::mat &cov_info, const bool case_b) {
    
    // par_index KEY: (0) zeta, (1) misclass, (2) delta, (3) tau2, (4) sigma2, 
    //                (5) gamma, (6) zeta_tilde, (7) sigma2_zeta
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    arma::vec B_return(b_curr.n_elem, arma::fill::zeros);
    
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(id == i);
        
        arma::vec b_i = b_curr.elem(sub_ind);
        
        int n_i = sub_ind.n_elem; 
        arma::vec y_1_sub = y_1.elem(sub_ind);

        bool start_run = false;
        
        // Initial state is always 1 (so can start k == 1, not k == 0)
        for (int k = 0; k < n_i - 1; k++) {
            // keeping state 1 fixed without error
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
                                                              case_b, y_1);

                double log_target = log_f_i_cpp_no_label(i, ii, pars, par_index,
                                                         t_pts, id, pr_B, y_2,
                                                         EIDs.n_elem, cov_info,
                                                         case_b, y_1);

                // Note that the proposal probs cancel in the MH ratio
                double diff_check = log_target - log_target_prev;
                double min_log = log(arma::randu(arma::distr_param(0,1)));
                if(diff_check > min_log){
                    b_i = pr_B;
                }   
            }
        }
        B_return.elem(sub_ind) = b_i;
    }
    
    return B_return;
}

// [[Rcpp::export]]
arma::mat state_space_sampler_no_label(const int steps, const int burnin, 
                              const arma::vec &EIDs, const arma::mat &pars,  
                              const arma::field<arma::uvec> &par_index,
                              const arma::vec &y_2, const arma::vec &id, 
                              const arma::vec &t, const arma::vec &y_1,
                              const arma::mat &cov_info, const bool case_b) {
    // par_index KEY: (0) zeta, (1) misclass, (2) delta, (3) tau2, (4) sigma2, 
    //                (5) gamma, (6) zeta_tilde, (7) sigma2_zeta
    
    arma::mat B_master(steps - burnin, y_2.n_elem, arma::fill::zeros);
    arma::vec prev_B(y_2.n_elem, arma::fill::ones);
    int max_ind = pars.n_rows;
    
    for(int ttt = 0; ttt < steps; ttt++) {
        
        Rcpp::Rcout << "---> " << ttt << std::endl;
        int sampled_index = arma::randi(arma::distr_param(1, max_ind));
        arma::vec pars_ttt = pars.row(sampled_index - 1).t();
        
        arma::vec curr_B = update_b_i_cpp_no_label(EIDs, pars_ttt, par_index, id, 
                                                   prev_B, y_2, y_1, cov_info,
                                                   case_b);
        
        if(ttt >= burnin)  B_master.row(ttt - burnin) = curr_B.t();
        prev_B = curr_B;
        
    }
    return B_master;
}

// [[Rcpp::export]]
arma::field<arma::vec> viterbi_alg(const arma::vec &EIDs, const arma::vec &pars,
                                   const arma::field<arma::uvec> &par_index,
                                   const arma::vec &id, const arma::vec &y_2, 
                                   const arma::vec &y_1,const arma::mat &cov_info, 
                                   const bool case_b) {

    arma::vec state_vec = {1,2,3};
    
    arma::field<arma::vec> most_likely_ss(EIDs.n_elem);
    
    arma::vec vec_zeta_content = pars.elem(par_index(0) - 1);
    arma::mat zeta = arma::reshape(vec_zeta_content, 6, 5); 
    
    arma::vec vec_misclass_content = pars.elem(par_index(1) - 1);
    arma::mat M = { {1, exp(vec_misclass_content(0)), exp(vec_misclass_content(1))},
                    {0, 1, exp(vec_misclass_content(2))},
                    {0, exp(vec_misclass_content(3)), 1}};
    arma::vec m_row_sums = arma::sum(M, 1);
    M = M.each_col() / m_row_sums;
    
    arma::vec delta = pars.elem(par_index(2) - 1);
    
    double log_tau2 = arma::as_scalar(pars.elem(par_index(3) - 1));
    double tau2 = exp(log_tau2);
    
    arma::vec log_sigma2 = pars.elem(par_index(4) - 1);
    arma::vec sigma_2_vec = {exp(log_sigma2(0)), exp(log_sigma2(1)), 
                             exp(log_sigma2(2))};
    
    arma::vec gamma = pars.elem(par_index(5) - 1);
    
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(id == i);
        int n_i = sub_ind.n_elem; 

        arma::mat T_1(3, n_i, arma::fill::zeros);
        arma::mat T_2(3, n_i, arma::fill::zeros);
        
        arma::vec y_1_i = y_1.elem(sub_ind);
        arma::vec y_2_i = y_2.elem(sub_ind);
        
        // Evaluating the probability transition matrix
        arma::mat cov_info_i = cov_info.rows(sub_ind);
        arma::colvec z_i = {1, cov_info_i(0,0), cov_info_i(0,1), 
                            cov_info_i(0,2), cov_info_i(0,3)};
        arma::colvec x_i = {cov_info_i(0,0), cov_info_i(0,1), 
                            cov_info_i(0,2), cov_info_i(0,3)};
        
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

        for (int j = 0; j < n_i; j++) {
            // State 1 is observed without error and so initial state prob
            // is pi = (1,0,0)
            if((y_1_i(j) == 1) && (y_1_i(j+1) != 1)) {
                T_1(0,j) = D_2_calc_fix(y_1_i(j), y_2_i(j), tau2, sigma_2_vec, 
                                      delta, x_i, gamma);
                if(!case_b) T_1(0,j) = T_1(0,j) * M(0, y_1_i(j)-1);
            }
            
            if(y_1_i(j) > 1) {
                for(int s = 0; s < 3; s++) {
                    int state_s = s+1;
                    double B_s_j = D_2_calc_fix(state_s, y_2_i(j), tau2, 
                                                sigma_2_vec, delta, x_i, gamma);
                    if(!case_b) B_s_j = B_s_j * M(s, y_1_i(j)-1);
                    
                    arma::uword curr_ind;
                    double curr_max;
                    
                    arma::vec T_temp(3);
                    
                    for(int k = 0; k < 3; k++) {
                        T_temp(k) = T_1(k, j-1) * P(k,s) * B_s_j;
                        if(k==0 || (T_temp(k) > curr_max)) {
                            curr_max = T_temp(k);
                            curr_ind = k;
                        }
                    }
                    T_1(s, j) = curr_max;
                    T_2(s, j) = curr_ind;
                }
            }
        }
        
        arma::vec z_t(n_i, arma::fill::zeros);
        arma::vec x_t(n_i, arma::fill::zeros);
        
        arma::uword ind_t = T_1.col(n_i - 1).index_max();
        z_t(n_i - 1) = ind_t;
        x_t(n_i - 1) = state_vec(ind_t);
        
        for(int j = n_i-1; j > 0; j--) {
            z_t(j-1) = T_2(z_t(j), j);
            x_t(j-1) = state_vec(z_t(j-1));
        }
        
        most_likely_ss(ii) = x_t;
    }
    
    return most_likely_ss;
}

// [[Rcpp::export]]
double test_functions(const arma::vec &pars, const arma::field<arma::vec> &prior_par, 
                    const arma::field<arma::uvec> &par_index) {
    
    // arma::mat x = {pars(0), pars(1), pars(2), pars(3), pars(4), pars(5),
    //                pars(6), pars(7), pars(8), pars(9), pars(10), pars(11),
    //                pars(12), pars(13), pars(14), pars(15), pars(16), pars(17),
    //                pars(18), pars(19), pars(20), pars(21), pars(22), pars(23),
    //                pars(24), pars(25), pars(26), pars(27), pars(28), pars(29),
    //                pars(30), pars(31), pars(34), pars(35),
    //                pars(36), pars(37), pars(38), pars(39)};
    // arma::mat x2 = arma::join_cols(pars.subvec(0, 31), pars.subvec(34, 39));
    // arma::mat x3 = pars;
    // Rcpp::Rcout << x << std::endl;
    // Rcpp::Rcout << x2 << std::endl;
    // Rcpp::Rcout << x3 << std::endl;
    
    // arma::mat A = {{2,0,0},{0,2,0},{0,0,2}};
    // arma::mat C;
    // bool test_suc = arma::inv(C, A);
    // arma::mat B = arma::inv(A);
    // Rcpp::Rcout << test_suc << std::endl;
    // Rcpp::Rcout << B << std::endl;
    // Rcpp::Rcout << C << std::endl;
    
    // for(int j = 0; j < 20; j++) {
    //     double ind_unif = arma::randu();
    //     Rcpp::Rcout << ind_unif << ", log " << log(ind_unif) << std::endl;
    // }
    
    // double inf_return = -1 * arma::datum::inf;
    // if(inf_return == -arma::datum::inf) {
    //     Rcpp::Rcout << "Inf Check1" << std::endl;
    // }
    
    // inf_return = arma::datum::inf;
    // if(inf_return == arma::datum::inf) {
    //     Rcpp::Rcout << "Inf Check2" << std::endl;
    // }
    
    // Rcpp::Rcout << 10 % 3 << std::endl;
    // Rcpp::Rcout << 9 % 3 << std::endl;
    
    // arma::vec t = {1.1192126, 1.1192068, 1.6840622, 0.2306889};
    // arma::vec t2 = {-0.2498609, 0.4671621, -0.2409424, -0.8548681};
    
    // arma::mat t3 = arma::join_horiz(t, t2);
    
    // Rcpp::Rcout << arma::cov(t, t2) << std::endl;
    // Rcpp::Rcout << arma::cov(t.t(), t2.t()) << std::endl;
    // Rcpp::Rcout << arma::cov(t3) << std::endl;
    // Rcpp::Rcout << arma::cov(t3.t()) << std::endl;
    
    // arma::vec a = {1,2,3,4};
    // arma::vec b = {1,2,3,4};
    
    // if(arma::accu(a == b) == a.n_elem) {
    //     Rcpp::Rcout << "equal" << std::endl;
    // }

    for(int i = 0; i < 10; i++) {
        Rcpp::Rcout << arma::randi(arma::distr_param(1, 3)) << std::endl;
    }
    
    return 0;
}
