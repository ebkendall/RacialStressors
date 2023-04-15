#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// Defining the Omega_List as a global variable when pre-compiling ----------
const arma::mat adj_mat = {{1, 1, 1},
                           {0, 1, 1},
                           {1, 1, 1}};


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

//  FUNCTIONS: ---------------------------------------------------------------

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
                   const arma::field<arma::vec> &prior_par, const arma::field<arma::uvec> &par_index,
                   const arma::vec &y_1, arma::vec t_pts, const arma::vec &id, 
                   const arma::vec &B, const arma::vec &y_2, arma::mat &V_i, const int n_sub) {
    // par_index KEY: (0) zeta, (1) misclass, (2) delta, (3) tau2, (4) sigma2, (5) delta_i
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    double in_value = 0;
    
    arma::vec eids = id;
    arma::uvec vec_zeta_ind = par_index(0);
    arma::uvec vec_misclass_ind = par_index(1);
    
    arma::vec vec_zeta_content = pars.elem(vec_zeta_ind - 1);
    arma::vec vec_misclass_content = pars.elem(vec_misclass_ind - 1);
    
    // Manually populate the matrix
    arma::mat zeta = arma::reshape(vec_zeta_content, 5, 1); // NO TIME COMPONENT
    arma::mat M = {{1, exp(vec_misclass_content(0)), exp(vec_misclass_content(1))},
                   {exp(vec_misclass_content(2)), 1, exp(vec_misclass_content(3))},
                   {exp(vec_misclass_content(4)), exp(vec_misclass_content(5)), 1}};
    arma::vec m_row_sums = arma::sum(M, 1);
    M = M.each_col() / m_row_sums;
    
    // The time-homogeneous probability transition matrix
    arma::uvec sub_ind = arma::find(eids == i);
    
    // Subsetting the data to relate only to this participant
    arma::mat b_i = B;
    arma::vec y_1_sub = y_1.elem(sub_ind);
    arma::mat y_2_sub = y_2.elem(sub_ind);

    // Full likelihood evaluation is not needed for updating pairs of b_i components
    int n_i = sub_ind.max() - sub_ind.min() + 1;
    if (any(t_pts == -1)) { t_pts = arma::linspace(0, n_i - 1, n_i - 1);}
    
    // Full likelihood evaluation is not needed for updating pairs of b_i components
    for(int w=0; w < t_pts.n_elem; ++w){
        
        int k = t_pts(w);
        if(k==0){
            int b_k = b_i(k);
            int y_1_k = y_1_sub(k);
            in_value = in_value + log(M(b_k - 1, y_1_k-1));
        } else{
            double k_scale = k / 100;         // scaling time by 10
            
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
            double q5_sub = arma::as_scalar(zeta.row(4) * z_i);
            double q5 = exp(q5_sub);
            
            arma::mat Q = { {  1,  q1,  q2},
            {  0,   1,  q3},
            { q4,  q5,   1}};
            
            arma::vec q_row_sums = arma::sum(Q, 1);
            arma::mat P_i = Q.each_col() / q_row_sums;
            
            int b_k_1 = b_i(k-1);
            int b_k = b_i(k);
            int y_1_k = y_1_sub(k);
            in_value = in_value + log(P_i( b_k_1 - 1, b_k - 1)) + log(M(b_k - 1, y_1_k-1));
        }
    }
    
    // Likelihood from the RSA response
    
    arma::vec delta_i_vec = pars.elem(par_index(5) - 1);
    arma::mat delta_i_big = arma::reshape(delta_i_vec, n_sub, 3);
    arma::vec delta_i = delta_i_big.row(ii).t();

    arma::vec mean_vec_y_2 = V_i * delta_i;
    arma::mat cov_y_2 = arma::eye(n_i, n_i);
    arma::vec vec_tau2 = pars.elem(par_index(3) - 1);
    double tau2 = vec_tau2(0);
    cov_y_2 = tau2 * cov_y_2;

    double log_dens = arma::as_scalar(dmvnorm(y_2_sub.t(), mean_vec_y_2, cov_y_2, true));

    in_value = in_value + log_dens;

    return in_value;
}

// [[Rcpp::export]]
double fn_log_post_continuous(const arma::vec &EIDs, const arma::vec &pars,  
                              const arma::field<arma::vec> &prior_par, const arma::field<arma::uvec> &par_index,
                              const arma::vec &y_1, const arma::vec &id, const arma::vec &y_2) {
    
    // par_index KEY: (0) zeta, (1) misclass, (2) delta, (3) tau2, (4) sigma2, (5) delta_i
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    arma::vec in_vals(EIDs.n_elem, arma::fill::zeros);
    
    // Initial state probabilities
    arma::vec init = {1, 0, 0};
    
    // Manually populate the misclassification probabilities
    arma::vec vec_misclass_content = pars.elem(par_index(1) - 1);
    arma::mat M = { {1, exp(vec_misclass_content(0)), exp(vec_misclass_content(1))},
                    {exp(vec_misclass_content(2)), 1, exp(vec_misclass_content(3))},
                    {exp(vec_misclass_content(4)), exp(vec_misclass_content(5)), 1}};
    arma::vec m_row_sums = arma::sum(M, 1);
    M = M.each_col() / m_row_sums;
    
    // Populate the Transition probability matrix (independent of time)
    arma::vec vec_zeta_content = pars.elem(par_index(0) - 1);
    arma::mat zeta = arma::reshape(vec_zeta_content, 5, 1); 
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
    double q5_sub = arma::as_scalar(zeta.row(4) * z_i);
    double q5 = exp(q5_sub);
    arma::mat Q = { {  1,  q1,  q2},
                    {  0,   1,  q3},
                    { q4,  q5,   1}};
    arma::vec q_row_sums = arma::sum(Q, 1);
    arma::mat P = Q.each_col() / q_row_sums;
    
    // Defining key terms
    arma::vec delta_i_vec = pars.elem(par_index(5) - 1);
    arma::mat delta_i_big = arma::reshape(delta_i_vec, EIDs.n_elem, 3);
    
    double tau2 = arma::as_scalar(pars.elem(par_index(3) - 1));
    double sigma2 = arma::as_scalar(pars.elem(par_index(4) - 1));
    
    double tau_sigma = sqrt(tau2 + sigma2);
    
    omp_set_num_threads(6);
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        int i = EIDs(ii);
        
        // Subsetting the data
        arma::uvec sub_ind = arma::find(id == i);
        arma::vec y_1_i = y_1.elem(sub_ind);
        arma::vec y_2_i = y_2.elem(sub_ind);
        arma::vec delta_i = delta_i_big.row(ii).t();
        
        // Likelihood component from y_1
        arma::vec misclass_fill = M.col(y_1_i(0) - 1);
        arma::mat D_i_1 = arma::diagmat(misclass_fill);
        
        // Likelihood component from y_2
        double mean_1 = delta_i(0);
        double mean_2 = delta_i(0) + delta_i(1);
        double mean_3 = delta_i(0) + delta_i(2);
        
        double d_1 = arma::normpdf(y_2_i(0), mean_1, tau_sigma);
        double d_2 = arma::normpdf(y_2_i(0), mean_2, tau_sigma);
        double d_3 = arma::normpdf(y_2_i(0), mean_3, tau_sigma);
        
        arma::vec d_fill = {d_1, d_2, d_3};
        arma::mat D_i_2 = arma::diagmat(d_fill);
        
        arma::mat init_transpose = init.t();
        
        arma::mat f_i = init_transpose * D_i_1 * D_i_2;
        // arma::mat f_i = init_transpose * D_i_2;
        
        for(int k = 1; k < y_1_i.n_elem; k++) {
            
            // Likelihood component from y_1
            arma::vec misclass_fill = M.col(y_1_i(k) - 1);
            arma::mat D_i_1 = arma::diagmat(misclass_fill);
            
            // Likelihood component from y_2
            double d_1 = arma::normpdf(y_2_i(k), mean_1, tau_sigma);
            double d_2 = arma::normpdf(y_2_i(k), mean_2, tau_sigma);
            double d_3 = arma::normpdf(y_2_i(k), mean_3, tau_sigma);
            
            arma::vec d_fill = {d_1, d_2, d_3};
            arma::mat D_i_2 = arma::diagmat(d_fill);
            
            arma::mat temp = f_i * P * D_i_1 * D_i_2;
            // arma::mat temp = f_i * P * D_i_2;
            f_i = temp;
        }
        in_vals(ii) = log(arma::accu(f_i));
    }
    
    double in_value = arma::accu(in_vals);
    
    // Likelihood components from the Metropolis priors
    arma::vec p_mean = prior_par(0);
    arma::mat p_sd = arma::diagmat(prior_par(1));

    arma::mat x = arma::join_cols(vec_zeta_content, vec_misclass_content);
    // arma::mat x = vec_zeta_content;
    double log_prior_dens = arma::as_scalar(dmvnorm(x.t(), p_mean, p_sd, true));
    
    in_value = in_value + log_prior_dens;
    
    return in_value;
}



// [[Rcpp::export]]
Rcpp::List update_b_i_cpp(const int t, const arma::vec &EIDs, const arma::vec &pars,
                          const arma::field<arma::vec> &prior_par, const arma::field<arma::uvec> &par_index,
                          const arma::vec &y_1, const arma::vec &id,
                          arma::field <arma::vec> &B, arma::field <arma::mat> &V_i,
                          const arma::vec &y_2) {

    // par_index KEY: (0) zeta, (1) misclass, (2) delta, (3) tau2, (4) sigma2, (5) delta_i
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    //  ALWAYS DOUBLE CHECK THESE INDICES. WE LOSE THE NAMES FEATURE

    arma::field<arma::vec> B_return(EIDs.n_elem);
    arma::field<arma::mat> V_return(EIDs.n_elem);

    omp_set_num_threads(t) ;
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(id == i);

        int n_i = sub_ind.n_elem; 
        arma::vec y_1_sub = y_1.elem(sub_ind);

        // Subsetting fields
        arma::vec B_temp = B(ii);
        arma::mat V_temp = V_i(ii);

        // The first state is always S1, therefore we start at 1 instead of 0
        for (int k = 1; k < n_i - 1; k++) {
            
            arma::vec t_pts;
            if (k == n_i - 2) {
                t_pts = arma::linspace(k, k+1, 2);
            } else {
                t_pts = arma::linspace(k, k+2, 3);
            }
            
            arma::vec pr_B = B_temp;
            arma::mat pr_V = V_temp;

            // Sample and update the two neighboring states
            arma::mat Omega_set = Omega_fun_cpp_new(k + 1, n_i, B_temp);
            
            int sampled_index = arma::randi(arma::distr_param(1, Omega_set.n_rows));

            pr_B.rows(k, k+1) = Omega_set.row(sampled_index-1).t();

            // Adding clinical review
            bool valid_prop = true;

            if(valid_prop) {
                double log_target_prev = log_f_i_cpp(i, ii, pars, prior_par,
                                                     par_index, y_1, t_pts, id,
                                                     B_temp, y_2, V_temp, EIDs.n_elem);

                arma::vec col1(pr_B.n_elem, arma::fill::ones);
                arma::vec col2(pr_B.n_elem, arma::fill::zeros);
                col2.elem(arma::find(pr_B == 2)).ones();
                arma::vec col3(pr_B.n_elem, arma::fill::zeros);
                col3.elem(arma::find(pr_B == 3)).ones();

                pr_V = arma::join_horiz(col1, arma::join_horiz(col2, col3));
                
                double log_target = log_f_i_cpp(i, ii, pars, prior_par,
                                                par_index, y_1, t_pts, id,
                                                pr_B, y_2, pr_V, EIDs.n_elem);

                // Note that the proposal probs cancel in the MH ratio
                double diff_check = log_target - log_target_prev;
                double min_log = log(arma::randu(arma::distr_param(0,1)));
                // Rcpp::Rcout << diff_check << "  " << min_log << std::endl;
                if(diff_check > min_log){
                    B_temp = pr_B;
                    V_temp = pr_V;
                }
            }
        }
        B_return(ii) = B_temp;
        V_return(ii) = V_temp;
    }

    List B_V = List::create(B_return, V_return);

    return B_V;
}

// [[Rcpp::export]]
arma::mat update_delta_i_cpp(const arma::vec &y_2, const arma::vec &pars,
                          const arma::field<arma::uvec> &par_index,
                          const arma::field<arma::mat> &V_i,
                          const arma::vec &EIDs, const arma::vec &id) {
    // par_index KEY: (0) zeta, (1) misclass, (2) delta, (3) tau2, (4) sigma2, (5) delta_i
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    arma::mat delta_i(EIDs.n_elem, 3);
    
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        int i = EIDs(ii);
        
        arma::mat V_small = V_i(ii);
        
        arma::vec vec_sigma2 = pars.elem(par_index(4) - 1);
        double sigma2 = vec_sigma2(0);
        
        arma::vec inv_diag_upsilon = {1/sigma2, 1/sigma2, 1/sigma2};
        arma::mat up_solve = arma::diagmat(inv_diag_upsilon);
        
        arma::uvec sub_ind = arma::find(id == i);
        arma::vec y_sub = y_2.elem(sub_ind); 
        
        arma::vec vec_tau2 = pars.elem(par_index(3) - 1);
        double tau2 = vec_tau2(0);
        
        arma::vec delta = pars.elem(par_index(2) - 1);
        
        arma::mat V_inv = (1/tau2) * (V_small.t() * V_small) + up_solve;
        arma::mat V = arma::inv(V_inv);
        
        arma::vec M = V * ((1/tau2) * (V_small.t() * y_sub) + up_solve * delta);
        
        arma::mat delta_i_small = rmvnorm(1, M, V);
        
        delta_i.row(ii) = delta_i_small;
        
    }
    
    return delta_i;
}



// [[Rcpp::export]]
void test_functions(const arma::vec &pars, const arma::field<arma::vec> &prior_par, 
                    const arma::field<arma::uvec> &par_index) {

    // int N = 3;
    // Rcpp::Rcout << "Case (c) Full" << std::endl;
    // for(int w=0; w < N; w++) {
    //   Rcpp::Rcout << "() -> () -> " << w+1 << std::endl;
    //   Rcpp::Rcout << Omega_List_GLOBAL(0)(w) << std::endl;
    // }
    
    // Rcpp::Rcout << "Case (b) Full" << std::endl;
    // for(int i = 0; i < N; i++) {
    //   for(int j = 0; j < N; j++) {
    //     Rcpp::Rcout << i+1 << "-->" << j+1 << std::endl;
    //     Rcpp::Rcout << Omega_List_GLOBAL(1)(i, j) << std::endl;
    //   }
    // }
    
    // Rcpp::Rcout << "Case (a) Full" << std::endl;
    // for(int w=0; w < N; w++) {
    //   Rcpp::Rcout << w + 1 << " -> () -> ()" << std::endl;
    //   Rcpp::Rcout << Omega_List_GLOBAL(2)(w) << std::endl;
    // }

    arma::mat M = {{1,0,0},
                   {0,2,0},
                   {0,0,3}};
    arma::mat W = arma::expmat(M);
    Rcpp::Rcout << W << std::endl;

    arma::mat temp = 1.3 * M;
    Rcpp::Rcout << temp << std::endl;

    W = arma::expmat(temp);
    Rcpp::Rcout << W << std::endl;

    // arma::vec temp = {1,2};
    // Rcpp::Rcout << temp + 1 << std::endl;
    // 
    // for(int i = 0; i < 30; i++) {
    //   int sampled_index = arma::randi(arma::distr_param(1, 5));
    //   Rcpp::Rcout << sampled_index << std::endl;
    // }

    // arma::vec t_pts = {1,2};
    // if(any(t_pts == -1)) {
    //     Rcpp::Rcout << "wrong" << std::endl;
    // } else{
    //     Rcpp::Rcout << "correct" << std::endl;
    // }
    
    // t_pts = {-1};
    // if(any(t_pts == -1)) {
    //     Rcpp::Rcout << "correct" << std::endl;
    // } else{
    //     Rcpp::Rcout << "wrong" << std::endl;
    // }
    // Multivariate Normal Check
    // arma::uvec vec_zeta_ind = par_index(0);
    // arma::uvec vec_misclass_ind = par_index(1);
    
    // arma::vec vec_zeta_content = pars.elem(vec_zeta_ind - 1);
    // arma::vec vec_misclass_content = pars.elem(vec_misclass_ind - 1);
    
    // arma::vec p_mean = prior_par(0);
    // arma::mat p_sd = arma::diagmat(prior_par(1));
    // arma::mat x = arma::join_cols(vec_zeta_content, vec_misclass_content);
    
    // arma::vec log_prior_dens = dmvnorm(x.t(), p_mean, p_sd, true);
    
    // Rcpp::Rcout << log_prior_dens << std::endl;
    
    // //  Sub setting check
    // arma::vec id = {1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,3};
    // arma::vec y_1 = 5 * id;
    // arma::uvec sub_ind = arma::find(id == 2);
    
    // int n_i = sub_ind.n_elem; // DOUBLE CHECK and compare to other method
    // arma::vec y_1_sub = y_1.elem(sub_ind);
    
    // Rcpp::Rcout << sub_ind << std::endl;
    // Rcpp::Rcout << y_1_sub << std::endl;
    
    // Rcpp::Rcout << "Attempt 1: " << n_i << std::endl;
    
    // int n_i_2 = sub_ind.max() - sub_ind.min() + 1;
    // Rcpp::Rcout << sub_ind.max() << " " << sub_ind.min() << std::endl;
    
    // arma::vec col1(id.n_elem, arma::fill::zeros);
    // col1.elem(arma::find(id == 1)).ones();
    // arma::vec col2(id.n_elem, arma::fill::zeros);
    // col2.elem(arma::find(id == 2)).ones();
    // arma::vec col3(id.n_elem, arma::fill::zeros);
    // col3.elem(arma::find(id == 3)).ones();
    
    // arma::mat pr_V = arma::join_horiz(col1, arma::join_horiz(col2, col3));
    
    // Rcpp::Rcout << pr_V << std::endl;
    
}