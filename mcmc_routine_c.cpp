#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// Defining the Omega_List as a global variable when pre-compiling ----------
const arma::mat adj_mat = {{1, 1, 0},
                           {1, 1, 1},
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
                   const arma::vec &B, arma::vec it_indices) {
    // par_index KEY: (0) beta, (1) misclass, (2) mu_tilde, (3) tau2, (4) upsilon, (5) mu_i
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    double in_value = 0;
    
    arma::vec eids = id;
    arma::uvec vec_beta_ind = par_index(0);
    arma::uvec vec_misclass_ind = par_index(1);
    
    arma::vec vec_beta_content = pars.elem(vec_beta_ind - 1);
    arma::vec vec_misclass_content = pars.elem(vec_misclass_ind - 1);
    
    // Manually populate the matrix
    arma::mat beta = arma::reshape(vec_beta_content, 5, 2);
    // arma::mat M = { {1, 0, 0},
    //                 {0, 1, exp(vec_misclass_content(0))},
    //                 {0, exp(vec_misclass_content(1)), 1}};
    // arma::vec m_row_sums = arma::sum(M, 1);
    // M = M.each_col() / m_row_sums; 
    
    // The time-homogeneous probability transition matrix
    arma::uvec sub_ind = arma::find(eids == i);
    
    // Subsetting the data to relate only to this participant
    arma::mat b_i = B;
    arma::vec t_pts_sub = t_pts.elem(sub_ind);
    arma::vec y_1_sub = y_1.elem(sub_ind);
    
    // Full likelihood evaluation is not needed for updating pairs of b_i components
    for(int w=0; w < it_indices.n_elem;++w){
        arma::colvec z_i = {1, t_pts_sub(it_indices(w))}; // using the current time point
        double q1_sub = arma::as_scalar(beta.row(0) * z_i);
        double q1 = exp(q1_sub);
        double q2_sub = arma::as_scalar(beta.row(1) * z_i);
        double q2 = exp(q2_sub);
        double q3_sub = arma::as_scalar(beta.row(2) * z_i);
        double q3 = exp(q3_sub);
        double q4_sub = arma::as_scalar(beta.row(3) * z_i);
        double q4 = exp(q4_sub);
        double q5_sub = arma::as_scalar(beta.row(4) * z_i);
        double q5 = exp(q5_sub);
        
        arma::mat Q = { {  1,  q1,   0},
                        { q2,   1,  q3},
                        { q4,  q5,   1}};
        
        arma::vec q_row_sums = arma::sum(Q, 1);
        arma::mat P_i = Q.each_col() / q_row_sums;
        // Q = (t_pts_sub(it_indices(w)) - t_pts_sub(it_indices(w) - 1)) * Q;
        // arma::mat P_i = arma::expmat(Q);
        
        int b_k_1 = b_i(it_indices(w)-1,0);
        int b_k = b_i(it_indices(w), 0);
        int y_1_k = y_1_sub(it_indices(w));
        in_value = in_value + log(P_i( b_k_1 - 1, b_k - 1)); // + log(M(b_k - 1, y_1_k-1));
    }
    
    // Likelihood components from the other parts
    arma::vec p_mean = prior_par(0);
    arma::mat p_sd = arma::diagmat(prior_par(1));
    arma::mat x = vec_beta_content;
    // arma::mat x = arma::join_cols(vec_beta_content, vec_misclass_content);
    double log_prior_dens = arma::as_scalar(dmvnorm(x.t(), p_mean, p_sd, true));

    in_value = in_value + log_prior_dens;
    
    return in_value;
}

// [[Rcpp::export]]
double log_f_i_cpp_total(const arma::vec &EIDs, const arma::vec &pars,  
                         const arma::field<arma::vec> &prior_par, const arma::field<arma::uvec> &par_index,
                         const arma::vec &y_1, arma::vec t_pts, const arma::vec &id,
                         const arma::field <arma::vec> &B) {
    
    // par_index KEY: (0) beta, (1) misclass, (2) mu_tilde, (3) tau2, (4) upsilon, (5) mu_i
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    arma::vec in_vals(EIDs.n_elem, arma::fill::zeros);
    
    omp_set_num_threads(16);
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(id == i);
        int n_i = sub_ind.max() - sub_ind.min();
        
        // Starting the it_indices at 1 instead of 0 because there is 
        //      no contribution to the likelihood at the initial state
        arma::vec it_indices = arma::linspace(1, n_i, n_i); 
        in_vals(ii) = log_f_i_cpp(i, ii, pars, prior_par, par_index, y_1, t_pts, id, B(ii), it_indices);
    }
    
    double in_value = arma::accu(in_vals);
    
    return in_value;
}

// [[Rcpp::export]]
Rcpp::List update_b_i_cpp(const int t, const arma::vec &EIDs, const arma::vec &pars,
                          const arma::field<arma::vec> &prior_par, const arma::field<arma::uvec> &par_index,
                          const arma::vec &y_1, arma::vec t_pts, const arma::vec &id,
                          arma::field <arma::vec> &B, arma::field <arma::mat> &V_i) {

    // par_index KEY: (0) beta, (1) misclass, (2) mu_tilde, (3) tau2, (4) upsilon, (5) mu_i
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

            arma::vec it_indices = {k, k+1};
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
                                                     B_temp, it_indices);

                double log_target = log_f_i_cpp(i, ii, pars, prior_par, 
                                                par_index, y_1, t_pts, id,
                                                pr_B, it_indices);
                
                arma::vec col1(pr_B.n_elem, arma::fill::zeros);
                col1.elem(arma::find(pr_B == 1)).ones(); 
                arma::vec col2(pr_B.n_elem, arma::fill::zeros);
                col2.elem(arma::find(pr_B == 2)).ones(); 
                arma::vec col3(pr_B.n_elem, arma::fill::zeros);
                col3.elem(arma::find(pr_B == 3)).ones(); 
                
                pr_V = arma::join_horiz(col1, arma::join_horiz(col2, col3));
                
                // Note that the proposal probs cancel in the MH ratio
                double diff_check = log_target - log_target_prev;
                double min_log = log(arma::randu(arma::distr_param(0,1)));
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
arma::mat update_mu_i_cpp(const arma::vec &y_2, const arma::vec &pars,
                          const arma::field<arma::uvec> &par_index,
                          const arma::field<arma::mat> &V_i,
                          const arma::vec &EIDs, const arma::vec &id) {
    // par_index KEY: (0) beta, (1) misclass, (2) mu_tilde, (3) tau2, (4) upsilon, (5) mu_i
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    arma::mat mu_i(EIDs.n_elem, 3);
    
    // omp_set_num_threads(8);
    // # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        int i = EIDs(ii);
        
        arma::mat V_small = V_i(ii);
        
        arma::uvec vec_upsilon_ind = par_index(4);
        arma::vec vec_upsilon_content = pars.elem(vec_upsilon_ind - 1);
        arma::mat upsilon = arma::reshape(vec_upsilon_content, 3, 3);
        
        arma::mat up_solve = arma::inv_sympd(upsilon);
        
        arma::uvec sub_ind = arma::find(id == i);
        arma::vec y_sub = y_2.elem(sub_ind); 
        
        arma::vec vec_tau2 = pars.elem(par_index(3) - 1);
        double tau2 = vec_tau2(0);
        
        arma::vec mu_tilde = pars.elem(par_index(2) - 1);
        
        arma::mat V_inv = (1/tau2) * (V_small.t() * V_small) + up_solve;
        arma::mat V = arma::inv(V_inv);
        
        arma::vec M = V * ((1/tau2) * (V_small.t() * y_sub) + up_solve * mu_tilde);
        
        arma::mat mu_i_small = rmvnorm(1, M, V);
        
        mu_i.row(ii) = mu_i_small;
        
    }
    
    return mu_i;
}



// [[Rcpp::export]]
void test_functions(const arma::vec &pars, const arma::field<arma::vec> &prior_par, 
                    const arma::field<arma::uvec> &par_index) {
    
    // Multivariate Normal Check
    arma::uvec vec_beta_ind = par_index(0);
    arma::uvec vec_misclass_ind = par_index(1);
    
    arma::vec vec_beta_content = pars.elem(vec_beta_ind - 1);
    arma::vec vec_misclass_content = pars.elem(vec_misclass_ind - 1);
    
    arma::vec p_mean = prior_par(0);
    arma::mat p_sd = arma::diagmat(prior_par(1));
    arma::mat x = arma::join_cols(vec_beta_content, vec_misclass_content);
    
    arma::vec log_prior_dens = dmvnorm(x.t(), p_mean, p_sd, true);
    
    Rcpp::Rcout << log_prior_dens << std::endl;
    
    //  Sub setting check
    arma::vec id = {1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,3};
    arma::vec y_1 = 5 * id;
    arma::uvec sub_ind = arma::find(id == 2);
    
    int n_i = sub_ind.n_elem; // DOUBLE CHECK and compare to other method
    arma::vec y_1_sub = y_1.elem(sub_ind);
    
    Rcpp::Rcout << sub_ind << std::endl;
    Rcpp::Rcout << y_1_sub << std::endl;
    
    Rcpp::Rcout << "Attempt 1: " << n_i << std::endl;
    
    int n_i_2 = sub_ind.max() - sub_ind.min() + 1;
    Rcpp::Rcout << sub_ind.max() << " " << sub_ind.min() << std::endl;
    
    arma::vec col1(id.n_elem, arma::fill::zeros);
    col1.elem(arma::find(id == 1)).ones();
    arma::vec col2(id.n_elem, arma::fill::zeros);
    col2.elem(arma::find(id == 2)).ones();
    arma::vec col3(id.n_elem, arma::fill::zeros);
    col3.elem(arma::find(id == 3)).ones();
    
    arma::mat pr_V = arma::join_horiz(col1, arma::join_horiz(col2, col3));
    
    Rcpp::Rcout << pr_V << std::endl;
    
}