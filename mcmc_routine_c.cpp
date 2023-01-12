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
    arma::mat M = { {1, 0, 0},
                    {0, 1, exp(vec_misclass_content(0))},
                    {0, exp(vec_misclass_content(0)), 1}};
    arma::vec m_row_sums = arma::sum(M, 1);
    M = M.each_col() / m_row_sums; // DOUBLE CHECK!!
    
    // The time-homogeneous probability transition matrix
    arma::uvec sub_ind = arma::find(eids == i);
    int n_i = sub_ind.max() - sub_ind.min() + 1;
    
    // Subsetting the data to relate only to this participant
    arma::mat b_i = B;
    arma::vec t_pts_sub = t_pts.elem(sub_ind);
    arma::vec y_1_sub = y_1.elem(sub_ind);
    
    // Full likelihood evaluation is not needed for updating pairs of b_i components
    for(int w=0; w < it_indices.n_elem;++w){
        arma::colvec z_i = {1, t_pts_sub(it_indices(w) - 1)}; // using the time point before
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
        
        arma::mat Q = { {-q1,         q1,         0},
                        { q2,   -q2 - q3,        q3},
                        { q4,         q5,  -q4 - q5}};
        Q = (t_pts_sub(it_indices(w)) - t_pts_sub(it_indices(w) - 1)) * Q;
        
        arma::mat P_i = arma::expmat(Q);
        int b_k_1 = b_i(it_indices(w)-1,0);
        int b_k = b_i(it_indices(w), 0);
        int y_1_k = y_1_sub(it_indices(w));
        in_value = in_value + log(P_i( b_k_1 - 1, b_k - 1)) + log(M(b_k - 1, y_1_k-1));
    }
    
    // Likelihood components from the other parts
    arma::vec p_mean = prior_par(0);
    arma::mat p_sd = arma::diagmat(prior_par(1));
    arma::mat x = arma::join_cols(vec_beta_content, vec_misclass_content);
    arma::vec log_prior_dens = dmvnorm(x, p_mean, p_sd, true);
    
    in_value = in_value + log_prior_dens(0);
    
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
        arma::vec it_indices = arma::linspace(1, n_i); // DOUBLE CHECK (always skip initial state)
        in_vals(ii) = log_f_i_cpp(i, ii, pars, prior_par, par_index, y_1, t_pts, id, B(ii), it_indices);
    }
    
    double in_value = arma::accu(in_vals);
    
    return in_value;
}

// [[Rcpp::export]]
Rcpp::List update_b_i_cpp(const int t, const arma::vec &EIDs, const arma::vec &pars,  
                          const arma::field<arma::vec> &prior_par, const arma::field<arma::uvec> &par_index,
                          const arma::vec &y_1, arma::vec t_pts, const arma::vec &id,
                          const arma::field <arma::vec> &B) {
    
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
        arma::uvec sub_ind = arma::find(eids == i);
        
        int n_i = sub_ind.n_elem;
        
        int rbc_rule = rbc_rule_vec(sub_ind.min());
        int clinic_rule = clinic_rule_vec(sub_ind.min());
        
        // Subsetting the data to only be a function of each ii is needed for parallel
        // Subsetting fields
        arma::vec B_temp = B(ii);
        arma::mat Dn_temp = Dn(ii);
        arma::vec A_temp = A(ii);
        arma::mat Xn_temp = Xn(ii);
        arma::sp_mat invKn_temp = invKn(ii);
        
        // Subsetting the remaining data
        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat z_temp = z.rows(sub_ind);
        
        // keep the last state at state 1 so we only iterate to n_i-1
        for (int k = 0; k < n_i - 2; k++) {
            
            arma::vec t_pts = arma::linspace(k+1, k+2, 2);
            arma::vec pr_B = B_temp;
            arma::mat pr_Dn = Dn_temp;
            
            // Sample and update the two neighboring states
            arma::mat Omega_set;
            if (clinic_rule >= 0) {
                Omega_set = Omega_fun_cpp_new(k + 1, n_i, B_temp, false);
            } else {
                Omega_set = Omega_fun_cpp_new(k + 1, n_i, B_temp, true);
            }
            
            int sampled_index = arma::randi(arma::distr_param(1, Omega_set.n_rows));
            
            pr_B.rows(k, k+1) = Omega_set.row(sampled_index-1).t();
            
            arma::vec b_i = pr_B;
            
            // Adding clinical review
            bool valid_prop = false;
            
            if(clinic_rule >= 0) {
                bool b_i_rule = arma::any(arma::vectorise(b_i)==2);
                if (clinic_rule == 1) {
                    if(b_i_rule) {valid_prop = true;}
                } else {
                    if (rbc_rule == 0 || (rbc_rule == 1 && b_i_rule)) {valid_prop = true;}
                }
            } else {  // evaluate the posterior regardless for clinic_rule=-1
                valid_prop = true;
            }
            
            if(valid_prop) {
                double log_target_prev = log_f_i_cpp(i, ii, t_pts, par, par_index,A_temp,B_temp,Y_temp,z_temp,Dn_temp,Xn_temp,invKn_temp);
                
                arma::vec twos(b_i.n_elem, arma::fill::zeros);
                arma::vec threes = twos; arma::vec fours = twos; arma::vec fives = twos;
                
                twos.elem(arma::find(b_i == 2)) += 1;
                threes.elem(arma::find(b_i == 3)) += 1;
                fours.elem(arma::find(b_i == 4)) += 1;
                fives.elem(arma::find(b_i == 5)) += 1;
                
                arma::vec ones(b_i.n_elem, arma::fill::ones);
                
                arma::mat bigB = arma::join_rows(ones, arma::cumsum(twos));
                bigB = arma::join_rows(bigB, arma::cumsum(threes));
                bigB = arma::join_rows(bigB, arma::cumsum(fours));
                bigB = arma::join_rows(bigB, arma::cumsum(fives));
                
                pr_Dn = arma::kron(arma::eye(4,4), bigB);
                double log_target = log_f_i_cpp( i,ii,t_pts,par,par_index,A_temp,pr_B,Y_temp,z_temp,pr_Dn,Xn_temp,invKn_temp);
                
                // Note that the proposal probs cancel in the MH ratio
                double diff_check = log_target - log_target_prev;
                double min_log = log(arma::randu(arma::distr_param(0,1)));
                if(diff_check > min_log){
                    B_temp = pr_B;
                    Dn_temp = pr_Dn;
                }
            }
        }
        B_return(ii) = B_temp;
        V_return(ii) = Dn_temp;
    }
    
    List B_Dn = List::create(B_return, V_return);
    
    return B_Dn;
}