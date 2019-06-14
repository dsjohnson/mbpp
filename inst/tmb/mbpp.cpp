#include <TMB.hpp>

using namespace density;
using std::sqrt;

// VARIANCE FUNCTIONS
template<class Type>
Type sigma_i(vector<int> omega, vector<int> omega_p,
             Type alpha1, Type alpha2, Type delta1, Type delta2){
  if(omega(0)==omega_p(0)){
    Type ret_val = alpha1 * alpha2 * delta1 * delta2;
    return ret_val;
  } else{
    return Type(0);
  }
}

template<class Type>
Type sigma_ij(vector<int> omega, vector<int> omega_p,
              Type alpha, Type delta1, Type delta2){
  if((omega(0)==omega_p(0)) & (omega(1)==omega_p(1))){
    Type ret_val = alpha *(1-alpha) * delta1 * delta2;
    return ret_val;
  } else{
    return Type(0);
  }
}

template<class Type>
Type sigma_ijk(vector<int> omega, vector<int> omega_p, Type alpha, Type delta){
  if((omega(0)==omega_p(0)) & (omega(1)==omega_p(1)) & (omega(2)==omega_p(2))){
    Type ret_val = alpha * delta *(1-delta);
    return ret_val;
  } else{
    return Type(0);
  }
}

// UTILITY
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}


//MAIN OBJECTIVE FUNCTION
template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA
  DATA_IMATRIX(omega);
  DATA_VECTOR(M);
  DATA_VECTOR(m);
  DATA_VECTOR(u);
  DATA_VECTOR(D);
  DATA_MATRIX(X_a);
  DATA_MATRIX(X_d);
  // DATA_MATRIX(C_a);
  // DATA_MATRIX(C_d);
  DATA_IVECTOR(sig_a_idx);
  DATA_IVECTOR(sig_d_idx);

  //PARAMETERS
  PARAMETER_VECTOR(log_lambda);
  PARAMETER(logit_xi);
  PARAMETER_VECTOR(logit_tau);
  // alpha
  PARAMETER_VECTOR(theta_a);
  PARAMETER_VECTOR(log_sig_a);
  // PARAMETER(log_gamma_a);
  //delta
  PARAMETER_VECTOR(theta_d);
  PARAMETER_VECTOR(log_sig_d);
  // PARAMETER(log_gamma_d);


  // DERIVED AND SIZE
  int n = m.size();
  int n_r = M.size();
  //int n_a = theta_a.size();
  //int n_d = theta_d.size();

  vector<Type> lambda = exp(log_lambda);
  Type xi = invlogit(logit_xi);
  vector<Type> sig_a = exp(log_sig_a);
  vector<Type> sig_d = exp(log_sig_d);
  // Type gamma_a = exp(log_gamma_a);
  // Type gamma_d = exp(log_gamma_d);
  vector<Type> tau(n_r);
  for(int i=0; i<n_r; i++){
    tau(i) = invlogit(logit_tau(i));
  }
  vector<Type> logit_alpha = X_a * theta_a;
  vector<Type> logit_delta = X_d * theta_d;
  matrix<Type> Cov_m(n,n); Cov_m.setZero();
  vector<Type> mu_m(n);
  matrix<Type> Cov_u(n,n); Cov_u.setZero();
  vector<Type> mu_u(n);
  Type si, sij, sijk, alpha1, alpha2, delta1, delta2, tau1, M1, lambda1;

  for(int k=0; k<n; k++){
    alpha1 = invlogit(logit_alpha(k));
    delta1 = invlogit(logit_delta(k));
    tau1 = tau(omega(k,0));
    M1 = M(omega(k,0));
    lambda1 = lambda(omega(k,0));
    mu_m(k) =  M1 * alpha1 * delta1;
    mu_u(k) = lambda1 * (1-xi) * (1-tau1) * alpha1 * delta1;
    si = sigma_i(omega.row(k), omega.row(k), alpha1, alpha1, delta1, delta1);
    sij = sigma_ij(omega.row(k), omega.row(k), alpha1, delta1, delta1);
    sijk = sigma_ijk(omega.row(k), omega.row(k), alpha1, delta1);
    Cov_m(k,k) = M1 *(sij + sijk);
    Cov_u(k,k) = lambda1 * (1-xi) * (1-tau1) * (si + sij + sijk);
    for(int l=k+1; l<n; l++){
      alpha2 = invlogit(logit_alpha(l));
      delta2 = invlogit(logit_delta(l));
      si = sigma_i(omega.row(k), omega.row(l), alpha1, alpha2, delta1, delta2);
      sij = sigma_ij(omega.row(k), omega.row(l), alpha1, delta1, delta2);
      sijk = Type(0);
      Cov_m(k,l) = M1 * (sij + sijk);
      Cov_m(l,k) = Cov_m(k,l);
      Cov_u(k,l) = lambda1 * (1-xi) * (1-tau1) * (si + sij + sijk);
      Cov_u(l,k) = Cov_u(k,l);
    }
  }

  //LIKELIHOOD
  Type nll = 0.0;
  // Observed data models
  vector<Type> res_m = m - mu_m;
  vector<Type> res_u = u - mu_u;
  nll += MVNORM(Cov_m)(res_m);
  nll += MVNORM(Cov_u)(res_u);

  // xi prior
  nll -= invlogit(logit_xi) * (1-invlogit(logit_xi));

  for(int i=0; i<n_r; i++){
    // lambda prior
    nll -= dgamma(lambda(i), Type(0.00001), Type(100000), 1) + log_lambda(i);
    // tau prior
    nll -= log(invlogit(logit_tau(i))) + log(1-invlogit(logit_tau(i)));
    // marked model (M)
    nll -= dpois(M(i), (1-xi)*tau(i)*lambda(i), 1);
    // dead model (D)
    if(isNA(D(i))){
      nll -= 0.0;
    }else{
      nll -= dpois(D(i), xi*lambda(i), 1);
    }

  }

  //PRIOR,PENALTY ON PARAMETERS
  // matrix<Type> Cov_a = gamma_a * C_a;
  // matrix<Type> Cov_d = gamma_d * C_d;
  // nll += MVNORM(Cov_a)(theta_a);
  // nll += MVNORM(Cov_d)(theta_d);

  for(int i=0; i<theta_a.size(); i++){
    nll -= dnorm(theta_a(i), Type(0), sig_a(sig_a_idx(i)), 1);
  }
  for(int i=0; i<theta_d.size(); i++){
    nll -= dnorm(theta_d(i), Type(0), sig_d(sig_d_idx(i)), 1);
  }
  for(int i=0; i<sig_a.size(); i++){
    nll -= dexp(sig_a(i), Type(1), 1) + log_sig_a(i);
  }
  for(int i=0; i<sig_d.size(); i++){
    nll -= dexp(sig_d(i), Type(1), 1) + log_sig_d(i);
  }

  vector<Type> N_est(n_r);
  vector<Type> miss(n_r);
  for(int i=0; i<n_r; i++){
    if(isNA(D(i))){
      N_est(i) = M(i) + xi*lambda(i) + (1-xi)*(1-tau(i))*lambda(i);
      miss(i) = xi*lambda(i) + (1-xi)*(1-tau(i))*lambda(i);
    }else{
      N_est(i) = M(i) + D(i) + (1-xi)*(1-tau(i))*lambda(i);
      miss(i) = (1-xi)*(1-tau(i))*lambda(i);
    }
  }


  REPORT(theta_a);
  REPORT(theta_d);
  ADREPORT(lambda);
  ADREPORT(xi);
  ADREPORT(tau);
  ADREPORT(sig_a);
  ADREPORT(sig_d);
  ADREPORT(N_est);
  ADREPORT(miss);
  return nll;

}
