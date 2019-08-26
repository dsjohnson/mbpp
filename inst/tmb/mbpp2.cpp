#include <TMB.hpp>

using namespace density;
using std::sqrt;

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
  PARAMETER(U);
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

  PARAMETER_MATRIX(M_avail);
  PARAMETER_MATRIX(U_avail);


  // DERIVED AND SIZE
  int n = m.size();
  int n_r = M.size();
  int n_o = M_avail.cols();
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

  for(int i=0; i<n_r; i++){
    for(int j; j<n_o; j++){
      nll -= dnorm(U*alpha)
    }
  }


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
