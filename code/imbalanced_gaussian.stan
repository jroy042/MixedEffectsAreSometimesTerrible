// Stan code for imbalanced gaussian mixed effects simulation

data {
  int<lower=2> N;  // number of observations
  int<lower=2> S;  // number of subjects

  int<lower=1> P;  // number of fixed effects
  int<lower=1,upper=P> QS;  // number of subject effects

  // sparse model matrix (CSR)
  int<lower=1> nz;  // number of non-zero elements in x
  vector[nz] x_w;  // non-zero elements in x
  int x_v[nz];  // column indices for x_w
  int x_u[N+1];  // row-start indices for x

  vector[N] y;  // continuous response
}

transformed data {
  int K;  // number of columns in x
  int SF;  // first subject effect column in x
  int SL;  // last subject effect column in x

  K = P + S * QS;
  SF = P + 1;
  SL = P + S * QS;
}

parameters {
  vector[P] beta_raw;
  real<lower=0> res_raw;

  matrix[QS,S] gamma_subj_raw;
  vector<lower=0>[QS] sigma_subj_raw;
  cholesky_factor_corr[QS] omega_subj_raw;
}

transformed parameters {
  vector[K] coef;  // all coefficients
  real<lower=0> res;  // residual standard error
  vector[N] y_hat;  // predicted log-odds

  // transform fixed effects
  coef[1:P] = beta_raw * 2;
  res = 0.5 * res_raw;

  // transform subject effects
  coef[SF:SL]
    = to_vector(rep_matrix(sigma_subj_raw,S)
      .* (omega_subj_raw * gamma_subj_raw));

  // y_hat = x * coef
  y_hat = csr_matrix_times_vector(N,K,x_w,x_v,x_u,coef);
}

model {
  beta_raw ~ normal(0,1);
  res_raw ~ normal(0,1);

  to_vector(gamma_subj_raw) ~ normal(0,1);
  sigma_subj_raw ~ normal(0,1);
  omega_subj_raw ~ lkj_corr_cholesky(2);

  y ~ normal(y_hat,res);  // linear model defined
}

generated quantities {
  real<lower=0> s0;
  real<lower=0> s1;
  real<lower=0> s2;
  real<lower=0> s3;
  real<lower=-1,upper=1> r01;
  real<lower=-1,upper=1> r02;
  real<lower=-1,upper=1> r03;
  real<lower=-1,upper=1> r12;
  real<lower=-1,upper=1> r13;
  real<lower=-1,upper=1> r23;

  s0 = sigma_subj_raw[1];
  s1 = sigma_subj_raw[2];
  s2 = sigma_subj_raw[3];
  s3 = sigma_subj_raw[4];
  {
    matrix[QS,QS] omega_subj;  // correlation in subject effects
    omega_subj = tcrossprod(omega_subj_raw);
    r01 = omega_subj[2,1];
    r02 = omega_subj[3,1];
    r03 = omega_subj[4,1];
    r12 = omega_subj[3,2];
    r13 = omega_subj[4,2];
    r23 = omega_subj[4,3];
  }
}
