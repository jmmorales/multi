functions { 
  /* compute the kronecker product
  * Copied from brms: Paul-Christian Bürkner (2018). 
  * Advanced Bayesian Multilevel Modeling with the R Package brms. 
  * The R Journal, 10(1), 395-411. <doi:10.32614/RJ-2018-017>
  * Args: 
    *   A,B: matrices 
  * Returns: 
    *   kronecker product of A and B
  */ 
    matrix kronecker(matrix A, matrix B) { 
      matrix[rows(A)*rows(B), cols(A)*cols(B)] kron; 
      for (i in 1:cols(A)) { 
        for (j in 1:rows(A)) { 
          kron[((j-1)*rows(B)+1):(j*rows(B)), ((i-1)*cols(B)+1):(i*cols(B))] = A[j,i] * B;
        } 
      } 
      return kron; 
    } 
  
  // copied from R. McElreath (2018). 
  // Statistical rethinking: A Bayesian course with examples in R and Stan. 
   matrix cov_GPL2(matrix x, real sq_alpha, real sq_rho, real delta) {
    int N = dims(x)[1];
    matrix[N, N] K;
    for (i in 1:(N-1)) {
      K[i, i] = sq_alpha + delta;
      for (j in (i + 1):N) {
        K[i, j] = sq_alpha * exp(-sq_rho * square(x[i,j]) );
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = sq_alpha + delta;
    return K;
  }    
}


data {
  int<lower = 2> N;         // number of posible outcomes (plant spp)
  int<lower = 1> N_J;       // number of groups (bird spp)
  int<lower = 1> L_J;       // number group-level predictors (bird traits)
  int<lower = 1> K;         // number of covariates
  array[N_J, N] int<lower=0> y;   // visits of birds to plants
  matrix[N, K] x;           // plant covariates
  matrix[N_J, L_J] TT;      // bird traits
  matrix[N_J, N_J] C;       // bird phylogenetic correlation matrix
  vector[N_J] ones;         // vector of 1s
}

// transformed data {
//   row_vector[D] zeros = rep_row_vector(0, D);
// }

parameters {
  corr_matrix[K] Omega;       // correlation matrix for regression parameters
  vector<lower=0>[K] tau;     // scale for var-cov of betas(parameters)
  vector[N_J * K] beta;         // regression coefficients
  real<lower=0,upper=1> rho;  // phylogenetic effect
  vector[L_J * K] z;            // coefficients for trait effects on reg pars
//  vector[N_1] r_1_1;        //  plant-level effects
  real<lower=0> etasqp;
  real<lower=0> rhosqp;
  real<lower=0> deltap;
  real<lower=0> phi;         //  NB overdispersion parameter
  // matrix[K-1, D] beta_raw;
  // matrix[D, N] beta;
}


transformed parameters { 
  matrix[K, K] Sigma = quad_form_diag(Omega, tau);
  matrix[N_J*K, N_J*K] S = kronecker(Sigma, rho * C + (1-rho) * diag_matrix(ones));
  matrix[L_J, K] Z = to_matrix(z, L_J, K);
  vector[N_J * K] m = to_vector(TT * Z); 
  matrix[N_J, K] b_m = to_matrix(beta, N_J, K);
}


model {
  matrix[N, N_J] x_beta = x * b_m';
  Omega ~ lkj_corr(2);
  tau ~ student_t(3,0,10); // cauchy(0, 2.5); // lognormal()
  beta ~ multi_normal(m, S);
  rho ~ beta(2,2);
  z ~ normal(0,5);
  
  to_vector(beta) ~ normal(0, 1);

  for (i in 1:N_J) {
    y[i,] ~ multinomial(softmax(x_beta[,i]));
  }
}
