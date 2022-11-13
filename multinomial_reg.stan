
data {
  int<lower = 2> K; // number of posible outcomes (plant spp)
  int<lower = 0> N; // bird spp
  int<lower = 1> D; // number of predictors
  int<lower = 0> y[N, K];
  matrix[K, D] x;
}

// transformed data {
//   row_vector[D] zeros = rep_row_vector(0, D);
// }

parameters {
  // matrix[K-1, D] beta_raw;
  matrix[D, N] beta;
}

// transformed parameters {
//   matrix[K, D] beta;
//   beta = append_row(beta_raw, zeros);
// }

model {
  matrix[K, N] x_beta = x * beta;
  
  to_vector(beta) ~ normal(0, 2);

  for (i in 1:N) {
    y[i,] ~ multinomial(softmax(x_beta[,i]));
  }
}