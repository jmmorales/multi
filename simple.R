library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)

# let's simulate some data
set.seed(123)

K = 5       # number of "boxes" (plant species)
N = 50      # number of observations
D = 3       # number of predictive variables
size = 100  # number of objects (interactions) put into K boxes

X = matrix(NA, K, D)
for(i in 1:D){
  m = seq(1:K) # mean of a predictor variable such as fruit size
  X[,i] = rnorm(K, m, 1) 
}

beta = matrix(runif(N * D) * 2, D, N)

x_beta = X %*% beta

Y <- matrix(NA, nrow = N, ncol = K)
for(i in 1:N){
  Y[i,] <- rmultinom(1, size=size, prob=exp(x_beta[,i]))
}

stan_dat <- list(
  K = K,
  N = N,
  D = D,
  y = Y,
  x = X
)

# pars <- c("Omega", "tau", "betas", "rho", "z",
#           "taup", "ps", "zp", "rhop", "sigmae",
#           "rhosq", "etasq", "sigmaee", "pa")


mod <- cmdstan_model('multinomial_reg.stan')

fit <- mod$sample(
  data = stan_dat, 
  seed = 123, 
  chains = 4, 
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  thin = 1,
  refresh = 200 # print update every 500 iters
)

fit_summary = fit$summary()
idx = seq(2, N*D+1, by= D)
plot(fit_summary$mean[idx], beta[1,])
plot(fit_summary$mean[idx + 1], beta[2,])
plot(fit_summary$mean[idx + 2], beta[3,])
abline(0,1)

