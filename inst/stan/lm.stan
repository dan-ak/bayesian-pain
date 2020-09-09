//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// Save this file as inst/stan/lm.stan
data {
  int<lower=1> N;
  vector[N] x;
  vector[N] y;
}
parameters {
  real intercept;
  real beta;
  real<lower=0> sigma;
}
model {
  // ... priors, etc.

  y ~ normal(intercept + beta * x, sigma);
}
