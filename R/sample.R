

exp_data <- BayesianPain::center_scale(data.frame(read.csv('data/exp_1_data.csv', header = TRUE)))
stan_data <- BayesianPain::create_stan_data(exp_data)
models <- BayesianPain::get_models(c(1,2,3,4,5))
output <- BayesianPain::fit_models(models, stan_data, iter=1000, chains = 4)
fit <- rstan::stan(file=models[[1]], data = stan_data, iter=1000, chains = 4)

util2 <- new.env()
source('R/stan_utility.R', local=util2)

util$check_all_diagnostics(fit)
