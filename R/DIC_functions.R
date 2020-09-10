#' Calculates DIC
#'
#' Uses the samples of the log likelihoods in a stan fit file to calculate the
#' deviance information criterion (DIC).
#'
#' @param stan_fit stan fit object with \code{log_lik_t} field
#' @return DIC value
#'
#' @export
calc_DIC <- function(stan_fit){

  MCMC_LL_mean <- mean(rstan::extract(stan_fit, 'log_lik_t')$log_lik_t)
  MLE_LL <- max(rstan::extract(stan_fit, 'log_lik_t')$log_lik_t)

  DIC <- -2*(2*MCMC_LL_mean - MLE_LL)

  return(DIC)
}



#' Find the Maximum Likelihood Estimate of Parameter
#'
#' Finds the Maximum Likelihood (MLE) of a given parameter in a stan_fit object
#'
#' @param stan_fit stan fit object
#' @param parameter parameter of interest
#' @return the MLE of the parameter
get_param_mle <- function(stan_fit, parameter){
  z <- density(rstan::extract(stan_fit, parameter)[[parameter]], bw="SJ")
  return(z$x[z$y==max(z$y)])
}

#' Add MLE Parameter Values
#'
#' Not sure yet.
#'
#' @param stan_fit stan fit object
#' @param run_data experimental data formated for stan
#' @return run_data obejct with mle parameters included
add_mle_param <- function(stan_fit, run_data){

  model_params <- stan_fit@model_pars
  run_data_mle <- run_data

  for(param in model_params){
    if ((substring(param, nchar(param)-1) == '_u' || substring(param, nchar(param)-1) == '_s')
        && param != 'log_lik_s'){
      run_data_mle[[param]] <- get_param_mle(stan_fit, param)
    }
  }
  return(run_data_mle)
}


integrate_log_lik <- function(stan_fit){
  log_lik = rstan::extract(stan_fit, 'log_lik_s')$log_lik_s
  max_log_lik = apply(log_lik, MARGIN = 2, function(x) max(x, na.rm=TRUE))
  avg_log_lik = log(colMeans(exp(sweep(log_lik,2,max_log_lik)))) + max_log_lik
  return(sum(avg_log_lik))
}


#' Calculates DIC for Hierarchical Models
#'
#' Uses the
#'
#' @param stan_fit with \code{log_lik_t} field
#' @return DIC value
#'
#' @export
calc_DIC_hierarchical <- function(stan_fit, stan_mle, run_data){

  run_data_mle <- add_mle_param(stan_fit, run_data)
  MLE_fit <-stan(file = stan_mle, data = run_data_mle, iter = 200, chains = 3)

  MCMC_LL_mean <- integrate_log_lik(stan_fit)
  MLE_LL <- integrate_log_lik(MLE_fit)

  DIC <- -2*(2*MCMC_LL_mean - MLE_LL)

  return(DIC)
}
