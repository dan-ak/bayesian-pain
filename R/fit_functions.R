#' Select Models from Hoskin2019
#'
#' Returns the stan model objects corresponding to the selected models from
#' the paper.
#'
#' @param model_numbers array of numbers corresponding to the desired models
#' @return an array of stan model objects corresponding to select models
#'
get_models <- function(model_numbers){

  models <- vector("list", length(model_numbers))

  model_names <- c("1_b", "2_br", "3_br", "4_bre", "5_bremn", "6_bremn_var")

  i <- 0
  for(m in model_numbers){
    i <- i+1
    models[i] <- system.file("stan", paste0("model_", model_names[m],".stan"),
                                      package = "BayesianPain")

  }
  return(models)
}



#' Fit Models
#'
#' Returns the stan model objects corresponding to the selected models from
#' the paper.
#'
#' @param model_numbers array of numbers corresponding to the desired models
#' @return an array of stan model objects corresponding to select models
#'
fit_models <- function(models, stan_data, iter=300, chains=3){

  fits <- vector("list", length(models))
  DICs <- vector("list", length(models))

  i <- 0
  for (m in models){
    i <- i+1
    fits[m] <- rstan::stan(file = m, data = stan_data, iter = iter, chains = chains)

    if(substring(m, nchar(m)-7) == 'var.stan'){

      model_mle <- paste0(substring(model_list[1], 1 ,nchar(model_list[1])-6),
                          "_mle.stan")

      DICs[m] <- BayesianPain::calc_DIC_hierarchical(fits[[m]], model_mle, stan_data)
    } else{

      DICs[m] <- BayesianPain::calc_DIC(fits[[m]])
    }
  }

  return(list("fits" = fits, "DICs" = DICs))
}
