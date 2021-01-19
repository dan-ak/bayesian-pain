#' Select Models from Hoskin2019
#'
#' Returns the stan model objects corresponding to the selected models from
#' the paper.
#'
#' @param model_numbers array of numbers corresponding to the desired models
#' @return an array of stan model objects corresponding to select models
#'
#' @export
get_models <- function(model_numbers){

  models <- c()

  model_names <- c("1_b", "2_b_multi", "3_br", "4_bre", "5_bremn", "6_bremn_var")

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
#'@export
fit_models <- function(models, stan_data, iter=300, chains=3){

  fits <- vector("list", length(models))
  DICs <- vector("list", length(models))

  i <- 0
  for (m in models){
    i <- i+1
    fits[[i]] <- rstan::stan(file = m, data = stan_data, iter = iter, chains = chains)

    if(substring(m, nchar(m)-7) == 'var.stan'){

      model_mle <- paste0(substring(m, 1 ,nchar(m)-5),
                          "_mle.stan")

      DICs[i] <- BayesianPain::calc_DIC_hierarchical(fits[[i]], model_mle, stan_data)
    } else{

      DICs[i] <- BayesianPain::calc_DIC(fits[[i]])
    }
  }

  return(list("fits" = fits, "DICs" = DICs))
}

