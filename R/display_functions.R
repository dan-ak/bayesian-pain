#' Display Model Parameters
#'
#' @param stanfit
#' @return none
#'
#'@export
get_parameters <- function(data){
  stanfits <- data[[1]]
  for (i_m in 1:length(stanfits)){
    if (!is.null(stanfits[[i_m]])){
      fit <- stanfits[[i_m]]

      model_name <- fit@model_name
      model_params <- fit@model_pars[1:20]
      num_params <- 0
      for (i in 1:20){
        par <- model_params[i]

        if (substr(par,(nchar(par)+1)-2, nchar(par)) == "_n"){
          num_params <- i-1
          break
        }
      }
      model_params <- fit@model_pars[1:num_params]

      param_df <- rstan::summary(fit, pars = model_params)$summary[,c("mean", "se_mean")]
      print(model_name)
      if (num_params == 1){
        print(model_params)
      }
      print(param_df)
    }
  }
}

