#' Centers and scales pain rating in experimental data.
#'
#' @param exp_data experimental data in the form of a dataframe object with fields (S,Q,Std,X,R).
#' S: subject id, Q: mean of presented que,
#' Std: standard deviation of presented que, X: administered stimulus, R: subject rating
#' @return dataframe object from \code{exp_data} with subject ratings that are centered and scaled
center_scale <- function(exp_data){

  clean_data = data.frame()

  for (s in unique(exp_data$S)){

    subject_data       = subset(exp_data, S == s)
    subject_fixed_data = subset(subject_data, Std == 0)

    mod                <- lm(X ~ R, data = subject_fixed_data)
    subject_data$R     <- predict(mod, data.frame(R = subject_data$R))
    clean_data         <- rbind(clean_data,subject_data)
  }

  return(clean_data)
}


#' Formats Data to Run with Stan Models
#'
#' @param exp_data experimental data in the form of a dataframe object with fields (S,Q,Std,X,R).
#' S: subject id, Q: mean of presented que,
#' Std: standard deviation of presented que, X: administered stimulus, R: subject rating
#' @return data to fit stan model
create_stan_data <- function(exp_data){

  stan_data <- list(N = length(exp_data$Q), subject= exp_data$S, q = exp_data$Q, std = exp_data$Std, X = exp_data$X, r = exp_data$R,
                  v = apply(do.call('rbind', strsplit(as.character(exp_data$V),',',fixed=TRUE)), 2, as.numeric),
                  w = apply(do.call('rbind', strsplit(as.character(exp_data$W),',',fixed=TRUE)), 2, as.numeric))
  return(stan_data)
}
