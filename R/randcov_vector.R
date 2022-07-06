#' Create a random effects covariance vector
#'
#' @param randcov_params A \code{cov_params} object
#' @param data Data
#' @param newdata Newdata (used for prediction)
#'
#' @return A random effects covariance vector
#'
#' @noRd
randcov_vector <- function(randcov_params = NULL, data, newdata,
                           reform_bar2_list = NULL, Z_index_data_list = NULL,
                           reform_bar1_list = NULL, Z_data_list = NULL) {
  if (is.null(randcov_params)) {
    randcov_vectors <- NULL
  } else {
    randcov_names <- names(randcov_params)
    randcov_vectors <- lapply(
      randcov_names, get_randcov_vectors, randcov_params, data, newdata,
      reform_bar2_list, Z_index_data_list, reform_bar1_list, Z_data_list
    )
    randcov_vectors <- Reduce("+", randcov_vectors)
  }
  randcov_vectors
}

get_randcov_vectors <- function(randcov_name, randcov_params, data, newdata, reform_bar2_list, Z_index_data_list,
                                reform_bar1_list, Z_val_data_list) {
  randcov_param <- randcov_params[randcov_name]
  bar_split <- unlist(strsplit(randcov_name, " | ", fixed = TRUE))
  if (is.null(reform_bar2_list)) {
    reform_bar2 <- reformulate(paste0("as.numeric(", bar_split[[2]], ")"), intercept = FALSE)
  } else {
    reform_bar2 <- reform_bar2_list[[randcov_name]]
  }
  if (is.null(Z_index_data_list)) {
    Z_index_data <- as.vector(model.matrix(reform_bar2, data))
  } else {
    Z_index_data <- Z_index_data_list[[randcov_name]]
  }
  # Z_index_newdata <- as.vector(model.matrix(reform_bar2, newdata))
  # Z_index <- vapply(Z_index_newdata, function(x) ifelse(x == Z_index_data, randcov_param, 0), numeric(length(Z_index_data)))
  Z_index_newdata <- as.vector(model.matrix(reform_bar2, model.frame(reform_bar2, newdata, na.action = na.pass)))
  Z_index <- vapply(Z_index_newdata, function(x) ifelse(x != Z_index_data | is.na(x), 0, randcov_param), numeric(length(Z_index_data))) # NA | TRUE = TRUE / NA | FALSE = NA
  Z_index <- Matrix(Z_index, sparse = TRUE)
  if (bar_split[[1]] != "1") {
    if (is.null(reform_bar1_list)) {
      reform_bar1 <- reformulate(paste0("as.numeric(", bar_split[[1]], ")"), intercept = FALSE)
    } else {
      reform_bar1 <- reform_bar1_list[[randcov_name]]
    }
    if (is.null(Z_val_data_list)) {
      Z_val_data <- as.vector(model.matrix(reform_bar1, data))
    } else {
      Z_val_data <- Z_val_data_list[[randcov_name]]
    }
    Z_val_newdata <- as.vector(model.matrix(reform_bar1, newdata))
    Z_halfcov <- sweep(Z_index, 2, Z_val_newdata, `*`)
    # above same as Z_index (cov where zero if diff grp) * Z_val_newdata (the value of the covariate in newdata)
    Z_cov <- sweep(Z_halfcov, 1, Z_val_data, `*`)
    # above same as Z_halfcov * Z_val_data (the value of the covariate ni data)
  } else {
    Z_cov <- Z_index
  }
  t(Z_cov)
}
