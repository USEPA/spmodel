#' Simulate a spatial binomial random variable
#'
#' @description Simulate a spatial binomial random variable with a specific
#'   mean and covariance structure.
#'
#' @param spcov_params An [spcov_params()] object.
#' @param size A numeric vector representing the sample size for each binomial trial.
#'   The default is \code{1}, which corresponds to a Bernoulli trial for each observation.
#' @param mean A numeric vector representing the mean. \code{mean} must have length 1
#'   (in which case it is recycled) or length equal
#'   to the number of rows in \code{data}. The default is \code{0}.
#' @param samples The number of independent samples to generate. The default
#'   is \code{1}.
#' @param data A data frame or \code{sf} object containing spatial information.
#' @param randcov_params A [randcov_params()] object.
#' @param partition_factor A formula indicating the partition factor.
#' @param ... Additional arguments passed to [sprnorm()].
#'
#' @details The values of \code{spcov_params}, \code{mean}, and \code{randcov_params}
#'   are assumed to be on the link scale. They are used to simulate a latent normal (Gaussian)
#'   response variable using [sprnorm()]. This latent variable is the
#'   conditional mean used with \code{dispersion} to simulate a binomial random variable.
#'
#' @return If \code{samples} is 1, a vector of random variables for each row of \code{data}
#'   is returned. If \code{samples} is greater than one, a matrix of random variables
#'   is returned, where the rows correspond to each row of \code{data} and the columns
#'   correspond to independent samples.
#'
#' @export
#'
#' @examples
#' spcov_params_val <- spcov_params("exponential", de = 0.2, ie = 0.1, range = 1)
#' sprbinom(spcov_params_val, data = caribou, xcoord = x, ycoord = y)
#' sprbinom(spcov_params_val, samples = 5, data = caribou, xcoord = x, ycoord = y)
sprbinom <- function(spcov_params, mean = 0, size = 1, samples = 1, data, randcov_params, partition_factor, ...) {
  n <- NROW(data)
  call_val <- match.call()
  call_val[[1]] <- as.symbol("sprnorm")
  call_list <- as.list(call_val)
  if ("size" %in% names(call_list)) {
    call_list <- call_list[-which(names(call_list) == "size")]
  }
  call_val <- as.call(call_list)
  sprnorm_val <- eval(call_val, envir = parent.frame())
  mu <- expit(sprnorm_val)

  if (is.matrix(mu)) {
    mu_list <- split(t(mu), seq_len(NCOL(mu)))
    sprbinom_val <- vapply(mu_list, function(x) rbinom(n, size, x), numeric(n))
  } else {
    sprbinom_val <- rbinom(n, size, mu)
  }
  sprbinom_val
}
