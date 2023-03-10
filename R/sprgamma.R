sprgamma <- function(spcov_params, dispersion = 1, mean = 0, samples = 1, data, randcov_params, partition_factor, ...) {

  n <- NROW(data)
  call_val <- match.call()
  call_val[[1]] <- as.symbol("sprnorm")
  call_list <- as.list(call_val)
  if ("dispersion" %in% names(call_list)) {
    call_list <- call_list[-which(names(call_list) == "dispersion")]
  }
  call_val <- as.call(call_list)
  sprnorm_val <- eval(call_val, envir = parent.frame())
  mu <- exp(sprnorm_val)

  if (is.matrix(mu)) {
    mu_list <- split(t(mu), seq_len(NCOL(mu)))
    sprgamma_val <- vapply(mu_list, function(x) rgamma(n, shape = dispersion, scale = x / dispersion), numeric(n))
  } else {
    sprgamma_val <- rgamma(n, shape = dispersion, scale = mu / dispersion)
  }
  sprgamma_val
}
