sprbinom <- function(spcov_params, mean = 0, size = 1, samples = 1, data, randcov_params, partition_factor, ...) {

  n <- NROW(data)
  call_val <- match.call()
  call_val[[1]] <- as.symbol("sprnorm")
  call_list <- as.list(call_val)
  if ("size" %in% names(call_list)) {
    call_list <- call_list[-which(names(call_list) == "size")]
  }
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
