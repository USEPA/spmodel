sprpois <- function(spcov_params, mean = 0, samples = 1, data, randcov_params, partition_factor, ...) {

  n <- NROW(data)
  call_val <- match.call()
  call_val[[1]] <- as.symbol("sprnorm")
  sprnorm_val <- eval(call_val, envir = parent.frame())
  mu <- exp(sprnorm_val)

  if (is.matrix(mu)) {
    mu_list <- split(t(mu), seq_len(NCOL(mu)))
    sprpois_val <- vapply(mu_list, function(x) rpois(n, x), numeric(n))
  } else {
    sprpois_val <- rpois(n, mu)
  }
  sprpois_val

}

