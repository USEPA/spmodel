sprinvgauss <- function(spcov_params, dispersion = 1, mean = 0, samples = 1, data, randcov_params, partition_factor, ...) {

  if (!requireNamespace("statmod", quietly = TRUE)) {
    stop("Install the statmod package before using sprinvgauss", call. = FALSE)
  } else {
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
      sprinvgauss_val <- vapply(mu_list, function(x) {
        dispersion_true <- 1 / (x * dispersion)
        statmod::rinvgauss(n, mean = x, dispersion = dispersion_true)
      }, numeric(n))
    } else {
      dispersion_true <- 1 / (mu * dispersion)
      sprinvgauss_val <- statmod::rinvgauss(n, mean = mu, dispersion = dispersion_true)
    }
    sprinvgauss_val
  }

}
