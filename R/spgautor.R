spgautor <- function(formula, family, data, spcov_type, spcov_initial, dispersion_initial,
                    estmethod = "reml", random, randcov_initial, partition_factor, W, row_st = TRUE, M, ...) {

  # set car as default if nothing specified
  if (missing(spcov_type) && missing(spcov_initial)) {
    spcov_type <- "car"
    message("No spatial covariance type provided. Assuming \"car\".")
  }

  if (!missing(spcov_type) && !missing(spcov_initial)) {
    message("Both spcov_type and spcov_initial provided. spcov_initial overriding spcov_type.")
  }

  # iterate if needed
  if (!missing(spcov_initial) && is.list(spcov_initial[[1]])) {
    call_list <- as.list(match.call())[-1]
    penv <- parent.frame()
    spgautor_out <- lapply(spcov_initial, function(x) {
      call_list$spcov_initial <- x
      do.call("spgautor", call_list, envir = penv)
    })
    names(spgautor_out) <- paste("spcov_initial", seq_along(spcov_initial), sep = "_")
    new_spgautor_out <- structure(spgautor_out, class = "spgautor_list")
    return(new_spgautor_out)
  } else if (!missing(spcov_type) && length(spcov_type) > 1) {
    call_list <- as.list(match.call())[-1]
    penv <- parent.frame()
    spgautor_out <- lapply(spcov_type, function(x) {
      call_list$spcov_type <- x
      do.call("spgautor", call_list, envir = penv)
    })
    names(spgautor_out) <- spcov_type
    new_spgautor_out <- structure(spgautor_out, class = "spgautor_list")
    return(new_spgautor_out)
  }

  # set dispersion initial
  if (missing(dispersion_initial)) dispersion_initial <- NULL else family <- class(dispersion_initial)

  # fix family
  if (is.symbol(substitute(family))) { # or is.language
    family <- deparse1(substitute(family))
  }

  # Call splm if necessary
  if (family == "gaussian") {
    call_val <- match.call()
    call_val[[1]] <- as.symbol("splm")
    call_list <- as.list(call_val)
    call_list <- call_list[-which(names(call_list) %in% c("family", "dispersion_initial"))]
    call_val <- as.call(call_list)
    object <- eval(call_val, envir = parent.frame())
    return(object)
  }

  # replace initial values with appropriate NA's
  if (missing(spcov_initial)) {
    spcov_initial <- spmodel::spcov_initial(spcov_type)
  }

  spgautor_checks(class(spcov_initial), !missing(W), data, estmethod)

  # set partition factor if necessary
  if (missing(W)) {
    W <- NULL
  }

  if (missing(M)) {
    M <- NULL
  }

  # set random NULL if necessary
  if (missing(random)) {
    random <- NULL
  }

  # set rancov_initial NULL if necessary
  if (missing(randcov_initial)) {
    randcov_initial <- NULL
  }

  # set partition factor if necessary
  if (missing(partition_factor)) {
    partition_factor <- NULL
  }

  # get data object
  data_object <- get_data_object_spgautor(
    formula, family, data, spcov_initial,
    estmethod, W, M, random, randcov_initial,
    partition_factor, row_st, ...
  )

  cov_est_object <- cov_estimate_laploglik_spgautor(data_object, formula,
                                                 spcov_initial, dispersion_initial, estmethod,
                                                 optim_dotlist = get_optim_dotlist(...))

  model_stats <- get_model_stats_spgautor(cov_est_object, data_object, estmethod)

  output <- list(
    coefficients = model_stats$coefficients,
    fitted = model_stats$fitted,
    hatvalues = model_stats$hatvalues,
    residuals = model_stats$residuals,
    cooks_distance = model_stats$cooks_distance,
    vcov = model_stats$vcov,
    deviance = model_stats$deviance,
    pseudoR2 = model_stats$pseudoR2,
    p = data_object$p,
    n = data_object$n,
    npar = model_stats$npar,
    formula = formula,
    terms = data_object$terms,
    call = match.call(),
    estmethod = estmethod,
    data = data_object$data,
    newdata = data_object$newdata,
    optim = cov_est_object$optim_output,
    is_known = cov_est_object$is_known,
    W = data_object$W,
    M = data_object$M,
    partition_factor = partition_factor,
    random = random,
    observed_index = data_object$observed_index,
    missing_index = data_object$missing_index,
    contrasts = data_object$contrasts,
    xlevels = data_object$xlevels,
    is_sf = data_object$is_sf,
    sf_column_name = data_object$sf_column_name,
    crs = data_object$crs,
    family = family,
    w = model_stats$w
  )

  new_output <- structure(output, class = "spgautor")
  new_output

}
