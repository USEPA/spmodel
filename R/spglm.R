spglm <- function(formula, family, data, spcov_type, xcoord, ycoord, spcov_initial,
                 dispersion_initial, estmethod = "reml", anisotropy = FALSE,
                 random, randcov_initial, partition_factor, local, ...) {

  # set exponential as default if nothing specified
  if (missing(spcov_type) && missing(spcov_initial)) {
    spcov_type <- "exponential"
    message("No spatial covariance type provided. Assuming \"exponential\".")
  }

  if (!missing(spcov_type) && !missing(spcov_initial)) {
    message("Both spcov_type and spcov_initial provided. spcov_initial overriding spcov_type.")
  }

  # iterate if needed
  if (!missing(spcov_initial) && is.list(spcov_initial[[1]])) {
    call_list <- as.list(match.call())[-1]
    penv <- parent.frame()
    spglm_out <- lapply(spcov_initial, function(x) {
      call_list$spcov_initial <- x
      do.call("spglm", call_list, envir = penv)
    })
    names(spglm_out) <- paste("spcov_initial", seq_along(spcov_initial), sep = "_")
    new_spglm_out <- structure(spglm_out, class = "spglm_list")
    return(new_spglm_out)
  } else if (!missing(spcov_type) && length(spcov_type) > 1) {
    call_list <- as.list(match.call())[-1]
    penv <- parent.frame()
    spglm_out <- lapply(spcov_type, function(x) {
      call_list$spcov_type <- x
      do.call("spglm", call_list, envir = penv)
    })
    names(spglm_out) <- spcov_type
    new_spglm_out <- structure(spglm_out, class = "spglm_list")
    return(new_spglm_out)
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

  # set spcov_initial
  if (missing(spcov_initial)) {
    spcov_initial <- spmodel::spcov_initial(spcov_type)
  }

  # perform checks to return errors
  spglm_checks(spcov_initial, !missing(xcoord), !missing(ycoord), estmethod, anisotropy, !missing(random))

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

  if (missing(local)) {
    local <- NULL
  }

  # non standard evaluation for x and y coordinates
  xcoord <- substitute(xcoord)
  ycoord <- substitute(ycoord)

  # get data object
  data_object <- get_data_object_spglm(
    formula, family, data, spcov_initial, xcoord, ycoord,
    estmethod, anisotropy, random, randcov_initial,
    partition_factor, local, ...
  )

  # parallel cluster if necessary
  if (data_object$parallel) {
    data_object$cl <- parallel::makeCluster(data_object$ncores)
    # invisible(clusterEvalQ(data_object$cl, library(Matrix)))
  }

  cov_est_object <- cov_estimate_laploglik_spglm(data_object, formula,
                                                 spcov_initial, dispersion_initial, estmethod,
                                                 optim_dotlist = get_optim_dotlist(...))

}
