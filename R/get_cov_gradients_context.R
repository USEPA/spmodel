# everything downstream (get_vcov_theta(), get_grad_g(), get_dSig_dtheta_cov())
# needs to treat the fitted covariance matrix Sigma as a function of its free
# parameters theta -- numerically differentiating the likelihood/covariance,
# or evaluating a closed-form derivative, at exactly the fitted optimum. This
# context object is a self-contained snapshot of every quantity that requires
# (rebuilt fresh from the fitted object rather than stored on it, since it is
# only needed for this satterthwaite()/emmeans-related machinery, not for
# ordinary model fitting/prediction)
get_cov_gradients_context <- function(object) {
  UseMethod("get_cov_gradients_context")
}

#' @exportS3Method
get_cov_gradients_context.splm <- function(object) {

  spcov_params_val <- coef(object, type = "spcov")
  spcov_type <- class(spcov_params_val)
  spcov_is_known <- object$is_known$spcov

  # rebuild an spcov_initial() object pinned at the fitted values (not
  # whatever starting values the original optimization used) -- this is the
  # point every downstream numerical derivative is taken around/near
  spcov_initial_val <- do.call(
    spcov_initial,
    c(
      list(spcov_type = spcov_type),
      as.list(spcov_params_val),
      list(known = names(spcov_is_known)[spcov_is_known])
    )
  )

  has_randcov <- !is.null(object$random)
  if (has_randcov) {
    randcov_params_val <- coef(object, type = "randcov")
    randcov_is_known <- object$is_known$randcov
    randcov_initial_val <- do.call(
      randcov_initial,
      c(
        as.list(randcov_params_val),
        list(known = names(randcov_is_known)[randcov_is_known])
        )
      )
  } else {
    randcov_params_val <- NULL
    randcov_is_known <- NULL
    randcov_initial_val <- NULL
  }

  data_object <- get_data_object_splm(
    formula = object$formula, data = object$obdata, spcov_initial = spcov_initial_val,
    xcoord = object$xcoord, ycoord = object$ycoord, estmethod = object$estmethod,
    anisotropy = object$anisotropy, random = object$random, randcov_initial = randcov_initial_val,
    partition_factor = object$partition_factor, local = FALSE, range_constrain = FALSE
  )

  if (object$anisotropy) {
    dist_matrix_list <- NULL
    dist_matrix <- NULL
  } else if (inherits(spcov_params_val, c("none", "ie")) && !has_randcov) {
    dist_matrix_list <- NULL
    dist_matrix <- NULL
  } else {
    dist_matrix_list <- lapply(data_object$obdata_list, function(d) spdist(d, data_object$xcoord, data_object$ycoord))
    dist_matrix <- as.matrix(dist_matrix_list[[1]])
  }


  # deprofiling step: during fitting, the overall spatial variance is often
  # profiled out of the optimization (solved for in closed form given the
  # other parameters) rather than optimized directly; spcov_profiled = FALSE
  # here forces the full (non-profiled) parameter set instead, since
  # get_vcov_theta()'s numeric branch needs to differentiate the likelihood
  # with respect to every free parameter, profiled ones included
  spcov_orig2optim_val <- spcov_orig2optim(spcov_initial_val, spcov_profiled = FALSE, data_object = data_object)
  randcov_orig2optim_val <- if (has_randcov) {
    randcov_orig2optim(randcov_initial_val, randcov_profiled = FALSE, spcov_initial = spcov_initial_val)
  } else {
    NULL
  }
  # eta: the fitted covariance parameters mapped onto the same unconstrained
  # ("optim") scale the original optimizer searched over (e.g. variances
  # log-transformed to be unconstrained) -- get_vcov_theta()'s numeric branch
  # differentiates the likelihood at this point, since a numeric Hessian is generally
  # well-behaved when away from the boundary of a constrained parameter space
  eta_val <- assemble_optim_par(spcov_orig2optim_val, randcov_orig2optim_val)

  spcov_names_free <- names(spcov_is_known)[!spcov_is_known]
  randcov_names_free <- if (has_randcov) names(randcov_is_known)[!randcov_is_known] else character(0)
  # order matters here: this concatenation (spcov names, then randcov names)
  # fixes the row/column order of every theta-indexed object built downstream
  # (vcov_theta, grad_g, the dSig_dtheta_cov() list) and is what
  # satterthwaite_core() relies on to split vcov_theta back into its spcov/
  # randcov blocks (see vcov(object, type = "spcov"/"randcov"))
  cov_names_free <- c(spcov_names_free, randcov_names_free)

  if (length(cov_names_free) == 0) {
    stop("All covariance parameters known. Satterthwaite not applicable.", call. = FALSE)
  }

  cov_val <- c(as.numeric(spcov_params_val), as.numeric(randcov_params_val))
  names(cov_val) <- c(names(spcov_params_val), names(randcov_params_val))
  cov_val_free <- cov_val[cov_names_free]

  # storing all the relevant context that will be used downstream
  context <- list(
    data_object = data_object,
    anisotropy = object$anisotropy, estmethod = object$estmethod,
    spcov_params = spcov_params_val, spcov_is_known = spcov_is_known,
    randcov_params = randcov_params_val, randcov_is_known = randcov_is_known,
    spcov_names_free = spcov_names_free,
    randcov_names_free = randcov_names_free, cov_names_free = cov_names_free,
    cov_val_free = cov_val_free,
    dist_matrix = dist_matrix,
    spcov_orig2optim = spcov_orig2optim_val, randcov_orig2optim = randcov_orig2optim_val,
    eta = eta_val
  )

  context
}

#' Build the covariance-gradient context for a fitted \code{spautor()} model
#'
#' Mirrors \code{get_cov_gradients_context.splm()}, but rebuilds the areal
#' (CAR/SAR) neighbor structure instead of a distance matrix. \code{W}/\code{M}
#' are reused exactly as already fitted (passed through with \code{row_st =
#' FALSE}) rather than recomputed, since re-row-standardizing an already
#' row-standardized \code{W} a second time would also silently (and
#' incorrectly) overwrite \code{M} -- see \code{build_car_neighbor_structure()}.
#'
#' @param object A fitted \code{spautor()} model object
#'
#' @return A context list, matching \code{get_cov_gradients_context.splm()}'s
#'   shape -- consumed by the \code{.spautor()} methods of
#'   \code{get_vcov_theta()}/\code{get_grad_g()}/\code{get_grad_gi()}/
#'   \code{get_dSig_dtheta_cov()}, which implement the areal (CAR/SAR)
#'   likelihood/covariance machinery where it genuinely differs from
#'   \code{splm()}'s
#'
#' @noRd
#' @exportS3Method
get_cov_gradients_context.spautor <- function(object) {

  spcov_params_val <- coef(object, type = "spcov")
  spcov_type <- class(spcov_params_val)
  spcov_is_known <- object$is_known$spcov

  spcov_initial_val <- do.call(
    spcov_initial,
    c(
      list(spcov_type = spcov_type),
      as.list(spcov_params_val),
      list(known = names(spcov_is_known)[spcov_is_known])
    )
  )

  has_randcov <- !is.null(object$random)
  if (has_randcov) {
    randcov_params_val <- coef(object, type = "randcov")
    randcov_is_known <- object$is_known$randcov
    randcov_initial_val <- do.call(
      randcov_initial,
      c(
        as.list(randcov_params_val),
        list(known = names(randcov_is_known)[randcov_is_known])
      )
    )
  } else {
    randcov_params_val <- NULL
    randcov_is_known <- NULL
    randcov_initial_val <- NULL
  }

  # range_positive is not stored on a fitted spautor object -- infer it from
  # the sign of the fitted autocorrelation parameter (the bound that must
  # have been active for the optimizer to land there). This only affects the
  # shape of the optim-scale transform used for the delta method below, not
  # the covariance matrix itself, which is built from the already-fitted W/M
  # regardless.
  range_positive_val <- isTRUE(spcov_params_val[["range"]] >= 0)

  # M is only meaningful for car (the CAR symmetry condition); sar ignores it
  # and build_car_neighbor_structure() warns if a non-NULL M is passed for a
  # sar model, so it is only forwarded for car here
  M_val <- if (inherits(spcov_params_val, "car")) object$M else NULL

  data_object <- get_data_object_spautor(
    formula = object$formula, data = object$data, spcov_initial = spcov_initial_val,
    estmethod = object$estmethod, W = object$W, M = M_val, random = object$random,
    randcov_initial = randcov_initial_val, partition_factor = object$partition_factor,
    row_st = FALSE, range_positive = range_positive_val, cutoff = NULL
  )

  spcov_orig2optim_val <- spcov_orig2optim(spcov_initial_val, spcov_profiled = FALSE, data_object = data_object)
  randcov_orig2optim_val <- if (has_randcov) {
    randcov_orig2optim(randcov_initial_val, randcov_profiled = FALSE, spcov_initial = spcov_initial_val)
  } else {
    NULL
  }
  eta_val <- assemble_optim_par(spcov_orig2optim_val, randcov_orig2optim_val)

  spcov_names_free <- names(spcov_is_known)[!spcov_is_known]
  randcov_names_free <- if (has_randcov) names(randcov_is_known)[!randcov_is_known] else character(0)
  cov_names_free <- c(spcov_names_free, randcov_names_free)

  if (length(cov_names_free) == 0) {
    stop("All covariance parameters known. Satterthwaite not applicable.", call. = FALSE)
  }

  cov_val <- c(as.numeric(spcov_params_val), as.numeric(randcov_params_val))
  names(cov_val) <- c(names(spcov_params_val), names(randcov_params_val))
  cov_val_free <- cov_val[cov_names_free]

  # storing all the relevant context that will be used downstream
  context <- list(
    data_object = data_object,
    anisotropy = FALSE, estmethod = object$estmethod,
    spcov_params = spcov_params_val, spcov_is_known = spcov_is_known,
    randcov_params = randcov_params_val, randcov_is_known = randcov_is_known,
    spcov_names_free = spcov_names_free,
    randcov_names_free = randcov_names_free, cov_names_free = cov_names_free,
    cov_val_free = cov_val_free,
    spcov_orig2optim = spcov_orig2optim_val, randcov_orig2optim = randcov_orig2optim_val,
    eta = eta_val
  )

  context
}
