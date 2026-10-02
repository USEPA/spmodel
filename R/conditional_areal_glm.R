#' @rdname conditional
#' @method conditional spgautor
#' @export
conditional.spgautor <- function(object, newdata, output = "newdata",
                                type = c("link", "response", "new"), samples = 1000, newdata_size, ...) {
  output <- check_conditional_areal(object, output, samples, list(...))
  type <- match.arg(type)
  if (missing(newdata_size)) newdata_size <- NULL
  context <- get_prediction_object_spgautor(object, newdata, dispersion = NULL,
    newdata_size = newdata_size, local = FALSE)
  newdata_size <- context$newdata_size
  if (object$family == "binomial") {
    if (!is.numeric(newdata_size) || !length(newdata_size) || anyNA(newdata_size) ||
        any(!is.finite(newdata_size) | newdata_size < 0 | newdata_size != floor(newdata_size)) ||
        !length(newdata_size) %in% c(1L, NROW(context$newdata))) {
      stop("newdata_size must contain nonnegative integers, with length one or one per prediction row.", call. = FALSE)
    }
    newdata_size <- rep(newdata_size, length.out = NROW(context$newdata))
  }
  observed <- fitted(object, type = "link")
  if (identical(output, "object")) {
    return(conditional_areal_snapshot(object, observed, samples))
  }

  cov_context <- get_conditional_areal_cov(object)
  joint <- get_conditional_glm_joint(object, cov_lowchol = cov_context$cov_lowchol)
  wts_residual <- backsolve(t(cov_context$cov_lowchol), cov_context$SqrtSigInv_c0)
  draws <- draw_conditional_glm_joint(joint, samples, residual = TRUE)
  new_val <- draw_conditional_areal(context, cov_context, draws$residual, draws$beta,
    samples, wts_residual = wts_residual)
  if (type != "link") {
    new_val <- invlink_conditional(new_val, type, context$dispersion_params_val,
      object$family, newdata_size)
  }
  conditional_areal_output(object, new_val, draws$beta, observed, output, samples)
}
