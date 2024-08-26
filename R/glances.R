#' Glance at many fitted model objects
#'
#' @description \code{glances()} repeatedly calls \code{glance()} on several
#'   fitted model objects and binds the output together, sorted by a column of interest.
#'
#' @param object A fitted model object from [splm()], [spautor()], [spglm()], or [spgautor()].
#' @param ... Additional fitted model objects. Ignored
#'   if \code{object} has class \code{splm_list}, \code{spautor_list}, \code{spglm_list}, or \code{spgautor_list}.
#' @param sort_by Sort by a \code{glance} statistic (i.e., the name of a column
#'   output from \code{glance()} or the order of model input (\code{sort_by = "order"}).
#'   The default is \code{"AICc"}.
#' @param decreasing Whether \code{sort_by} should sort by decreasing order? The default is \code{FALSE}.
#' @param warning Whether a warning is displayed when model comparisons violate certain rules.
#'   The default is \code{TRUE}.
#' @return A tibble where each row represents the output of \code{glance()} for
#'   each fitted model object.
#'
#' @order 1
#' @export
#'
#' @examples
#' lmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "none"
#' )
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' glances(lmod, spmod)
#' glances(lmod, spmod, sort_by = "logLik", decreasing = TRUE)
glances <- function(object, ...) {
  UseMethod("glances", object)
}

#' @rdname glances
#' @method glances splm
#' @order 2
#' @export
glances.splm <- function(object, ..., sort_by = "AICc", decreasing = FALSE, warning = TRUE) {
  model_list <- c(list(object), list(...))
  if (any(!(vapply(model_list, function(x) class(x), character(1)) %in% c("splm", "spautor")))) {
    stop("All models must be of class splm or spautor", call. = FALSE)
  }
  model_list_names <- c(as.character(as.list(substitute(list(object)))[-1]), as.character(as.list(substitute(list(...)))[-1]))
  if (warning && length(model_list) > 1) {
    check_likstat_use(model_list)
  }
  model_glance <- lapply(model_list, function(x) glance(x))
  model_bind <- do.call(rbind, model_glance)
  model_bind <- cbind(data.frame(model = model_list_names), model_bind)
  if (sort_by == "order") {
    model_bind <- model_bind[order(seq_len(NROW(model_bind)), decreasing = decreasing), , drop = FALSE]
  } else {
    model_bind <- model_bind[order(model_bind[[substitute(sort_by)]], decreasing = decreasing), , drop = FALSE]
  }
  tibble::as_tibble(model_bind)
}

#' @rdname glances
#' @method glances spautor
#' @order 3
#' @export
glances.spautor <- glances.splm

#' @rdname glances
#' @method glances splm_list
#' @order 4
#' @export
glances.splm_list <- function(object, ..., sort_by = "AICc", decreasing = FALSE, warning = TRUE) {
  model_glance <- lapply(object, function(x) glance(x))
  model_bind <- do.call(rbind, model_glance)
  model_bind <- cbind(data.frame(model = names(model_glance), model_bind))
  if (sort_by == "order") {
    model_bind <- model_bind[order(seq_len(NROW(model_bind)), decreasing = decreasing), , drop = FALSE]
  } else {
    model_bind <- model_bind[order(model_bind[[substitute(sort_by)]], decreasing = decreasing), , drop = FALSE]
  }
  tibble::as_tibble(model_bind)
}

#' @rdname glances
#' @method glances spautor_list
#' @order 5
#' @export
glances.spautor_list <- glances.splm_list


check_likstat_use <- function(model_list) {


  est_methods <- vapply(model_list, function(x) x$estmethod, character(1))
  n <- vapply(model_list, function(x) x$n, numeric(1))

  if (any(est_methods %in% c("ml", "reml"))) {
    if (length(unique(n)) > 1) {
      warning('Likelihood-based comparisons (e.g., AIC, AICc, BIC) should only be used to compare models that have the same response variable values (and sample size).', call. = FALSE)
    }
    # probably should also check that the response vectors (sorted) are actually equal
    # e.g., any(sort(model.response(model.frame(model1))) != sort(model.response(model.frame(model2))))
    # any problems here we don't foresee?
  }

  if (any("reml" %in% est_methods) && any("ml" %in% est_methods)) {
    warning('Likelihood-based comparisons (e.g., AIC, AICc, BIC) should not be used to compare a model fit with estmethod = "ml" to a model with estmethod = "reml".', call. = FALSE)
  }
  reml_model_list <- model_list[est_methods == "reml"]
  if (length(reml_model_list) > 1) {
    mm_names <- lapply(reml_model_list, function(x) sort(colnames(model.matrix(x))))
    if (any(!duplicated(mm_names)[-1])) { # drop first as it is always FALSE
      warning('Likelihood-based comparisons (e.g., AIC, AICc, BIC) should not be used to compare models fit with estmethod = "reml" when the models have distinct explanatory variable structures (i.e., distinct formula arguments).', call. = FALSE)
    }
  }
  if (inherits(model_list[[1]], c("spglm", "spgautor"))) {
    check_wrong_family(model_list)
  }
  # NULL
}

check_wrong_family <- function(model_list) {

  families <- vapply(model_list, function(x) x$family, character(1))

  wrong_family <- 0

  # binomial warning
  is_family_in <- families %in% "binomial"
  if (any(is_family_in) && any(!is_family_in)) {
    wrong_family <- wrong_family + 1
  }

  # beta warning
  is_family_in <- families %in% "beta"
  if (any(is_family_in) && any(!is_family_in)) {
    wrong_family <- wrong_family + 1
  }

  # count warning
  is_family_in <- families %in% c("poisson", "nbinomial")
  if (any(is_family_in) && any(!is_family_in)) {
    wrong_family <- wrong_family + 1
  }

  # skewed warning
  is_family_in <- families %in% c("Gamma", "inverse.gaussian")
  if (any(is_family_in) && any(!is_family_in)) {
    wrong_family <- wrong_family + 1
  }

  if (wrong_family > 0) {
    warning('Likelihood-based comparisons (e.g., AIC, AICc, BIC) should only be used to compare models fit with the same response variable support.', call. = FALSE)
  }

}
