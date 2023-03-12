#' Glance at many fitted model objects
#'
#' @description \code{glances()} repeatedly calls \code{glance()} on several
#'   fitted model objects and binds the output together, sorted by a column of interest.
#'
#' @param object A fitted model object from [splm()], [spautor()], [spglm()], or [spgautor()].
#' @param ... Additional fitted model objects. Ignored
#'   if \code{object} has class \code{splm_list} or \code{spautor_list}.
#' @param sort_by Sort by a \code{glance} statistic (i.e., the name of a column
#'   output from \code{glance()} or the order of model input (\code{sort_by = "order"}).
#'   The default is \code{"AICc"}.
#' @param decreasing Should \code{sort_by} be decreasing or not? The default is \code{FALSE}.
#'
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
glances <- function(object, ..., sort_by = "AICc", decreasing = FALSE) {
  UseMethod("glances", object)
}

#' @rdname glances
#' @method glances splm
#' @order 2
#' @export
glances.splm <- function(object, ..., sort_by = "AICc", decreasing = FALSE) {
  model_list <- c(list(object), list(...))
  if (any(!(vapply(model_list, function(x) class(x), character(1)) %in% c("splm", "spautor")))) {
    stop("All models must be of class splm or spautor", call. = FALSE)
  }
  model_list_names <- c(as.character(as.list(substitute(list(object)))[-1]), as.character(as.list(substitute(list(...)))[-1]))
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
glances.splm_list <- function(object, ..., sort_by = "AICc", decreasing = FALSE) {
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
