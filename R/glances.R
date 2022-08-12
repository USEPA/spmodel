#' Glance at many fitted model objects
#'
#' @description \code{glances()} repeatedly calls \code{glance()} on several
#'   fitted model objects and binds the output together, sorted by a column of interest.
#'
#' @param ... Fitted model objects from [splm()] or [spautor()].
#' @param sort_by Sort by a \code{glance} statistic. The default is \code{AICc}.
#' @param decreasing Should \code{sort_by} be decreasing or not? The default is \code{FALSE}.
#'
#' @return A tibble where each row represents the output of \code{glance()} for
#'   each fitted model object.
#'
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
glances <- function(..., sort_by = "AICc", decreasing = FALSE) {
  model_list <- list(...)
  model_list_names <- as.character(as.list(substitute(list(...)))[-1])
  model_glance <- lapply(model_list, function(x) glance(x))
  model_bind <- do.call(rbind, model_glance)
  model_bind <- cbind(data.frame(model = model_list_names), model_bind)
  model_bind <- model_bind[order(model_bind[[substitute(sort_by)]], decreasing = decreasing), , drop = FALSE]
  tibble::as_tibble(model_bind)
}
