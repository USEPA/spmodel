#' Create a random effects covariance parameter object
#'
#' @description Create a random effects covariance parameter object for use with other
#'   functions.
#'
#' @param ... A named vector (or vectors) whose names represent the name of each random
#'   effect and whose values represent the variance of each random effect. If
#'   unnamed, \code{nm} is used to set names.
#' @param nm A character vector of names to assign to \code{...}.
#'
#' @details Names of the random effects should match eligible names given to
#'   \code{random} in modeling functions. While with the \code{random}
#'   argument to these functions, an
#'   intercept is implicitly assumed, with \code{randcov_params}, an intercept must be
#'   explicitly specified. That is, while with \code{random}, \code{x | group}
#'   is shorthand for \code{(1 | group) + (x | group)}, with \code{randcov_params},
#'   \code{x | group} implies just \code{x | group}, which means that if \code{1 | group}
#'   is also desired, it must be explicitly specified.
#'
#' @return A named numeric vector of random effect covariance parameters.
#'
#' @export
#'
#' @examples
#' randcov_params(group = 1, subgroup = 2)
#' randcov_params(1, 2, nm = c("group", "subgroup"))
#' # same as
#' randcov_params("1 | group" = 1, "1 | subgroup" = 2)
randcov_params <- function(..., nm) {
  dotlist <- list(...)
  names_dotlist <- names(dotlist)
  if (length(names_dotlist) > 0) {
    # remove unnecessary names if they exist
    dotlist <- lapply(seq_along(names_dotlist), function(x) {
      if (length(dotlist[[x]]) > 1) {
        stop("Each random effect must have only one variance parameter.", call. = FALSE)
      }
      if (!is.null(names_dotlist[[x]]) && names_dotlist[[x]] != "") {
        unname(dotlist[[x]])
        names(dotlist[[x]]) <- names_dotlist[[x]]
      }
      dotlist[[x]]
    })
  }
  randcov_params_val <- unlist(dotlist)
  if (!missing(nm)) {
    names(randcov_params_val) <- nm
  }
  randcov_params_val
}
