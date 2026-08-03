#' Row group labels derived from a one-sided formula's model.matrix()
#'
#' @param reform A one-sided formula for the grouping variable(s)
#' @param data The data
#' @param xlev Optional factor levels to enforce (from \code{.getXlevels()}),
#'   so this call spans the same levels as a different (e.g. full) data set
#' @param na_pass If \code{TRUE}, use \code{na.action = na.pass} (for newdata-side
#'   calls, where a missing predictor value should surface as \code{NA} rather
#'   than error or drop the row); if \code{FALSE}, use the default \code{na.action}
#'   (for data/obdata-side calls, matching existing behavior there)
#'
#' @details \code{model.matrix()} one-hot encodes the (possibly multi-variable)
#'   grouping formula into dummy columns named with its own "varname + level"
#'   convention; for each row, \code{which()} finds the single column that's 1,
#'   and its column name becomes that row's group label. This is the single
#'   source of truth for that naming convention on the \code{model.matrix()}
#'   side -- used both for random effect grouping variables and the partition
#'   factor, at both fitting/observed-data time and newdata/prediction time.
#'
#' @return A character vector, one group label per row of \code{data}
#'
#' @noRd
model_matrix_group_labels <- function(reform, data, xlev = NULL, na_pass = FALSE) {
  mf <- if (na_pass) {
    model.frame(reform, data, na.action = na.pass, xlev = xlev)
  } else if (is.null(xlev)) {
    model.frame(reform, data)
  } else {
    model.frame(reform, data, xlev = xlev)
  }
  mx <- model.matrix(reform, mf)
  names_mx <- colnames(mx)
  split_rows <- split(mx, seq_len(NROW(mx)))
  names_mx[vapply(split_rows, function(y) which(as.logical(y)), numeric(1))]
}

#' Canonical "varname + level" column names for a single factor's levels
#'
#' @param varname The variable's name as it appears in the data/formula
#' @param levels The factor's levels, in level order (e.g. \code{levels(factor(...))})
#'
#' @details Matches \code{model.matrix()}'s own naming convention
#'   (\code{model_matrix_group_labels()} above), without building a
#'   (potentially large, dense) \code{model.matrix()} -- for use alongside
#'   \code{Matrix::fac2sparse()}-based sparse construction, where a dense
#'   \code{model.matrix()} call is exactly what's being avoided.
#'
#' @return A character vector of column names, one per level, in level order
#'
#' @noRd
get_factor_level_names <- function(varname, levels) {
  paste0(varname, levels)
}
