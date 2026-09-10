#' Get relevant L lists for anova
#'
#' @param assign_index A single assign value from the model matrix
#' @param assign_indices The assign values from the model matrix
#'
#' @return L lists for anova
#'
#' @noRd
get_L_list <- function(assign_index, assign_indices) {
  # assign_indices maps each model matrix column to the model term it belongs
  # to (as produced by attr(model.matrix(...), "assign")); find every column
  # belonging to this particular term
  assign_vals <- which(assign_indices == assign_index)
  # build one indicator row per column, then stack them into the L matrix
  # used to test this term's joint contrast (L %*% beta = 0) in anova()
  L_vectors <- lapply(assign_vals, get_L_vector, assign_indices)
  do.call(rbind, L_vectors)
}

#' Get a single indicator row vector for anova's L matrix
#'
#' @param assign_val The column index (into the model matrix) to indicate
#' @param assign_indices The assign values from the model matrix
#'
#' @return A single-row matrix of zeros with a one in column \code{assign_val}
#'
#' @noRd
get_L_vector <- function(assign_val, assign_indices) {
  L_vector <- matrix(0, nrow = 1, ncol = length(assign_indices))
  # a single 1 at position assign_val picks out that coefficient when
  # multiplied against the full coefficient vector (L %*% beta)
  L_vector[, assign_val] <- 1
  L_vector
}

#' Build an L list for anova based on a model object.
#'
#' @param L The L matrix (created if null)
#' @param Terms Model terms (created if null)
#' @param object Fitted model object
#'
#' @return A list of L matrices
#'
#' @noRd
get_L <- function(L, Terms, object) {
  # build a hypothesis matrix L (or list of them) and run a general linear
  # hypothesis test (GLHT) L*beta = 0 for each set of terms
  if (is.null(L)) {
    # "assign" attribute maps each column of the model matrix to the model
    # term that generated it (0 = intercept), used to group coefficients
    # belonging to the same term (e.g. all dummy columns of a factor)
    assign_indices <- attr(model.matrix(object), "assign") + 1
    # attr(model.matrix(object), "assign") if centering at zero
    if (is.null(Terms)) {
      # default: test each term separately (type III / marginal tests)
      assign_index <- unique(assign_indices)
      L <- lapply(assign_index, get_L_list, assign_indices)
      label <- labels(object)
      if (attr(terms(object), "intercept") == 1) {
        label <- c("(Intercept)", label)
      }
      names(L) <- label
    } else {
      # Terms specified: build one L testing the listed terms jointly
      if (is.character(Terms)) {
        Terms <- which(c("(Intercept)", labels(object)) %in% Terms) # - 1 if centering at zero
      }
      L <- list(do.call(rbind, lapply(Terms, get_L_list, assign_indices)))
      label <- c("(Intercept)", labels(object))
      label <- label[Terms] # label[Terms + 1] if centering at zero
      names(L) <- paste(label, collapse = ", ")
    }
  } else {
    # user supplied custom contrast matrix/matrices directly
    if (!is.list(L)) {
      L <- list(L)
    }
    names(L) <- paste("contrast", seq_along(L), sep = "")
  }
  L
}
