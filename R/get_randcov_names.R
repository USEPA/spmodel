#' Find the names of random effects and coerce them (if needed) to consistent
#'   structure (1 | ranef) or (x | ranef)
#'
#' @param random  A random effects formula (one sided e.g., ~ random effects)
#'   where random intercepts are specified by ~ group or ~ (1 | group) and
#'   random slopes are specified by ~ (x | group)
#'
#' @return Names of random effects
#'
#' @noRd
get_randcov_names <- function(random = NULL) {
  if (is.null(random)) {
    new_labels <- NULL
  } else {
    # get the random formula and turn it into a named vector here
    labels_initial <- labels(terms(random))
    new_labels <- unlist(lapply(labels_initial, get_randcov_name))
  }
  new_labels
}

#' Coerce a single random effect term label to \code{"(x | group)"} form
#'
#' @param label A single term label from \code{labels(terms(random))}
#'
#' @return The coerced label(s) (more than one if \code{label} contains a
#'   nesting operator \code{/}, which expands to one label per nesting level)
#'
#' @noRd
get_randcov_name <- function(label) {
  # a bare grouping variable (no "|") means a random intercept, so prepend "1 |"
  if (grepl("|", label, fixed = TRUE)) {
    new_label <- label
  } else {
    new_label <- paste("1", label, sep = " | ")
  }
  # "/" denotes nesting (e.g. group1/group2); expand it into one "x | group"
  # label per nesting level, with each level's grouping key built by
  # concatenating itself and all coarser levels with ":"
  if (grepl("/", new_label, fixed = TRUE)) {
    bar_split <- unlist(strsplit(new_label, " | ", fixed = TRUE))
    dash_split <- unlist(strsplit(bar_split[[2]], "/", fixed = TRUE))
    front <- bar_split[[1]]
    backs <- dash_split
    new_label <- lapply(seq_along(backs), function(x) paste(front, paste(backs[seq(from = 1, to = x, by = 1)], collapse = ":"), sep = " | "))
  }
  new_label <- unlist(lapply(new_label, function(x) get_randcov_label(x)))
  new_label
}

#' Expand a \code{"(x | group)"} label's left-hand side into individual terms
#'
#' @param label A single \code{"x | group"} random effect label
#'
#' @return One label per term on the left-hand side of \code{|} (including the
#'   intercept, if present), each still paired with the same grouping variable
#'
#' @noRd
get_randcov_label <- function(label) {
  strsplits <- strsplit(label, " | ", fixed = TRUE)
  # reformulate + terms() parses the left-hand side ("x") as its own formula
  # so multi-term slopes (e.g. "x1 + x2 | group") are split into one term each
  terms_fronts <- terms(reformulate(strsplits[[1]][[1]]))
  labels_fronts <- labels(terms_fronts)
  # terms() drops the intercept from labels(), so add it back explicitly if present
  if (attr(terms_fronts, "intercept") == 1) {
    labels_fronts <- c("1", labels_fronts)
  }
  form_fronts <- lapply(labels_fronts, function(x) paste(x, strsplits[[1]][[2]], sep = " | "))
}

#' Find the names of partition factor terms
#'
#' @param partition_factor A partition factor formula (one-sided, e.g. \code{~ group})
#'
#' @return Names of partition factor terms, coerced (if needed) to \code{"1 | group"} form
#'
#' @noRd
get_partition_names <- function(partition_factor) {
  labels_initial <- labels(terms(partition_factor))
  unlist(lapply(labels_initial, get_partition_name))
}

#' Coerce a single partition factor term label to \code{"1 | group"} form
#'
#' @param label A single term label from \code{labels(terms(partition_factor))}
#'
#' @return The coerced label
#'
#' @noRd
get_partition_name <- function(label) {
  if (grepl("|", label, fixed = TRUE)) {
    new_label <- label
  } else {
    new_label <- paste("1", label, sep = " | ")
  }
  new_label
}
