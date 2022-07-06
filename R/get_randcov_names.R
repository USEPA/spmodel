#' Find the names of random effects and coerce them (if needed) to consistent
#'   structure (1 | ranef) or (x | ranef)
#'
#' @param random  A random effects formula (one sided e.g., ~ random effects)
#'   where random intercepts are specified by ~ group or ~ (1 | group) and
#'   random slosep are specified by ~ (x | group)
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

get_randcov_name <- function(label) {
  if (grepl("|", label, fixed = TRUE)) {
    new_label <- label
  } else {
    new_label <- paste("1", label, sep = " | ")
  }
  if (grepl("/", new_label, fixed = TRUE)) {
    bar_split <- unlist(strsplit(new_label, " | ", fixed = TRUE))
    dash_split <- unlist(strsplit(bar_split[[2]], "/", fixed = TRUE))
    front <- bar_split[[1]]
    backs <- dash_split
    new_label <- lapply(seq_along(backs), function(x) paste(front, paste(backs[seq(from = 1, to = x, by = 1)], collapse = ":"), sep = " | "))
    # fronts <- labels(terms(reformulate(bar_split[[1]])))
    # backs <- labels(terms(reformulate(bar_split[[2]])))
    # new_label <- unlist(lapply(backs, function(x) paste(front, x, sep = " | ")))
  }
  new_label <- unlist(lapply(new_label, function(x) get_randcov_label(x)))
  new_label
}

get_randcov_label <- function(label) {
  strsplits <- strsplit(label, " | ", fixed = TRUE)
  terms_fronts <- terms(reformulate(strsplits[[1]][[1]]))
  labels_fronts <- labels(terms_fronts)
  if (attr(terms_fronts, "intercept") == 1) {
    labels_fronts <- c("1", labels_fronts)
  }
  form_fronts <- lapply(labels_fronts, function(x) paste(x, strsplits[[1]][[2]], sep = " | "))
}
# next in get_rand_Zs need to separate to get left and right of | (removing whitespace)
# then can make appropriate random effect design matrices by
# filling a matrix with the left values then subsetting by the design matrix right side (which only includes zeros and ones)




# could use old version for partition names
get_partition_names <- function(partition_factor) {
  # get the partition formula and turn it into a named vector here
  labels_initial <- labels(terms(partition_factor))
  new_labels <- unlist(lapply(labels_initial, get_partition_name))
}

get_partition_name <- function(label) {
  if (grepl("|", label, fixed = TRUE)) {
    new_label <- label
  } else {
    new_label <- paste("1", label, sep = " | ")
  }
  # if (grepl("/", new_label, fixed = TRUE)) { # patition factor can't have / here so comment out
  #   bar_split <- unlist(strsplit(new_label, " | ", fixed = TRUE))
  #   fronts <- labels(terms(reformulate(bar_split[[1]])))
  #   backs <- labels(terms(reformulate(bar_split[[2]])))
  #   new_label <- unlist(lapply(fronts, function(x) paste(x, backs, sep = " | ")))
  # }
  new_label
}
