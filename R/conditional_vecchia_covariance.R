# Prepare covariance and neighbor-search
get_conditional_vecchia_covariance <- function(object, data) {
  spcov <- coef(object, type = "spcov")
  coords <- as.matrix(data[, c(object$xcoord, object$ycoord), drop = FALSE])
  if (object$anisotropy) {
    transformed <- transform_anis(data, object$xcoord, object$ycoord, spcov[["rotate"]], spcov[["scale"]])
    coords <- cbind(transformed$xcoord_val, transformed$ycoord_val)
  }
  covariance <- object[c("coefficients", "random", "partition_factor", "anisotropy", "xcoord", "ycoord", "dim_coords", "diagtol")]
  class(covariance) <- class(object)
  # Reuse prediction's group/slope encodings, including nested random effects.
  # All candidate rows are already in data, so each site's group can be looked up
  # without rebuilding a model frame during the sequential searches.
  random <- get_extra_randcov_list(object, data, data)$randcov_terms
  partition <- get_extra_partition_list(object, data, data)$partition_index_obdata
  list(data = data, coords = coords, spcov = spcov, object = covariance,
    x = coords[, 1L], y = coords[, 2L], random = random,
    randcov = coef(object, type = "randcov"), partition = partition,
    structured = !is.null(object$random) || !is.null(object$partition_factor))
}

# Build the small block
conditional_vecchia_covariance <- function(context, rows) {
  if (context$structured) {
    object <- context$object
    object$obdata <- context$data[rows, , drop = FALSE]
    return(as.matrix(covmatrix(object)))
  }
  xy <- context$coords[rows, , drop = FALSE]
  distance <- spdist_vectors2(xy[, 1L], xy[, 2L], xy[, 1L], xy[, 2L], sparse = FALSE)
  as.matrix(cov_matrix2(context$spcov, dist_matrix = distance, diagtol = context$object$diagtol))
}

# Compute a score vector for site i; both observed neighborhoods reuse it.
conditional_vecchia_scores <- function(context, i, method) {
  distance2 <- (context$x - context$x[i])^2 + (context$y - context$y[i])^2
  if (method == "distance") {
    score <- distance2
  } else {
    random <- NULL
    if (length(context$random)) {
      random <- numeric(length(distance2))
      for (name in names(context$random)) {
        term <- context$random[[name]]
        same <- term$level_index_map[[term$group_label[i]]]
        contribution <- as.numeric(context$randcov[name])
        if (!is.null(term$slope_val)) contribution <- contribution * term$slope_val[i] * term$slope_val[same]
        random[same] <- random[same] + contribution
      }
    }
    score <- -abs(as.numeric(cov_vector(context$spcov, sqrt(distance2), random)))
  }
  if (!is.null(context$partition)) {
    score[context$partition$group_label != context$partition$group_label[i]] <- Inf
  }
  score[i] <- Inf
  score
}

# Select from eligible candidate indices, with deterministic ties by row index.
# A cached score vector avoids recomputing distances/covariances for N and K.
# method = "all" retains the full candidate set for exactness comparisons.
conditional_vecchia_neighbors <- function(context, i, pool, size, method, scores = NULL) {
  if (size == 0L && method != "all") return(integer())
  if (method == "all") return(pool[pool != i])
  if (is.null(scores)) scores <- conditional_vecchia_scores(context, i, method)
  score <- scores[pool]
  eligible <- is.finite(score)
  pool <- pool[eligible]
  score <- score[eligible]
  if (length(pool) <= size) return(pool)
  cutoff <- sort(score, partial = size)[size]
  candidates <- which(score <= cutoff)
  pool[candidates[order(score[candidates], pool[candidates])][seq_len(size)]]
}

# Order unique locations, then restore coincident observations as separate rows.
# Apply the requested method separately to observed and prediction sites.
conditional_vecchia_order <- function(coords, ordering = "maxmin") {
  n <- NROW(coords)
  if (n <= 1L) return(seq_len(n))
  if (!ordering %in% c("maxmin", "grts")) return(get_decorrelate_order(ordering, coords[, 1L], coords[, 2L])$order)
  key <- paste(format(coords[, 1L], digits = 17), format(coords[, 2L], digits = 17), sep = ":")
  unique_rows <- which(!duplicated(key))
  group <- match(key, key[unique_rows])
  order_unique <- if (length(unique_rows) == 1L) 1L else
    get_decorrelate_order(ordering, coords[unique_rows, 1L], coords[unique_rows, 2L])$order
  order(match(group, order_unique), seq_len(n))
}
