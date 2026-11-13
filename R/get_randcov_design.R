# Encode information without allocating unused combinations of factor levels.
randcov_design_key <- function(frame) {
  if (!NCOL(frame)) return(rep("", NROW(frame)))
  parts <- lapply(frame, function(x) {
    x <- as.character(x)
    if (anyNA(x)) stop("Random-effect groups and partitions must not be missing.", call. = FALSE)
    paste0(nchar(x, type = "bytes"), ":", x)
  })
  do.call(paste, c(parts, sep = "|"))
}

get_randcov_design <- function(observed, randcov_params, partition_factor = NULL, newdata = NULL, index = NULL) {
  if (is.null(index)) index <- seq_len(NROW(observed))
  if (is.null(newdata)) newdata <- observed[FALSE, , drop = FALSE]
  n <- NROW(observed)
  m <- NROW(newdata)
  terms <- names(randcov_params)
  Zo <- Zn <- list()
  info <- list()
  variance <- numeric()
  partition <- if (is.null(partition_factor)) data.frame(row.names = seq_len(n + m)) else {
    f <- reformulate(labels(terms(partition_factor)), intercept = FALSE)
    rbind(model.frame(f, observed, na.action = na.pass), model.frame(f, newdata, na.action = na.pass))
  }
  for (term in terms) {
    parts <- strsplit(term, " | ", fixed = TRUE)[[1L]]
    f <- reformulate(parts[2L], intercept = FALSE)
    groups <- rbind(model.frame(f, observed, na.action = na.pass), model.frame(f, newdata, na.action = na.pass))
    key <- randcov_design_key(cbind(groups, partition))
    labels <- unique(key[c(index, n + seq_len(m))])
    first <- match(labels, key)
    old <- which(first <= n)
    # Interaction levels vary fastest in their first factor; partitions follow levels.
    codes <- function(frame) rev(lapply(frame, function(x) {
      levels <- if (is.factor(x)) levels(x) else sort(unique(as.character(x[seq_len(n)])))
      match(as.character(x[first[old]]), levels)
    }))
    ordering <- do.call(order, c(codes(groups), codes(partition)))
    labels <- labels[c(old[ordering], which(first > n))]
    column <- match(key, labels)
    values <- if (parts[1L] == "1") rep(1, n + m) else {
      f <- reformulate(parts[1L], intercept = FALSE)
      c(as.numeric(model.matrix(f, model.frame(f, observed))),
        if (m) as.numeric(model.matrix(f, model.frame(f, newdata))) else numeric())
    }
    if (length(values) != n + m || any(!is.finite(values))) stop("Invalid random-effect slope values.", call. = FALSE)
    selected <- which(!is.na(column))
    Z <- Matrix::sparseMatrix(i = selected, j = column[selected], x = values[selected], dims = c(n + m, length(labels)))
    Zo[[term]] <- Z[index, , drop = FALSE]
    Zn[[term]] <- Z[n + seq_len(m), , drop = FALSE]
    first <- match(labels, key)
    level_names <- function(frame) Reduce(function(x, y) paste(x, y, sep = ":"),
      Map(get_factor_level_names, names(frame), lapply(frame[first, , drop = FALSE], as.character)))
    group_labels <- level_names(groups)
    if (NCOL(partition)) group_labels <- paste(group_labels, level_names(partition), sep = ":")
    partition_labels <- if (NCOL(partition)) apply(partition[first, , drop = FALSE], 1L, paste, collapse = ":") else rep(NA_character_, length(labels))
    row <- data.frame(term = rep(term, length(labels)), level = group_labels,
      partition = partition_labels, observed = labels %in% key[seq_len(n)],
      conditioned = as.numeric(Matrix::colSums(abs(Zo[[term]]))) > 0,
      stringsAsFactors = FALSE)
    row$level_tuple <- I(lapply(first, function(i) as.list(groups[i, , drop = FALSE])))
    row$partition_tuple <- I(lapply(first, function(i) as.list(partition[i, , drop = FALSE])))
    info[[term]] <- row
    variance <- c(variance, rep(as.numeric(randcov_params[term]), length(labels)))
  }
  if (!length(terms)) {
    metadata <- data.frame(term = character(), level = character(), partition = character(), observed = logical(), conditioned = logical())
    metadata$level_tuple <- I(list())
    metadata$partition_tuple <- I(list())
    return(list(observed = Matrix::Matrix(0, length(index), 0, sparse = TRUE),
      newdata = Matrix::Matrix(0, m, 0, sparse = TRUE), variance = numeric(), info = metadata))
  }
  metadata <- do.call(rbind, info)
  ids <- make.unique(paste0(metadata$term, ":", metadata$level))
  rownames(metadata) <- ids
  list(observed = do.call(cbind, Zo), newdata = do.call(cbind, Zn), variance = variance, info = metadata)
}

