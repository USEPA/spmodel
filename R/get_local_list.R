#' Build the big-data local list used for covariance/fixed effect estimation
#'
#' @param local The user-supplied \code{local} argument (logical or list)
#' @param data A data frame or \code{sf} object
#' @param xcoord The x-coordinate name
#' @param ycoord The y-coordinate name
#' @param n The sample size
#' @param partition_factor A partition factor formula (or \code{NULL})
#'
#' @return A fully-specified \code{local} list, with defaults filled in for
#'   \code{index} (the partition assignment for each observation), \code{method},
#'   \code{var_adjust}, \code{parallel}, and \code{ncores} as needed. \code{size}
#'   can be set directly, or the number of \code{groups} can be set instead
#'   (in which case \code{size} is derived from it); \code{index}, if supplied,
#'   overrides both.
#'
#' @noRd
get_local_list_estimation <- function(local, data, xcoord, ycoord, n, partition_factor) {
  if (is.logical(local)) {
    if (local) {
      local <- list()
    } else {
      if (is.null(partition_factor)) {
        local <- list(index = rep(1, n))
      } else {
        # when local estimation is off but a partition factor is given, use the
        # partition factor's levels directly as the local index so estimation
        # is still split (but not approximated) along those groups
        index <- unname(model.response(model.frame(reformulate("1", response = labels(terms(partition_factor))), data = data)))
        index <- as.character(index) # turn into character if factor (this will also remove unused factor levels if there are any)
        local <- list(index = index)
        # resetting partition factor as NULL because it is in index but saving
        partition_factor <- NULL
      }
    }
  }

  names_local <- names(local)

  # errors
  if (!"index" %in% names_local && "method" %in% names_local) {
    if (!local$method %in% c("random", "kmeans")) {
      stop("Invalid local method. Local method must be \"random\" or \"kmeans\".", call. = FALSE)
    }
  }

  if (!"index" %in% names_local && "var_adjust" %in% names_local) {
    if (!local$var_adjust %in% c("none", "theoretical", "empirical", "pooled")) {
      stop("Invalid local var_adjust. Local var_adjust must be \"none\", \"theoretical\", \"empirical\", or \"pooled\".", call. = FALSE)
    }
  }

  if ("index" %in% names_local) {
    # if index is a factor and there are levels in the factor not in the observed
    # data, the code will fail. Storing as character prevents this (acts as droplevels)
    if (is.factor(local$index)) {
      local$index <- as.character(local$index)
    }
    local$size <- NULL
    local$groups <- NULL
    local$method <- NULL
  } else {
    # size (observations per partition) and groups (number of partitions) are
    # two ways to specify the same split; whichever one the user gave is used
    # to derive the other so both are always available downstream
    if (!"size" %in% names_local) {
      if ("groups" %in% names_local) {
        local$size <- ceiling(n / local$groups)
      } else {
        local$size <- 100
        local$groups <- ceiling(n / local$size)
      }
    } else {
      local$groups <- ceiling(n / local$size)
    }
    if (!"method" %in% names_local) {
      local$method <- "kmeans"
    }
    local$index <- get_local_estimation_index(local, data, xcoord, ycoord, n)
  }

  # setting var adjust
  if (!"var_adjust" %in% names_local) {
    if (n <= 100000) {
      local$var_adjust <- "theoretical"
    } else {
      message('var_adjust was not specified and the sample size exceeds 100,000, so the default var_adjust value is being changed from "theoretical" to "none". To override this behavior, rerun and set var_adjust in local. Be aware that setting var_adjust to "theoretical" may result in exceedingly long computational times.')
      local$var_adjust <- "none"
    }
  } # "none", "empirical", "theoretical", and "pooled"

  # setting partition factor
  local$partition_factor <- partition_factor

  if (!"parallel" %in% names_local) {
    local$parallel <- FALSE
    local$ncores <- NULL
  }

  if (local$parallel) {
    # never request more cores than there are partitions to farm out, or more
    # than the machine actually has available
    n_index <- length(unique(local$index))
    if ("ncores" %in% names_local) {
      cores_available <- parallel::detectCores()
      local$ncores <- min(n_index, local$ncores, cores_available)
    } else {
      local$ncores <- parallel::detectCores()
      local$ncores <- min(n_index, local$ncores)
    }
  }

  local
}

#' Assign each observation to a local-estimation partition index
#'
#' @param local A \code{local} list (must contain \code{method} and \code{groups})
#' @param data A data frame or \code{sf} object
#' @param xcoord The x-coordinate name
#' @param ycoord The y-coordinate name
#' @param n The sample size
#'
#' @return An integer/cluster vector of length \code{n} assigning each
#'   observation to one of \code{local$groups} partitions, via either random
#'   assignment or k-means clustering on the coordinates
#'
#' @noRd
get_local_estimation_index <- function(local, data, xcoord, ycoord, n) {
  if (local$method == "random") {
    # cycle group labels 1:groups enough times to cover n observations, then
    # shuffle -- gives partitions of (approximately) equal size at random
    index <- sample(rep(seq_len(local$groups), times = local$size)[seq_len(n)])
  } else if (local$method == "kmeans") {
    # any extra elements in local (beyond the reserved names below) are passed
    # straight through to kmeans(), e.g. to control algorithm or nstart
    kmeans_arg_names <- setdiff(names(local), c("size", "groups", "method", "index", "parallel", "ncores", "var_adjust"))
    kmeans_args <- local[kmeans_arg_names]
    # cluster on spatial coordinates so each partition is a compact neighborhood
    x <- cbind(data[[xcoord]], data[[ycoord]])
    index <- do.call("kmeans", c(list(x = x, centers = local$groups, iter.max = 30), kmeans_args))$cluster
  } else {
    stop("local$method must be random (the default) or kmeans")
  }
  index
}

#' Build the big-data local list used for point prediction
#'
#' @param local The user-supplied \code{local} argument (logical or list)
#'
#' @return A fully-specified \code{local} list, with defaults filled in for
#'   \code{method} (\code{"all"} for all data, \code{"distance"} for local
#'   distance neighborhoods, or \code{"covariance"} for local covariance
#'   neighborhoods), \code{size}, \code{byrow_threshold} (see
#'   \code{predict.spmodel()}), \code{parallel}, and \code{ncores} as needed
#'
#' @noRd
get_local_list_prediction <- function(local) {
  if (is.logical(local)) {
    if (local) {
      local <- list(method = "covariance", size = 100)
    } else {
      local <- list(method = "all")
    }
  }

  names_local <- names(local)

  # errors
  if ("method" %in% names_local) {
    if (!local$method %in% c("all", "covariance", "distance")) {
      stop("Invalid local method. Local method must be \"all\", \"covariance\", or \"distance\".", call. = FALSE)
    }
  }

  if ("byrow_threshold" %in% names_local) {
    if (!is.numeric(local$byrow_threshold) || length(local$byrow_threshold) != 1 || local$byrow_threshold < 0) {
      stop("local$byrow_threshold must be a single non-negative number.", call. = FALSE)
    }
  }


  if (!"method" %in% names_local) {
    local$method <- "covariance"
  }

  if (local$method %in% c("distance", "covariance") && !"size" %in% names_local) {
    local$size <- 100
  }

  # 10,000 x 10,000 matrix is about as large as we want to hold in memory
  if (!"byrow_threshold" %in% names_local) {
    local$byrow_threshold <- 10000^2
  }

  if (!"parallel" %in% names_local) {
    local$parallel <- FALSE
    local$ncores <- NULL
  }

  if (local$parallel) {
    if (!"ncores" %in% names_local) {
      local$ncores <- parallel::detectCores()
    }
  }

  local
}

#' Build the big-data local list used for block prediction
#'
#' @param local The user-supplied \code{local} argument (logical or list)
#'
#' @return A fully-specified \code{local} list, with defaults filled in for
#'   \code{method}, \code{size} (defaulting larger than point prediction's,
#'   since block prediction needs relatively more neighbors for the same
#'   accuracy), \code{parallel}, and \code{ncores} as needed
#'
#' @noRd
get_local_list_prediction_block <- function(local) {
  if (is.logical(local)) {
    if (local) {
      local <- list(method = "covariance", size = 4000)
    } else {
      local <- list(method = "all")
    }
  }

  names_local <- names(local)

  # errors
  if ("method" %in% names_local) {
    if (!local$method %in% c("all", "covariance", "distance")) {
      stop("Invalid local method. Local method must be \"all\", \"covariance\", or \"distance\".", call. = FALSE)
    }
  }


  if (!"method" %in% names_local) {
    local$method <- "covariance"
  }

  if (local$method %in% c("distance", "covariance") && !"size" %in% names_local) {
    local$size <- 4000
  }

  if (!"parallel" %in% names_local) {
    local$parallel <- FALSE
    local$ncores <- NULL
  }

  if (local$parallel) {
    if (!"ncores" %in% names_local) {
      local$ncores <- parallel::detectCores()
    }
  }

  local
}

#' Build the big data approximation settings for \code{\link{sprnorm}()}
#'
#' Fills in defaults for (and validates) the \code{local} argument used by
#' the unconditional simulation big data approximation: a random or
#' spatially-balanced (GRTS) ordering is used to draw a "base" sample, the
#' remaining locations are split into blocks (optionally via k-means on
#' coordinates), and \code{\link{get_conditional_new_from_base}()} later
#' simulates each block conditional on the base sample alone.
#'
#' @param local A logical or list; see the \code{local} argument to
#'   \code{\link{sprnorm}()}.
#' @param n The total number of locations to simulate.
#' @param data A data frame containing \code{...xcoord...}/\code{...ycoord...}
#'   columns (only used when \code{reorder = "grts"} or \code{kmeans = TRUE}).
#'
#' @return A list with the resolved \code{local} settings, including
#'   \code{index = list(base = ..., new = ...)} when \code{method != "all"}.
#'
#' @noRd
get_local_list_simulation <- function(local, n, data) {

  if (is.null(local)) {
    if (n > 5000) {
      local <- TRUE
      message("Because the desired simulation size exceeds 5,000, we are setting local = TRUE to perform computationally efficient approximations. To override this behavior and compute the exact solution, rerun with local = FALSE. Be aware that setting local = FALSE may result in exceedingly long computational times.")
    } else {
      local <- FALSE
    }
  }

  if (is.logical(local)) {
    if (local) {
      local <- list()
    } else {
      local <- list(method = "all")
    }
  }

  names_local <- names(local)

  if (!"method" %in% names_local) local$method <- "base"
  if (!"size_base" %in% names_local) local$size_base <- 3000
  if (!"size_new" %in% names_local) local$size_new <- 500
  if (!"reorder" %in% names_local) local$reorder <- "grts"
  if (!"kmeans" %in% names_local) {
    if (local$reorder == "none") {
      local$kmeans <- FALSE
    } else {
      local$kmeans <- TRUE
    }
  }

  if (!local$reorder %in% c("none", "random", "grts")) {
    stop("method must be \"random\", \"grts\", or \"none\".", call. = FALSE)
  }


  if (local$size_base >= n) {
    local <- list(method = "all")
  }

  if (local$method != "all") {

    if (local$size_base > 10000) {
      warning("size_base exceeds 10,000, which may result in exceedingly long computational times. Consider reducing size_base.", call. = FALSE)
    }

    if (local$size_new > 5000) {
      warning("size_new exceeds 5,000, which may result in exceedingly long computational times. Consider reducing size_new.", call. = FALSE)
    }

  }


  if (local$method != "all") {

    index <- seq(1, n)

    if (local$reorder == "random") {
      index <- sample(index)
    } else if (local$reorder == "grts") {
      if (!requireNamespace("spsurvey", quietly = TRUE)) {
        stop("Install the spsurvey package before using local method \"grts\".", call. = FALSE)
      } else {
        data_sf <- st_as_sf(data, coords = c("...xcoord...", "...ycoord..."), crs = NA)
        data_sf$...index... <- index
        samp <- spsurvey::grts(data_sf, n_base = n, projcrs_check = FALSE)
        index <- samp$sites_base$...index...
      }
    }


    index_base <- index[seq(1, local$size_base)]
    index_new <- index[-seq(1, local$size_base)]
    n_index_new <- length(index_new)
    groups <- ceiling(n_index_new / local$size_new) # consider adding groups as an argument

    if (local$kmeans) {
      # any extra local list elements beyond the recognized settings are
      # forwarded to kmeans() (e.g. nstart, algorithm), letting advanced
      # users tune the clustering without a dedicated argument for each
      kmeans_arg_names <- setdiff(names(local), c("method", "size_base", "size_new", "reorder", "kmeans", "parallel", "ncores"))
      kmeans_args <- local[kmeans_arg_names]
      x <- cbind(data[index_new, "...xcoord..."], data[index_new, "...ycoord..."])
      index_new <- split(index_new, do.call("kmeans", c(list(x = x, centers = groups, iter.max = 30), kmeans_args))$cluster)
    } else {
      # non-kmeans grouping: assign the (already ordered) remaining indices
      # to `groups` blocks of roughly equal size, distributing the
      # n_index_new %% groups leftover observations one-per-group among the
      # first few groups rather than dumping them all in the last group
      index_new <- split(index_new, rep(seq(1, groups), times = c(rep(n_index_new %/% groups + 1, n_index_new %% groups), rep(n_index_new %/% groups, groups - n_index_new %% groups))))
    }

    local$index <- list(base = index_base, new = index_new)

    if (!"parallel" %in% names_local) {
      local$parallel <- FALSE
      local$ncores <- NULL
    }

    if (local$parallel) {
      n_index <- length(unique(local$index))
      if ("ncores" %in% names_local) {
        cores_available <- parallel::detectCores()
        local$ncores <- min(n_index, local$ncores, cores_available)
      } else {
        local$ncores <- parallel::detectCores()
        local$ncores <- min(n_index, local$ncores)
      }
    }
  }


  local

}

#' Build the big data approximation settings for \code{\link{conditional}()}
#'
#' Analog of \code{\link{get_local_list_simulation}()} for conditional
#' simulation, which needs two independent big-data decisions: how to
#' subsample the *observed* data down to a base sample
#' (\code{method_base}/\code{size_base}/\code{reorder_base}), and how to
#' split the *prediction* locations into blocks
#' (\code{method_new}/\code{size_new}/\code{reorder_new}/\code{kmeans_new}).
#' Unlike \code{\link{get_local_list_simulation}()}, the base sample and
#' newdata blocks are independent of each other, so each gets its own
#' \code{method_}/\code{size_} settings.
#'
#' @param local A logical or list; see the \code{local} argument to
#'   \code{\link{conditional}()}.
#' @param object A fitted \code{splm} or \code{spglm} model object.
#' @param newdata A data frame or \code{sf} object of prediction locations.
#'
#' @return A list with the resolved \code{local} settings, always including
#'   \code{index = list(base = ..., new = ...)} (defaulting to the full index
#'   on whichever side, base or new, its \code{method_*} is \code{"all"}).
#'
#' @noRd
get_local_list_conditional <- function(local, object, newdata) {

  n <- object$n
  n_pred <- NROW(newdata)

  if (is.null(local)) {
    if (n > 5000 || n_pred > 5000) {
      local <- TRUE
      message("Because the data size or number of conditional simulations exceeds 5,000, we are setting local = TRUE to perform computationally efficient approximations. To override this behavior and compute the exact solution, rerun with local = FALSE. Be aware that setting local = FALSE may result in exceedingly long computational times.")
    } else {
      local <- FALSE
    }
  }

  if (is.logical(local)) {
    if (local) {
      local <- list()
    } else {
      local <- list(method_base = "all", method_new = "all")
    }
  }

  names_local <- names(local)

  if (!"method_base" %in% names_local) local$method_base <- "base"
  if (!"method_new" %in% names_local) local$method_new <- "base"
  if (!"size_base" %in% names_local) local$size_base <- 3000
  if (!"size_new" %in% names_local) local$size_new <- 500
  if (!"reorder_base" %in% names_local) local$reorder_base <- "grts"
  if (!"reorder_new" %in% names_local) local$reorder_new <- "random"

  if (!"kmeans_new" %in% names_local) {
    if (local$reorder_new == "none") {
      local$kmeans_new <- FALSE
    } else {
      local$kmeans_new <- TRUE
    }
  }

  if (!local$reorder_base %in% c("none", "random", "grts")) {
    stop("method must be \"grts\", \"random\", or \"none\".", call. = FALSE)
  }
  if (!local$reorder_new %in% c("none", "random")) {
    stop("method must be \"random\", or \"none\".", call. = FALSE)
  }


  if (local$size_base >= n) {
    local$method_base <- "all"
  }

  if (local$size_new >= n_pred) {
    local$method_new <- "all"
  }

  # default to the full index for whichever side (base/new) ends up not
  # needing subsetting, so local$index below is always well-formed regardless
  # of which of method_base/method_new (independently) is "all"
  index_base <- seq(1, n)
  index_new <- seq(1, n_pred)

  if (local$method_base != "all") {

    if (local$reorder_base == "random") {
      index_base <- sample(index_base)
    } else if (local$reorder_base == "grts") {
      if (!requireNamespace("spsurvey", quietly = TRUE)) {
        stop("Install the spsurvey package before using local method \"grts\".", call. = FALSE)
      } else {
        obdata_sf <- st_as_sf(object$obdata, coords = c(object$xcoord, object$ycoord), crs = NA)
        obdata_sf$.index_base <- index_base
        samp <- spsurvey::grts(obdata_sf, n_base = n, projcrs_check = FALSE)
        index_base <- samp$sites_base$.index_base
      }
    }

    index_base <- index_base[seq(1, local$size_base)]
  }

  if (local$method_new != "all") {

    if (local$reorder_new == "random") {
      index_new <- sample(index_new)
    }

    groups <- ceiling(n_pred / local$size_new) # consider adding groups as an argument

    if (local$kmeans_new) {
      kmeans_arg_names <- setdiff(names(local), c("method_base", "method_new", "size_base", "size_new", "reorder_base", "reorder_new", "kmeans_new", "parallel", "ncores"))
      kmeans_args <- local[kmeans_arg_names]

      # kmeans() needs plain x/y coordinate columns; if newdata is an sf
      # object (possibly polygons), fall back to its centroids and rename the
      # resulting geometry-derived columns to match object's coordinate names
      if (inherits(newdata, "sf")) {
        newdata <- suppressWarnings(sf::st_centroid(newdata))
        newdata <- sf_to_df(newdata)
        names(newdata)[[which(names(newdata) == ".xcoord")]] <- as.character(object$xcoord) # only relevant if newdata is sf data is not
        names(newdata)[[which(names(newdata) == ".ycoord")]] <- as.character(object$ycoord) # only relevant if newdata is sf data is not
      }
      x <- cbind(newdata[[object$xcoord]], newdata[[object$ycoord]])[index_new, ]
      index_new <- split(index_new, do.call("kmeans", c(list(x = x, centers = groups, iter.max = 30), kmeans_args))$cluster)
    } else {
      index_new <- split(index_new, rep(seq(1, groups), times = c(rep(n_pred %/% groups + 1, n_pred %% groups), rep(n_pred %/% groups, groups - n_pred %% groups))))
    }
  }

  # set unconditionally (not just when method_new != "all") so that
  # local$index$base is always available whenever method_base != "all" --
  # previously this was nested inside the method_new != "all" block, so
  # local$index was left unset (NULL) whenever a large observed data set
  # needed subsetting but a small newdata did not
  local$index <- list(base = index_base, new = index_new)

  if (!"parallel" %in% names_local) {
    local$parallel <- FALSE
    local$ncores <- NULL
  }

  if (local$parallel) {
    n_index <- length(unique(local$index))
    if ("ncores" %in% names_local) {
      cores_available <- parallel::detectCores()
      local$ncores <- min(n_index, local$ncores, cores_available)
    } else {
      local$ncores <- parallel::detectCores()
      local$ncores <- min(n_index, local$ncores)
    }
  }

  local

}
