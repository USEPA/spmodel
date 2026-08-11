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
#' Dispatches on \code{local$approximation} to one of two big-data
#' approximations for unconditional simulation (see \code{\link{sprnorm}()}'s
#' \code{local} argument for the full description of each):
#' \itemize{
#'   \item \code{"low-rank"} (the default): \code{\link{get_local_list_simulation_lowrank}()}.
#'     A random or spatially-balanced (GRTS) ordering is used to draw a
#'     "base" sample, the remaining locations are split into blocks
#'     (optionally via k-means on coordinates), and
#'     \code{\link{get_conditional_new_from_base}()} simulates each block
#'     conditional on the base sample alone (blocks are conditionally
#'     independent given the base).
#'   \item \code{"vecchia"}: \code{\link{get_local_list_simulation_vecchia}()}.
#'     Every location is simulated sequentially, each conditional on every
#'     earlier-simulated location, optionally truncated to a nearest/most-
#'     correlated neighbor subset (\code{method}/\code{size} which is the same
#'     \code{"all"}/\code{"distance"}/\code{"covariance"} neighbor-selection
#'     convention used elsewhere in the package) (see
#'     \code{\link{get_sprnorm_vecchia}()}). There is no "base" sample to
#'     subsample at all here (contrast \code{"low-rank"}'s \code{size_base})
#'     (see \code{\link{get_local_list_conditional_vecchia}()}'s equivalent
#'     note for \code{\link{conditional}()}).
#' }
#'
#' @param local A logical or list; see the \code{local} argument to
#'   \code{\link{sprnorm}()}.
#' @param n The total number of locations to simulate.
#' @param data A data frame containing \code{...xcoord...}/\code{...ycoord...}
#'   columns.
#'
#' @return A list with the resolved \code{local} settings (shape depends on
#'   \code{approximation} (see \code{\link{get_local_list_simulation_lowrank}()}/
#'   \code{\link{get_local_list_simulation_vecchia}()})).
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
      local <- list(approximation = "low-rank", method_base = "all")
    }
  }

  names_local <- names(local)

  if (!"approximation" %in% names_local) local$approximation <- "low-rank"
  if (!local$approximation %in% c("low-rank", "vecchia")) {
    stop("local$approximation must be \"low-rank\" or \"vecchia\".", call. = FALSE)
  }

  if (local$approximation == "vecchia") {
    local <- get_local_list_simulation_vecchia(local, n, data)
  } else {
    local <- get_local_list_simulation_lowrank(local, n, data)
  }

  local
}

#' Build the \code{"low-rank"} (base+block) big data approximation settings
#' for \code{\link{sprnorm}()}
#'
#' @param local The partially-resolved \code{local} list (already has
#'   \code{approximation == "low-rank"}).
#' @param n The total number of locations to simulate.
#' @param data A data frame containing \code{...xcoord...}/\code{...ycoord...}
#'   columns (only used when \code{reorder_base = "grts"} or \code{kmeans_new = TRUE}).
#'
#' @return \code{local}, with every \code{"low-rank"} default filled in and
#'   \code{index = list(base = ..., new = ...)} set when
#'   \code{method_base != "all"}.
#'
#' @noRd
get_local_list_simulation_lowrank <- function(local, n, data) {

  names_local <- names(local)

  if (!"method_base" %in% names_local) local$method_base <- "base"
  if (!"size_base" %in% names_local) local$size_base <- 5000
  if (!"size_new" %in% names_local) local$size_new <- 1000
  if (!"reorder_base" %in% names_local) local$reorder_base <- "grts"
  if (!"kmeans_new" %in% names_local) {
    if (local$reorder_base == "none") {
      local$kmeans_new <- FALSE
    } else {
      local$kmeans_new <- TRUE
    }
  }

  if (!local$reorder_base %in% c("none", "random", "grts")) {
    stop("method must be \"random\", \"grts\", or \"none\".", call. = FALSE)
  }


  if (local$size_base >= n) {
    local <- list(approximation = "low-rank", method_base = "all")
  }

  if (local$method_base != "all") {

    if (local$size_base > 10000) {
      warning("size_base exceeds 10,000, which may result in exceedingly long computational times. Consider reducing size_base.", call. = FALSE)
    }

    if (local$size_new > 5000) {
      warning("size_new exceeds 5,000, which may result in exceedingly long computational times. Consider reducing size_new.", call. = FALSE)
    }

  }


  if (local$method_base != "all") {

    index <- seq(1, n)

    if (local$reorder_base == "random") {
      index <- sample(index)
    } else if (local$reorder_base == "grts") {
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

    if (local$kmeans_new) {
      # any extra local list elements beyond the recognized settings are
      # forwarded to kmeans() (e.g. nstart, algorithm), letting advanced
      # users tune the clustering without a dedicated argument for each
      kmeans_arg_names <- setdiff(names(local), c("approximation", "method_base", "size_base", "size_new", "reorder_base", "kmeans_new", "parallel", "ncores"))
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

    if (!"parallel" %in% names(local)) {
      local$parallel <- FALSE
      local$ncores <- NULL
    }

    if (local$parallel) {
      n_index <- length(unique(local$index))
      if ("ncores" %in% names(local)) {
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

#' Build the \code{"vecchia"} big data approximation settings for
#' \code{\link{sprnorm}()}
#'
#' Unlike \code{"low-rank"}, there is no base sample at all (every location
#' is simulated, none are treated as already known -- contrast
#' \code{\link{get_local_list_conditional_vecchia}()}, where observed data
#' *is* already known), so the only settings needed here are the simulation
#' order (over all \code{n} locations) and the per-location neighbor
#' truncation (\code{method}/\code{size}).
#'
#' @param local The partially-resolved \code{local} list (already has
#'   \code{approximation == "vecchia"}).
#' @param n The total number of locations to simulate.
#' @param data A data frame containing \code{...xcoord...}/\code{...ycoord...}
#'   columns.
#'
#' @return \code{local}, with every \code{"vecchia"} default filled in and
#'   \code{order}/\code{inv_order} set (the simulation order from
#'   \code{\link{get_decorrelate_order}()}, and its inverse).
#'
#' @noRd
get_local_list_simulation_vecchia <- function(local, n, data) {

  names_local <- names(local)

  if (!"size" %in% names_local) local$size <- 30
  if (!"method" %in% names_local) local$method <- "covariance"
  if (!local$method %in% c("all", "distance", "covariance")) {
    stop("local$method must be \"all\", \"distance\", or \"covariance\".", call. = FALSE)
  }
  if (!"ordering" %in% names_local) local$ordering <- "maxmin"
  if (!local$ordering %in% c("middleout", "outsidein", "coordinate", "maxmin", "grts", "random", "none")) {
    stop("local$ordering must be \"maxmin\", \"middleout\", \"outsidein\", \"coordinate\", \"grts\", \"random\", or \"none\".", call. = FALSE)
  }

  if ("parallel" %in% names_local && isTRUE(local$parallel)) {
    warning("local$parallel is not used when local$approximation = \"vecchia\" -- the simulation is inherently sequential (each location can depend on earlier-simulated ones), so there is no block-level work to parallelize. Ignoring.", call. = FALSE)
  }

  # method = "all" disables neighbor truncation entirely, so every location's
  # conditioning pool grows to include every earlier-simulated location,
  # unlike the rest of "vecchia" (whose per-location cost is capped by size,
  # independent of n), this makes the total cost scale roughly like n^4
  # (each of n sequential steps factors a covariance matrix up to n x n).
  # It exists to numerically verify the exactness identity against
  # local = FALSE on modest sample sizes, not as a scalable configuration.
  if (local$method == "all" && n > 2000) {
    warning("local$method = \"all\" disables neighbor truncation, so every location's conditioning set grows to include all previously-simulated locations. Unlike local$method = \"distance\"/\"covariance\", this does not scale well (cost grows roughly like n^4) and can be extremely slow for more than a few thousand locations. Consider local$method = \"distance\" or \"covariance\" instead, or reserve local$method = \"all\" for exactness checks on modest sample sizes.", call. = FALSE)
  }

  ord <- get_decorrelate_order(local$ordering, data[["...xcoord..."]], data[["...ycoord..."]])
  local$order <- ord$order
  local$inv_order <- ord$inv_order

  local
}

#' Build the big data approximation settings for \code{\link{conditional}()}
#'
#' Dispatches on \code{local$approximation} to one of two big-data
#' approximations for conditional simulation (see \code{\link{conditional}()}'s
#' \code{local} argument for the full description of each):
#' \itemize{
#'   \item \code{"low-rank"} (the default): \code{\link{get_local_list_conditional_lowrank}()}.
#'     Two independent big-data decisions: how to subsample the *observed*
#'     data down to a base sample (\code{method_base}/\code{size_base}/
#'     \code{reorder_base}), and how to split the *prediction* locations into
#'     blocks (\code{method_new}/\code{size_new}/\code{reorder_new}/
#'     \code{kmeans_new}), treated as conditionally independent given the base.
#'   \item \code{"vecchia"}: \code{\link{get_local_list_conditional_vecchia}()}.
#'     Every \code{newdata} location is simulated sequentially, each
#'     conditional on all observed data plus every earlier-simulated
#'     \code{newdata} location, optionally truncated to a nearest/most-
#'     correlated neighbor subset (\code{method}/\code{size} -- matching the
#'     \code{"all"}/\code{"distance"}/\code{"covariance"} neighbor-selection
#'     convention used elsewhere in the package, e.g. \code{predict()}'s own
#'     \code{local$method}) (see \code{\link{get_conditional_vecchia}()}).
#' }
#' \code{approximation} is a new top-level key with no analog before
#' \code{"vecchia"} existed (every big-data \code{local} list in the package
#' previously had exactly one strategy, so nothing picked between strategies);
#' it is deliberately not named \code{method} or \code{type} to avoid
#' colliding with vecchia's own \code{method} (the neighbor-selection rule,
#' same name/meaning as \code{predict()}/\code{decorrelate()}'s
#' \code{local$method}) or \code{conditional()}'s own top-level \code{type}
#' argument (the \code{spglm()} link/response/new scale).
#'
#' @param local A logical or list; see the \code{local} argument to
#'   \code{\link{conditional}()}.
#' @param object A fitted \code{splm} or \code{spglm} model object.
#' @param newdata A data frame or \code{sf} object of prediction locations.
#'
#' @return A list with the resolved \code{local} settings (shape depends on
#'   \code{approximation}; see \code{\link{get_local_list_conditional_lowrank}()}/
#'   \code{\link{get_local_list_conditional_vecchia}()}).
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
      local <- list(approximation = "low-rank", method_base = "all", method_new = "all")
    }
  }

  names_local <- names(local)

  if (!"approximation" %in% names_local) local$approximation <- "low-rank"
  if (!local$approximation %in% c("low-rank", "vecchia")) {
    stop("local$approximation must be \"low-rank\" or \"vecchia\".", call. = FALSE)
  }

  if (local$approximation == "vecchia") {
    local <- get_local_list_conditional_vecchia(local, object, newdata)
  } else {
    local <- get_local_list_conditional_lowrank(local, object, newdata, n, n_pred)
  }

  if (!"parallel" %in% names(local)) {
    local$parallel <- FALSE
    local$ncores <- NULL
  }

  # vecchia's simulation loop is inherently sequential (each newdata location
  # can depend on earlier-simulated ones), so there is no block-level index
  # to parallelize over the way "low-rank" has. get_local_list_conditional_vecchia()
  # already warns if the user set parallel = TRUE, so this is just a guard
  if (local$approximation == "low-rank" && local$parallel) {
    n_index <- length(unique(local$index))
    if ("ncores" %in% names(local)) {
      cores_available <- parallel::detectCores()
      local$ncores <- min(n_index, local$ncores, cores_available)
    } else {
      local$ncores <- parallel::detectCores()
      local$ncores <- min(n_index, local$ncores)
    }
  }

  local

}

#' Build the \code{"low-rank"} (base+block) big data approximation settings
#' for \code{\link{conditional}()}
#'
#' The base-sample settings (\code{method_base}/\code{size_base}/
#' \code{reorder_base}) and the \code{newdata}-blocking settings
#' (\code{method_new}/\code{size_new}/\code{reorder_new}/\code{kmeans_new})
#' are independent of one another, so each gets its own \code{method_}/
#' \code{size_} settings.
#'
#' @param local The partially-resolved \code{local} list (already has
#'   \code{approximation == "low-rank"}).
#' @param object A fitted \code{splm} or \code{spglm} model object.
#' @param newdata A data frame or \code{sf} object of prediction locations.
#' @param n The observed sample size.
#' @param n_pred The number of \code{newdata} rows.
#'
#' @return \code{local}, with every \code{"low-rank"} default filled in and
#'   \code{index = list(base = ..., new = ...)} set (defaulting to the full
#'   index on whichever side, base or new, its \code{method_*} is \code{"all"}).
#'
#' @noRd
get_local_list_conditional_lowrank <- function(local, object, newdata, n, n_pred) {

  names_local <- names(local)

  if (!"method_base" %in% names_local) local$method_base <- "base"
  if (!"method_new" %in% names_local) local$method_new <- "base"
  if (!"size_base" %in% names_local) local$size_base <- 5000
  if (!"size_new" %in% names_local) local$size_new <- 1000
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
      kmeans_arg_names <- setdiff(names(local), c("approximation", "method_base", "method_new", "size_base", "size_new", "reorder_base", "reorder_new", "kmeans_new", "parallel", "ncores"))
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

  local
}

#' Build the \code{"vecchia"} big data approximation settings for
#' \code{\link{conditional}()}
#'
#' Unlike \code{"low-rank"}, the observed data is never subsampled (see
#' \code{\link{get_conditional_vecchia}()}), so the only settings needed here
#' are the \code{newdata} simulation order and the per-location neighbor
#' truncation (\code{method}/\code{size} -- \code{"all"}/\code{"distance"}/
#' \code{"covariance"}, the same neighbor-selection convention \code{predict()}/
#' \code{decorrelate()} already use for their own \code{local$method}).
#'
#' @param local The partially-resolved \code{local} list (already has
#'   \code{approximation == "vecchia"}).
#' @param object A fitted \code{splm} or \code{spglm} model object.
#' @param newdata A data frame or \code{sf} object of prediction locations.
#'
#' @return \code{local}, with every \code{"vecchia"} default filled in and
#'   \code{order}/\code{inv_order} set (the \code{newdata} simulation order
#'   from \code{\link{get_decorrelate_order}()}, and its inverse).
#'
#' @noRd
get_local_list_conditional_vecchia <- function(local, object, newdata) {

  names_local <- names(local)

  if (!"size" %in% names_local) local$size <- 30
  if (!"method" %in% names_local) local$method <- "covariance"
  if (!local$method %in% c("all", "distance", "covariance")) {
    stop("local$method must be \"all\", \"distance\", or \"covariance\".", call. = FALSE)
  }
  if (!"ordering" %in% names_local) local$ordering <- "maxmin"
  if (!local$ordering %in% c("middleout", "outsidein", "coordinate", "maxmin", "grts", "random", "none")) {
    stop("local$ordering must be \"maxmin\", \"middleout\", \"outsidein\", \"coordinate\", \"grts\", \"random\", or \"none\".", call. = FALSE)
  }

  if ("parallel" %in% names_local && isTRUE(local$parallel)) {
    warning("local$parallel is not used when local$approximation = \"vecchia\". Ignoring.", call. = FALSE)
  }

  # method = "all" disables neighbor truncation entirely, so every newdata
  # location's conditioning pool always includes ALL observed data (never
  # subsampled for "vecchia") plus every earlier-simulated newdata location --
  # unlike the rest of "vecchia" (whose per-location cost is capped by size,
  # independent of the observed sample size), this makes the total cost scale
  # roughly like n_pred * n_obs^3 (each of n_pred sequential steps factors a
  # covariance matrix close to n_obs x n_obs in size). It exists to
  # numerically verify the exactness identity against local = FALSE on modest
  # sample sizes, not as a scalable configuration.
  if (local$method == "all" && object$n > 2000) {
    warning("local$method = \"all\" should not be used with \"vecchia\" for large sample sizes because of exceedingly long computational times.", call. = FALSE)
  }

  # newdata coordinates for ordering only -- sf objects fall back to
  # centroids, matching the "low-rank" path's kmeans_new handling
  if (inherits(newdata, "sf")) {
    newdata <- suppressWarnings(sf::st_centroid(newdata))
    newdata <- sf_to_df(newdata)
    names(newdata)[[which(names(newdata) == ".xcoord")]] <- as.character(object$xcoord)
    names(newdata)[[which(names(newdata) == ".ycoord")]] <- as.character(object$ycoord)
  }
  xcoord_new <- newdata[[object$xcoord]]
  ycoord_new <- newdata[[object$ycoord]]

  ord <- get_decorrelate_order(local$ordering, xcoord_new, ycoord_new)
  local$order <- ord$order
  local$inv_order <- ord$inv_order

  local
}
