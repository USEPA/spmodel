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
    kmeans_args <- setdiff(names(local), c("size", "groups", "method", "index", "parallel", "ncores", "var_adjust"))
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
