get_local_list_estimation <- function(local, data, xcoord, ycoord, n, partition_factor) {

  # size can be an integer and sets group size
  # alternatively, set the number of groups
  # index overrides size and groups
  # set var_adjust as "none", "theoretical", "empirical", and "pooled"

  if (is.logical(local)) {
    if (local) {
      local <- list()
    } else {
      if (is.null(partition_factor)) {
        local <- list(index = rep(1, n))
      } else {
        index <- unname(model.response(model.frame(reformulate("1", response = labels(terms(partition_factor))), data = data)))
        # index <- data[[labels(terms(partition_factor))]]
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
    local$size <- NULL
    local$groups <- NULL
    local$method <- NULL
  } else {
    if (!"size" %in% names_local) {
      if ("groups" %in% names_local) {
        local$size <- ceiling(n / local$groups)
      } else {
        local$size <- 50
        local$groups <- ceiling(n / local$size)
      }
    } else {
      local$groups <- ceiling(n / local$size)
    }
    if (!"method" %in% names_local) {
      local$method <- "random"
    }
    local$index <- get_local_estimation_index(local, data, xcoord, ycoord, n)
  }

  # setting var adjust
  if (!"var_adjust" %in% names_local) {
    local$var_adjust <- "theoretical"
  } # "none", "empirical", "theoretical", and "pooled"

  # setting partition factor
  local$partition_factor <- partition_factor

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

get_local_estimation_index <- function(local, data, xcoord, ycoord, n) {
  if (local$method == "random") {
    index <- sample(rep(seq_len(local$groups), times = local$size)[seq_len(n)])
  } else if (local$method == "kmeans") {
    kmeans_args <- setdiff(names(local), c("size", "groups", "method", "index", "parallel", "ncores", "var_adjust"))
    x <- cbind(data[[xcoord]], data[[ycoord]])
    index <- do.call("kmeans", c(list(x = x, centers = local$groups), kmeans_args))$cluster
  } else {
    stop("local$method must be random (the default) or kmeans")
  }
  index
}

get_local_list_prediction <- function(local) {

  # set local neighborhood size
  # method can be "all" (for all data), "distance" (for local distance neighborhoods)
  # or "covariance" (for local covariance neighborhoods)

  if (is.logical(local)) {
    if (local) {
      local <- list(method = "covariance", size = 50)
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
    # local$method <- "all"
    local$method <- "covariance"
  }

  if (local$method %in% c("distance", "covariance") && !"size" %in% names_local) {
    local$size <- 50
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
