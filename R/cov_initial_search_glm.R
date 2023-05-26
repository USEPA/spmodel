#' Search for initial covariance parameters
#'
#' @param spcov_initial_NA A spatial initial NA object
#' @param dispersion_initial_NA A dispersion initial NA object
#' @param ... Additional arguments passed to other methods
#' @noRd
cov_initial_search_glm <- function(spcov_initial_NA, dispersion_initial_NA, ...) {
  UseMethod("cov_initial_search_glm", spcov_initial_NA)
}
#' @export
cov_initial_search_glm.exponential <- function(spcov_initial_NA, dispersion_initial_NA, estmethod, data_object,
                                               dist_matrix_list, weights,
                                               randcov_initial_NA = NULL, esv_dotlist, ...) {


  # find ols sample variance
  X <- do.call("rbind", data_object$X_list)
  y <- do.call("rbind", data_object$y_list)
  # # same as
  # betahat <- as.vector(solve(R_val) %*% t(qr.Q(qr_val)) %*% y_trans)
  # # same as
  # betahat <- solve(t(X) %*% X) %*% t(X) %*% y_trans
  s2 <- data_object$s2
  ns2 <- 1.2 * s2
  # # old way
  # lm_out <- lm(data_object$formula, do.call("rbind", data_object$obdata_list))
  # if (inherits(dispersion_initial_NA, "binomial") && length(summary(lm_out)) == 2) {
  #   s2 <- summary(lm_out)[[1]]$sigma^2
  # } else {
  #   s2 <- summary(lm_out)$sigma^2
  # }
  # ns2 <- log(1 + 1.2 * s2) # so no negative variances

  # find sets of starting values
  de <- c(0.1, 0.5, 0.9)
  # ie
  ie <- c(0.1, 0.5, 0.9)
  ## range
  range <- get_initial_range(class(spcov_initial_NA), data_object$max_halfdist) * c(0.5, 1.5)
  if (data_object$anisotropy) {
    ## rotate
    rotate <- c(0, 30 * pi / 180, 60 * pi / 180)
    ## scale
    scale <- c(0.25, 0.75, 1)
  } else {
    rotate <- 0
    scale <- 1
  }


  # find starting spatial grid
  spcov_grid <- expand.grid(de = de, ie = ie, range = range, rotate = rotate, scale = scale)
  spcov_grid <- spcov_grid[spcov_grid$de + spcov_grid$ie == 1, , drop = FALSE]
  spcov_grid[, c("de", "ie")] <- ns2 * spcov_grid[, c("de", "ie")]

  # save initial state (used with random effects)
  spcov_grid_init <- spcov_grid

  # replace with initial values
  for (x in names(spcov_grid)) {
    if (!is.na(spcov_initial_NA$initial[[x]])) {
      spcov_grid[, x] <- spcov_initial_NA$initial[[x]]
    }
  }


  # take unique rows
  spcov_grid <- unique(spcov_grid)

  # perform search if no random effects
  if (is.null(randcov_initial_NA)) {

    # add dispersion
    spcov_grid$dispersion <- 1
    if (!is.na(dispersion_initial_NA$initial)) {
      spcov_grid[, "dispersion"] <- dispersion_initial_NA$initial
    }


    # split
    cov_grid_splits <- split(spcov_grid, seq_len(NROW(spcov_grid)))
    # apply (add family)
    objvals <- vapply(
      X = cov_grid_splits, FUN = eval_grid_glm, FUN.VALUE = numeric(1),
      data_object = data_object, spcov_type = class(spcov_initial_NA),
      family = class(dispersion_initial_NA), estmethod = estmethod, dist_matrix_list = dist_matrix_list
    )

    # find minimum value and record parameters
    min_params <- unlist(cov_grid_splits[[which.min(objvals)]])
    spcov_params <- min_params[c("de", "ie", "range", "rotate", "scale")]
    spcov_initial_NA$initial <- spcov_params
    dispersion_params <- min_params[c("dispersion")]
    dispersion_initial_NA$initial <- dispersion_params

    # return the best parameters
    best_params <- list(
      spcov_initial_val = spcov_initial_NA,
      dispersion_initial_val = dispersion_initial_NA, randcov_initial_val = NULL
    )
  } else {

    # randcov names
    randcov_names <- data_object$randcov_names
    # find number of random effects
    nvar_randcov <- length(randcov_names)

    # spatially dominant grid
    ## spatial components 90% of variance
    spcov_grid[, c("de", "ie")] <- 0.9 * spcov_grid[, c("de", "ie")]

    # replace with initial
    for (x in names(spcov_grid)) {
      if (!is.na(spcov_initial_NA$initial[[x]])) {
        spcov_grid[, x] <- spcov_initial_NA$initial[[x]]
      }
    }
    ## random effects 10% of variance and evenly spread
    for (x in randcov_names) {
      if (!is.na(randcov_initial_NA$initial[[x]])) {
        spcov_grid[, x] <- randcov_initial_NA$initial[[x]]
      } else {
        spcov_grid[, x] <- 0.1 * ns2 / nvar_randcov
      }
    }

    # Evenly dominated grid
    ## find number of spatial variance parameters
    nvar_spcov <- 2
    ## find all variance parameters
    nvar_cov <- nvar_spcov + nvar_randcov
    ## keep only spatial cases with variance spread out
    evencov_grid <- spcov_grid_init[spcov_grid_init$de == spcov_grid_init$ie, , drop = FALSE]
    ## scale to incorporate overall variance
    evencov_grid[, c("de", "ie")] <- nvar_spcov / nvar_cov * evencov_grid[, c("de", "ie")]

    # replace with initial values
    ## spatial
    for (x in names(evencov_grid)) {
      if (!is.na(spcov_initial_NA$initial[[x]])) {
        evencov_grid[, x] <- spcov_initial_NA$initial[[x]]
      }
    }
    ## random
    for (x in randcov_names) {
      if (!is.na(randcov_initial_NA$initial[[x]])) {
        evencov_grid[, x] <- randcov_initial_NA$initial[[x]]
      } else {
        evencov_grid[, x] <- 1 / nvar_cov * ns2
      }
    }
    # find unique combinations
    evencov_grid <- unique(evencov_grid)

    # random dominated grid
    ## spatial component 10%
    randcov_grid_spcov <- spcov_grid_init[spcov_grid_init$de == spcov_grid_init$ie & spcov_grid_init$range == min(spcov_grid_init$range), , drop = FALSE]
    randcov_grid_spcov[, c("de", "ie")] <- 0.1 * randcov_grid_spcov[, c("de", "ie")]

    # replace spatial initial values
    for (x in names(randcov_grid_spcov)) {
      if (!is.na(spcov_initial_NA$initial[[x]])) {
        randcov_grid_spcov[, x] <- spcov_initial_NA$initial[[x]]
      }
    }
    # find unique values
    randcov_grid_spcov <- unique(randcov_grid_spcov)

    # random dominant grid
    randcov_grid_randcov <- as.data.frame(as.list(rep(1 / nvar_randcov, nvar_randcov)))
    names(randcov_grid_randcov) <- randcov_names
    ## if there is more than one random effect, split it up into relevant scenarios
    if (nvar_randcov > 1) {
      ## set 0.1 for all proportions
      extra_grid_randcov <- lapply(seq_len(nvar_randcov), function(x) rep(0.1 * 1 / (nvar_randcov - 1), nvar_randcov))
      extra_grid_randcov <- do.call("rbind", extra_grid_randcov)
      ## fill in 0.9 for one proportion in each row
      diag(extra_grid_randcov) <- 0.9
      extra_grid_randcov <- as.data.frame(extra_grid_randcov)
      names(extra_grid_randcov) <- randcov_names
    } else {
      extra_grid_randcov <- NULL
    }
    ## give the random effects 90% of the variance
    randcov_grid_randcov <- 0.9 * ns2 * rbind(randcov_grid_randcov, extra_grid_randcov)
    for (x in randcov_names) {
      if (!is.na(randcov_initial_NA$initial[[x]])) {
        randcov_grid_randcov[, x] <- randcov_initial_NA$initial[[x]]
      }
    }

    ## bind together and replicate
    randcov_grid_spcov_rep <- do.call("rbind", replicate(NROW(randcov_grid_randcov), randcov_grid_spcov, simplify = FALSE))
    randcov_grid_randcov_rep <- do.call("rbind", replicate(NROW(randcov_grid_spcov), randcov_grid_randcov, simplify = FALSE))
    randcov_grid <- cbind(randcov_grid_spcov_rep, randcov_grid_randcov_rep)

    # bind together all grids
    cov_grid <- rbind(spcov_grid, evencov_grid, randcov_grid)
    cov_grid <- unique(cov_grid)

    # add dispersion
    cov_grid$dispersion <- 1
    if (!is.na(dispersion_initial_NA$initial)) {
      cov_grid[, "dispersion"] <- dispersion_initial_NA$initial
    }

    cov_grid_splits <- split(cov_grid, seq_len(NROW(cov_grid)))
    # apply
    objvals <- vapply(
      X = cov_grid_splits, FUN = eval_grid_glm, FUN.VALUE = numeric(1),
      data_object = data_object, spcov_type = class(spcov_initial_NA),
      family = class(dispersion_initial_NA), estmethod = estmethod, dist_matrix_list = dist_matrix_list
    )
    # find minimum value and record parameters
    min_params <- unlist(cov_grid_splits[[which.min(objvals)]])
    spcov_params <- min_params[c("de", "ie", "range", "rotate", "scale")]
    # return the spatial parameters
    spcov_initial_NA$initial <- spcov_params
    # return the dispersion parameter
    dispersion_params <- min_params[c("dispersion")]
    dispersion_initial_NA$initial <- dispersion_params
    # return the random effect parameters
    randcov_initial_NA$initial <- randcov_params(min_params[randcov_names])
    # return the best parameters
    best_params <- list(
      spcov_initial_val = spcov_initial_NA,
      dispersion_initial_val = dispersion_initial_NA, randcov_initial_val = randcov_initial_NA
    )
  }

  # return the best parameters
  best_params
}

#' @export
cov_initial_search_glm.spherical <- cov_initial_search_glm.exponential
#' @export
cov_initial_search_glm.gaussian <- cov_initial_search_glm.exponential
#' @export
cov_initial_search_glm.triangular <- cov_initial_search_glm.exponential
#' @export
cov_initial_search_glm.circular <- cov_initial_search_glm.exponential
#' @export
cov_initial_search_glm.cubic <- cov_initial_search_glm.exponential
#' @export
cov_initial_search_glm.pentaspherical <- cov_initial_search_glm.exponential
#' @export
cov_initial_search_glm.cosine <- cov_initial_search_glm.exponential
#' @export
cov_initial_search_glm.wave <- cov_initial_search_glm.exponential
#' @export
cov_initial_search_glm.jbessel <- cov_initial_search_glm.exponential
#' @export
cov_initial_search_glm.gravity <- cov_initial_search_glm.exponential
#' @export
cov_initial_search_glm.rquad <- cov_initial_search_glm.exponential
#' @export
cov_initial_search_glm.magnetic <- cov_initial_search_glm.exponential

#' @export
cov_initial_search_glm.none <- function(spcov_initial_NA, dispersion_initial_NA, estmethod, data_object,
                                        dist_matrix_list, weights,
                                        randcov_initial_NA = NULL, esv_dotlist, ...) {
  # find ols sample variance
  X <- do.call("rbind", data_object$X_list)
  y <- do.call("rbind", data_object$y_list)
  s2 <- data_object$s2
  ns2 <- 1.2 * s2

  # find sets of starting values
  ## de
  de <- 0
  ## ie
  ie <- 1
  ## range
  range <- get_initial_range(class(spcov_initial_NA), NULL)
  ## rotate
  rotate <- 0
  ## scale
  scale <- 1


  # find starting spatial grid
  spcov_grid <- expand.grid(de = de, ie = ie, range = range, rotate = rotate, scale = scale)
  spcov_grid <- spcov_grid[spcov_grid$de + spcov_grid$ie == 1, , drop = FALSE]
  spcov_grid[, c("de", "ie")] <- ns2 * spcov_grid[, c("de", "ie")]

  # save initial state (used with random effects)
  spcov_grid_init <- spcov_grid

  # replace with initial values
  for (x in names(spcov_grid)) {
    if (!is.na(spcov_initial_NA$initial[[x]])) {
      spcov_grid[, x] <- spcov_initial_NA$initial[[x]]
    }
  }

  # take unique rows
  spcov_grid <- unique(spcov_grid)

  # perform search if no random effects
  if (is.null(randcov_initial_NA)) {
    # add dispersion
    spcov_grid$dispersion <- 1
    if (!is.na(dispersion_initial_NA$initial)) {
      spcov_grid[, "dispersion"] <- dispersion_initial_NA$initial
    }

    # split
    cov_grid_splits <- split(spcov_grid, seq_len(NROW(spcov_grid)))
    # apply (add family)
    objvals <- vapply(
      X = cov_grid_splits, FUN = eval_grid_glm, FUN.VALUE = numeric(1),
      data_object = data_object, spcov_type = class(spcov_initial_NA),
      family = class(dispersion_initial_NA), estmethod = estmethod, dist_matrix_list = dist_matrix_list
    )

    # find minimum value and record parameters
    min_params <- unlist(cov_grid_splits[[which.min(objvals)]])
    spcov_params <- min_params[c("de", "ie", "range", "rotate", "scale")]
    spcov_initial_NA$initial <- spcov_params
    dispersion_params <- min_params[c("dispersion")]
    dispersion_initial_NA$initial <- dispersion_params

    # return the best parameters
    best_params <- list(
      spcov_initial_val = spcov_initial_NA,
      dispersion_initial_val = dispersion_initial_NA, randcov_initial_val = NULL
    )
  } else {

    # randcov vars names
    randcov_names <- data_object$randcov_names
    # find number of random effects
    nvar_randcov <- length(randcov_names)

    # spatially dominant grid
    ## spatial components 90% of variance
    spcov_grid[, c("de", "ie")] <- 0.9 * spcov_grid[, c("de", "ie")]

    # replace with initial
    for (x in names(spcov_grid)) {
      if (!is.na(spcov_initial_NA$initial[[x]])) {
        spcov_grid[, x] <- spcov_initial_NA$initial[[x]]
      }
    }

    ## random effects 10% of variance and evenly spread
    for (x in randcov_names) {
      if (!is.na(randcov_initial_NA$initial[[x]])) {
        spcov_grid[, x] <- randcov_initial_NA$initial[[x]]
      } else {
        spcov_grid[, x] <- 0.1 * ns2 / nvar_randcov
      }
    }

    # Evenly dominated grid
    ## find number of spatial variance parameters
    nvar_spcov <- 1
    ## find all variance parameters
    nvar_cov <- nvar_spcov + nvar_randcov
    ## keep only spatial cases with variance spread out
    evencov_grid <- spcov_grid_init
    ## scale to incorporate overall variance
    evencov_grid[, c("de", "ie")] <- nvar_spcov / nvar_cov * evencov_grid[, c("de", "ie")]

    # replace with initial values
    ## spatial
    for (x in names(evencov_grid)) {
      if (!is.na(spcov_initial_NA$initial[[x]])) {
        evencov_grid[, x] <- spcov_initial_NA$initial[[x]]
      }
    }
    ## random
    for (x in randcov_names) {
      if (!is.na(randcov_initial_NA$initial[[x]])) {
        evencov_grid[, x] <- randcov_initial_NA$initial[[x]]
      } else {
        evencov_grid[, x] <- 1 / nvar_cov * ns2
      }
    }
    # find unique combinations
    evencov_grid <- unique(evencov_grid)

    # random dominated grid
    ## spatial component 10%
    randcov_grid_spcov <- spcov_grid_init
    randcov_grid_spcov[, c("de", "ie")] <- 0.1 * randcov_grid_spcov[, c("de", "ie")]

    # replace spatial initial values
    for (x in names(randcov_grid_spcov)) {
      if (!is.na(spcov_initial_NA$initial[[x]])) {
        randcov_grid_spcov[, x] <- spcov_initial_NA$initial[[x]]
      }
    }
    # find unique values
    randcov_grid_spcov <- unique(randcov_grid_spcov)

    # random dominant grid
    randcov_grid_randcov <- as.data.frame(as.list(rep(1 / nvar_randcov, nvar_randcov)))
    names(randcov_grid_randcov) <- randcov_names
    ## if there is more than one random effect, split it up into relevant scenarios
    if (nvar_randcov > 1) {
      ## set 0.1 for all proportions
      extra_grid_randcov <- lapply(seq_len(nvar_randcov), function(x) rep(0.1 * 1 / (nvar_randcov - 1), nvar_randcov))
      extra_grid_randcov <- do.call("rbind", extra_grid_randcov)
      ## fill in 0.9 for one proportion in each row
      diag(extra_grid_randcov) <- 0.9
      extra_grid_randcov <- as.data.frame(extra_grid_randcov)
      names(extra_grid_randcov) <- randcov_names
    } else {
      extra_grid_randcov <- NULL
    }
    ## give the random effects 90% of the variance
    randcov_grid_randcov <- 0.9 * ns2 * rbind(randcov_grid_randcov, extra_grid_randcov)
    for (x in randcov_names) {
      if (!is.na(randcov_initial_NA$initial[[x]])) {
        randcov_grid_randcov[, x] <- randcov_initial_NA$initial[[x]]
      }
    }

    ## bind together and replicate
    randcov_grid_spcov_rep <- do.call("rbind", replicate(NROW(randcov_grid_randcov), randcov_grid_spcov, simplify = FALSE))
    randcov_grid_randcov_rep <- do.call("rbind", replicate(NROW(randcov_grid_spcov), randcov_grid_randcov, simplify = FALSE))
    randcov_grid <- cbind(randcov_grid_spcov_rep, randcov_grid_randcov_rep)

    # bind together all grids
    cov_grid <- rbind(spcov_grid, evencov_grid, randcov_grid)
    cov_grid <- unique(cov_grid)

    # add dispersion
    cov_grid$dispersion <- 1
    if (!is.na(dispersion_initial_NA$initial)) {
      cov_grid[, "dispersion"] <- dispersion_initial_NA$initial
    }

    cov_grid_splits <- split(cov_grid, seq_len(NROW(cov_grid)))
    # apply
    objvals <- vapply(
      X = cov_grid_splits, FUN = eval_grid_glm, FUN.VALUE = numeric(1),
      data_object = data_object, spcov_type = class(spcov_initial_NA),
      family = class(dispersion_initial_NA), estmethod = estmethod, dist_matrix_list = dist_matrix_list
    )
    # find minimum value and record parameters
    min_params <- unlist(cov_grid_splits[[which.min(objvals)]])
    spcov_params <- min_params[c("de", "ie", "range", "rotate", "scale")]
    # return the spatial parameters
    spcov_initial_NA$initial <- spcov_params
    # return the dispersion parameter
    dispersion_params <- min_params[c("dispersion")]
    dispersion_initial_NA$initial <- dispersion_params
    # return the random effect parameters
    randcov_initial_NA$initial <- randcov_params(min_params[randcov_names])
    # return the best parameters
    best_params <- list(
      spcov_initial_val = spcov_initial_NA,
      dispersion_initial_val = dispersion_initial_NA, randcov_initial_val = randcov_initial_NA
    )
  }
  # return the best parameters
  best_params
}

#' @export
cov_initial_search_glm.matern <- function(spcov_initial_NA, dispersion_initial_NA, estmethod, data_object,
                                          dist_matrix_list, weights,
                                          randcov_initial_NA = NULL, esv_dotlist, ...) {

  # find ols sample variance
  X <- do.call("rbind", data_object$X_list)
  y <- do.call("rbind", data_object$y_list)
  s2 <- data_object$s2
  ns2 <- 1.2 * s2



  # find sets of starting values
  de <- c(0.1, 0.5, 0.9)
  # ie
  ie <- c(0.1, 0.5, 0.9)
  ## range
  range <- get_initial_range(class(spcov_initial_NA), data_object$max_halfdist) * c(0.5, 1.5)
  ## extra
  extra <- get_initial_extra(class(spcov_initial_NA)) * c(0.5, 2)
  if (data_object$anisotropy) {
    ## rotate
    rotate <- c(0, 30 * pi / 180, 60 * pi / 180)
    ## scale
    scale <- c(0.25, 0.75, 1)
  } else {
    rotate <- 0
    scale <- 1
  }


  # find starting spatial grid
  spcov_grid <- expand.grid(de = de, ie = ie, range = range, extra = extra, rotate = rotate, scale = scale)
  spcov_grid <- spcov_grid[spcov_grid$de + spcov_grid$ie == 1, , drop = FALSE]
  spcov_grid[, c("de", "ie")] <- ns2 * spcov_grid[, c("de", "ie")]

  # save initial state (used with random effects)
  spcov_grid_init <- spcov_grid

  # replace with initial values
  for (x in names(spcov_grid)) {
    if (!is.na(spcov_initial_NA$initial[[x]])) {
      spcov_grid[, x] <- spcov_initial_NA$initial[[x]]
    }
  }

  # take unique rows
  spcov_grid <- unique(spcov_grid)

  # perform search if no random effects
  if (is.null(randcov_initial_NA)) {
    # add dispersion
    spcov_grid$dispersion <- 1
    if (!is.na(dispersion_initial_NA$initial)) {
      spcov_grid[, "dispersion"] <- dispersion_initial_NA$initial
    }

    # split
    cov_grid_splits <- split(spcov_grid, seq_len(NROW(spcov_grid)))
    # apply (add family)
    objvals <- vapply(
      X = cov_grid_splits, FUN = eval_grid_glm, FUN.VALUE = numeric(1),
      data_object = data_object, spcov_type = class(spcov_initial_NA),
      family = class(dispersion_initial_NA), estmethod = estmethod, dist_matrix_list = dist_matrix_list
    )

    # find minimum value and record parameters
    min_params <- unlist(cov_grid_splits[[which.min(objvals)]])
    spcov_params <- min_params[c("de", "ie", "range", "extra", "rotate", "scale")]
    spcov_initial_NA$initial <- spcov_params
    dispersion_params <- min_params[c("dispersion")]
    dispersion_initial_NA$initial <- dispersion_params

    # return the best parameters
    best_params <- list(
      spcov_initial_val = spcov_initial_NA,
      dispersion_initial_val = dispersion_initial_NA, randcov_initial_val = NULL
    )
  } else {
    # randcov names
    randcov_names <- data_object$randcov_names
    # find number of random effects
    nvar_randcov <- length(randcov_names)

    # spatially dominant grid
    ## spatial components 90% of variance
    spcov_grid[, c("de", "ie")] <- 0.9 * spcov_grid[, c("de", "ie")]

    # replace with initial
    for (x in names(spcov_grid)) {
      if (!is.na(spcov_initial_NA$initial[[x]])) {
        spcov_grid[, x] <- spcov_initial_NA$initial[[x]]
      }
    }

    ## random effects 10% of variance and evenly spread
    for (x in randcov_names) {
      if (!is.na(randcov_initial_NA$initial[[x]])) {
        spcov_grid[, x] <- randcov_initial_NA$initial[[x]]
      } else {
        spcov_grid[, x] <- 0.1 * ns2 / nvar_randcov
      }
    }



    # Evenly dominated grid
    ## find number of spatial variance parameters
    nvar_spcov <- 2
    ## find all variance parameters
    nvar_cov <- nvar_spcov + nvar_randcov
    ## keep only spatial cases with variance spread out
    evencov_grid <- spcov_grid_init[spcov_grid_init$de == spcov_grid_init$ie, , drop = FALSE]
    ## scale to incorporate overall variance
    evencov_grid[, c("de", "ie")] <- nvar_spcov / nvar_cov * evencov_grid[, c("de", "ie")]

    # replace with initial values
    ## spatial
    for (x in names(evencov_grid)) {
      if (!is.na(spcov_initial_NA$initial[[x]])) {
        evencov_grid[, x] <- spcov_initial_NA$initial[[x]]
      }
    }
    ## random
    for (x in randcov_names) {
      if (!is.na(randcov_initial_NA$initial[[x]])) {
        evencov_grid[, x] <- randcov_initial_NA$initial[[x]]
      } else {
        evencov_grid[, x] <- 1 / nvar_cov * ns2
      }
    }
    # find unique combinations
    evencov_grid <- unique(evencov_grid)

    # random dominated grid
    ## spatial component 10%
    randcov_grid_spcov <- spcov_grid_init[spcov_grid_init$de == spcov_grid_init$ie & spcov_grid_init$range == min(spcov_grid_init$range) & spcov_grid_init$extra == min(spcov_grid_init$extra), , drop = FALSE]
    randcov_grid_spcov[, c("de", "ie")] <- 0.1 * randcov_grid_spcov[, c("de", "ie")]

    # replace spatial initial values
    for (x in names(randcov_grid_spcov)) {
      if (!is.na(spcov_initial_NA$initial[[x]])) {
        randcov_grid_spcov[, x] <- spcov_initial_NA$initial[[x]]
      }
    }
    # find unique values
    randcov_grid_spcov <- unique(randcov_grid_spcov)

    # random dominant grid
    randcov_grid_randcov <- as.data.frame(as.list(rep(1 / nvar_randcov, nvar_randcov)))
    names(randcov_grid_randcov) <- randcov_names
    ## if there is more than one random effect, split it up into relevant scenarios
    if (nvar_randcov > 1) {
      ## set 0.1 for all proportions
      extra_grid_randcov <- lapply(seq_len(nvar_randcov), function(x) rep(0.1 * 1 / (nvar_randcov - 1), nvar_randcov))
      extra_grid_randcov <- do.call("rbind", extra_grid_randcov)
      ## fill in 0.9 for one proportion in each row
      diag(extra_grid_randcov) <- 0.9
      extra_grid_randcov <- as.data.frame(extra_grid_randcov)
      names(extra_grid_randcov) <- randcov_names
    } else {
      extra_grid_randcov <- NULL
    }
    ## give the random effects 90% of the variance
    randcov_grid_randcov <- 0.9 * ns2 * rbind(randcov_grid_randcov, extra_grid_randcov)
    for (x in randcov_names) {
      if (!is.na(randcov_initial_NA$initial[[x]])) {
        randcov_grid_randcov[, x] <- randcov_initial_NA$initial[[x]]
      }
    }

    ## bind together and replicate
    randcov_grid_spcov_rep <- do.call("rbind", replicate(NROW(randcov_grid_randcov), randcov_grid_spcov, simplify = FALSE))
    randcov_grid_randcov_rep <- do.call("rbind", replicate(NROW(randcov_grid_spcov), randcov_grid_randcov, simplify = FALSE))
    randcov_grid <- cbind(randcov_grid_spcov_rep, randcov_grid_randcov_rep)

    # bind together all grids
    cov_grid <- rbind(spcov_grid, evencov_grid, randcov_grid)
    cov_grid <- unique(cov_grid)
    # add dispersion
    cov_grid$dispersion <- 1
    if (!is.na(dispersion_initial_NA$initial)) {
      cov_grid[, "dispersion"] <- dispersion_initial_NA$initial
    }

    cov_grid_splits <- split(cov_grid, seq_len(NROW(cov_grid)))
    # apply
    objvals <- vapply(
      X = cov_grid_splits, FUN = eval_grid_glm, FUN.VALUE = numeric(1),
      data_object = data_object, spcov_type = class(spcov_initial_NA),
      family = class(dispersion_initial_NA), estmethod = estmethod, dist_matrix_list = dist_matrix_list
    )
    # find minimum value and record parameters
    min_params <- unlist(cov_grid_splits[[which.min(objvals)]])
    spcov_params <- min_params[c("de", "ie", "range", "extra", "rotate", "scale")]
    # return the spatial parameters
    spcov_initial_NA$initial <- spcov_params
    # return the dispersion parameter
    dispersion_params <- min_params[c("dispersion")]
    dispersion_initial_NA$initial <- dispersion_params
    # return the random effect parameters
    randcov_initial_NA$initial <- randcov_params(min_params[randcov_names])
    # return the best parameters
    best_params <- list(
      spcov_initial_val = spcov_initial_NA,
      dispersion_initial_val = dispersion_initial_NA, randcov_initial_val = randcov_initial_NA
    )
  }
  # return the best parameters
  best_params
}

#' @export
cov_initial_search_glm.cauchy <- cov_initial_search_glm.matern
#' @export
cov_initial_search_glm.pexponential <- cov_initial_search_glm.matern

#' @export
cov_initial_search_glm.car <- function(spcov_initial_NA, dispersion_initial_NA, estmethod, data_object,
                                       dist_matrix_list, randcov_initial_NA = NULL, ...) {

  # find ols sample variance
  X <- data_object$X
  y <- data_object$y
  s2 <- data_object$s2
  ns2 <- 1.2 * s2

  # store W as dist_matrix
  # MAKE IT CLEAR DIST_MATRIX_LIST IS NOT A LIST
  W <- dist_matrix_list

  # find sets of starting values
  ## de
  de <- c(0.1, 0.5, 0.9)
  ## ie
  ie <- c(0.1, 0.5, 0.9)
  ## range
  rho_length <- data_object$rho_ub - data_object$rho_lb
  range <- c(data_object$rho_lb + 0.25 * rho_length, data_object$rho_ub - 0.25 * rho_length)

  # find starting spatial grid
  spcov_grid <- expand.grid(de = de, ie = ie, range = range)
  spcov_grid <- spcov_grid[spcov_grid$de + spcov_grid$ie == 1, , drop = FALSE]
  spcov_grid[, c("de", "ie")] <- ns2 * spcov_grid[, c("de", "ie")]
  spcov_grid$extra <- spcov_grid$de

  # save initial state (used with random effects)
  spcov_grid_init <- spcov_grid

  # replace with initial values
  for (x in names(spcov_grid)) {
    if (!is.na(spcov_initial_NA$initial[[x]])) {
      spcov_grid[, x] <- spcov_initial_NA$initial[[x]]
    }
  }

  # take unique rows
  spcov_grid <- unique(spcov_grid)

  # perform search if no random effects
  if (is.null(randcov_initial_NA)) {
    # add dispersion
    spcov_grid$dispersion <- 1
    if (!is.na(dispersion_initial_NA$initial)) {
      spcov_grid[, "dispersion"] <- dispersion_initial_NA$initial
    }

    # split
    cov_grid_splits <- split(spcov_grid, seq_len(NROW(spcov_grid)))
    # apply (add family)
    objvals <- vapply(
      X = cov_grid_splits, FUN = eval_grid_glm, FUN.VALUE = numeric(1),
      data_object = data_object, spcov_type = class(spcov_initial_NA),
      family = class(dispersion_initial_NA), estmethod = estmethod, dist_matrix_list = W
    )

    # find minimum value and record parameters
    min_params <- unlist(cov_grid_splits[[which.min(objvals)]])
    spcov_params <- min_params[c("de", "ie", "range", "extra")]
    spcov_initial_NA$initial <- spcov_params
    dispersion_params <- min_params[c("dispersion")]
    dispersion_initial_NA$initial <- dispersion_params
    # return the best parameters
    best_params <- list(
      spcov_initial_val = spcov_initial_NA,
      dispersion_initial_val = dispersion_initial_NA, randcov_initial_val = NULL
    )
  } else {

    # randcov vars names
    randcov_names <- data_object$randcov_names
    # find number of random effects
    nvar_randcov <- length(randcov_names)

    # spatially dominant grid
    ## spatial components 90% of variance
    spcov_grid[, c("de", "ie", "extra")] <- 0.9 * spcov_grid[, c("de", "ie", "extra")]

    # replace with initial
    for (x in names(spcov_grid)) {
      if (!is.na(spcov_initial_NA$initial[[x]])) {
        spcov_grid[, x] <- spcov_initial_NA$initial[[x]]
      }
    }

    ## random effects 10% of variance and evenly spread
    for (x in randcov_names) {
      if (!is.na(randcov_initial_NA$initial[[x]])) {
        spcov_grid[, x] <- randcov_initial_NA$initial[[x]]
      } else {
        spcov_grid[, x] <- 0.1 * ns2 / nvar_randcov
      }
    }

    # Evenly dominated grid
    ## find number of spatial variance parameters
    nvar_spcov <- 2
    ## find all variance parameters
    nvar_cov <- nvar_spcov + nvar_randcov
    ## keep only spatial cases with variance spread out
    evencov_grid <- spcov_grid_init[spcov_grid_init$de == spcov_grid_init$ie, , drop = FALSE]
    ## scale to incorporate overall variance
    evencov_grid[, c("de", "ie", "extra")] <- nvar_spcov / nvar_cov * evencov_grid[, c("de", "ie", "extra")]

    # replace with initial values
    ## spatial
    for (x in names(evencov_grid)) {
      if (!is.na(spcov_initial_NA$initial[[x]])) {
        evencov_grid[, x] <- spcov_initial_NA$initial[[x]]
      }
    }
    ## random
    for (x in randcov_names) {
      if (!is.na(randcov_initial_NA$initial[[x]])) {
        evencov_grid[, x] <- randcov_initial_NA$initial[[x]]
      } else {
        evencov_grid[, x] <- 1 / nvar_cov * ns2
      }
    }
    # find unique combinations
    evencov_grid <- unique(evencov_grid)

    # random dominated grid
    ## spatial component 10%
    randcov_grid_spcov <- spcov_grid_init[spcov_grid_init$de == spcov_grid_init$ie & spcov_grid_init$range == min(spcov_grid_init$range), , drop = FALSE]
    randcov_grid_spcov[, c("de", "ie", "extra")] <- 0.1 * randcov_grid_spcov[, c("de", "ie", "extra")]

    # replace spatial initial values
    for (x in names(randcov_grid_spcov)) {
      if (!is.na(spcov_initial_NA$initial[[x]])) {
        randcov_grid_spcov[, x] <- spcov_initial_NA$initial[[x]]
      }
    }
    # find unique values
    randcov_grid_spcov <- unique(randcov_grid_spcov)

    # random dominant grid
    randcov_grid_randcov <- as.data.frame(as.list(rep(1 / nvar_randcov, nvar_randcov)))
    names(randcov_grid_randcov) <- randcov_names
    ## if there is more than one random effect, split it up into relevant scenarios
    if (nvar_randcov > 1) {
      ## set 0.1 for all proportions
      extra_grid_randcov <- lapply(seq_len(nvar_randcov), function(x) rep(0.1 * 1 / (nvar_randcov - 1), nvar_randcov))
      extra_grid_randcov <- do.call("rbind", extra_grid_randcov)
      ## fill in 0.9 for one proportion in each row
      diag(extra_grid_randcov) <- 0.9
      extra_grid_randcov <- as.data.frame(extra_grid_randcov)
      names(extra_grid_randcov) <- randcov_names
    } else {
      extra_grid_randcov <- NULL
    }
    ## give the random effects 90% of the variance
    randcov_grid_randcov <- 0.9 * ns2 * rbind(randcov_grid_randcov, extra_grid_randcov)
    for (x in randcov_names) {
      if (!is.na(randcov_initial_NA$initial[[x]])) {
        randcov_grid_randcov[, x] <- randcov_initial_NA$initial[[x]]
      }
    }

    ## bind together and replicate
    randcov_grid_spcov_rep <- do.call("rbind", replicate(NROW(randcov_grid_randcov), randcov_grid_spcov, simplify = FALSE))
    randcov_grid_randcov_rep <- do.call("rbind", replicate(NROW(randcov_grid_spcov), randcov_grid_randcov, simplify = FALSE))
    randcov_grid <- cbind(randcov_grid_spcov_rep, randcov_grid_randcov_rep)

    # bind together all grids
    cov_grid <- rbind(spcov_grid, evencov_grid, randcov_grid)
    cov_grid <- unique(cov_grid)

    # add dispersion
    cov_grid$dispersion <- 1
    if (!is.na(dispersion_initial_NA$initial)) {
      cov_grid[, "dispersion"] <- dispersion_initial_NA$initial
    }

    cov_grid_splits <- split(cov_grid, seq_len(NROW(cov_grid)))

    # apply
    objvals <- vapply(
      X = cov_grid_splits, FUN = eval_grid_glm, FUN.VALUE = numeric(1),
      data_object = data_object, spcov_type = class(spcov_initial_NA),
      family = class(dispersion_initial_NA), estmethod = estmethod, dist_matrix_list = W
    )
    # find minimum value and record parameters
    min_params <- unlist(cov_grid_splits[[which.min(objvals)]])
    spcov_params <- min_params[c("de", "ie", "range", "extra")]
    # return the spatial parameters
    spcov_initial_NA$initial <- spcov_params
    # return the dispersion parameter
    dispersion_params <- min_params[c("dispersion")]
    dispersion_initial_NA$initial <- dispersion_params
    # return the random effect parameters
    randcov_initial_NA$initial <- randcov_params(min_params[randcov_names])
    # return the best parameters
    best_params <- list(
      spcov_initial_val = spcov_initial_NA,
      dispersion_initial_val = dispersion_initial_NA, randcov_initial_val = randcov_initial_NA
    )
  }
  # return the best parameters
  best_params
}

#' @export
cov_initial_search_glm.sar <- cov_initial_search_glm.car

eval_grid <- function(cov_grid_split, data_object, spcov_type,
                      estmethod, dist_matrix_list,
                      weights, esv, dist_vector, residual_vector2) {

  # convert list structure to a vector
  cov_grid <- unlist(cov_grid_split)

  # find spatial covariance parameter vector
  spcov_grid <- cov_grid[c("de", "ie", "range", "extra", "rotate", "scale")]
  spcov_grid <- spcov_grid[!is.na(spcov_grid)]
  spcov_params_val <- do.call("spcov_params", c(list(spcov_type = spcov_type), as.list(spcov_grid)))

  # find REML or ML objective function value
  if (estmethod %in% c("reml", "ml")) {

    # incorporate random effects if necessary
    if (is.null(data_object$randcov_initial)) {
      randcov_params_val <- NULL
    } else {
      randcov_names <- data_object$randcov_names
      randcov_params_val <- randcov_params(cov_grid[randcov_names], nm = randcov_names)
    }

    # incorporate anisotropy if necessary
    if (data_object$anisotropy) {
      new_coords_list_q1 <- lapply(data_object$obdata_list, transform_anis, data_object$xcoord, data_object$ycoord,
        rotate = spcov_params_val[["rotate"]], scale = spcov_params_val[["scale"]]
      )
      dist_matrix_list_q1 <- lapply(new_coords_list_q1, function(x) spdist(xcoord_val = x$xcoord_val, ycoord_val = x$ycoord_val))

      # compute relevant products
      gll_prods_q1 <- gloglik_products(
        spcov_params_val, data_object, estmethod,
        dist_matrix_list_q1, randcov_params_val
      )

      # find -2loglik
      objval_q1 <- get_minustwologlik(gll_prods_q1, estmethod, data_object$n, data_object$p, spcov_profiled = FALSE, randcov_profiled = FALSE)


      new_coords_list_q2 <- lapply(data_object$obdata_list, transform_anis, data_object$xcoord, data_object$ycoord,
        rotate = abs(pi - spcov_params_val[["rotate"]]), scale = spcov_params_val[["scale"]]
      )
      dist_matrix_list_q2 <- lapply(new_coords_list_q2, function(x) spdist(xcoord_val = x$xcoord_val, ycoord_val = x$ycoord_val))

      # compute relevant products
      gll_prods_q2 <- gloglik_products(
        spcov_params_val, data_object, estmethod,
        dist_matrix_list_q2, randcov_params_val
      )

      # find -2loglik
      objval_q2 <- get_minustwologlik(gll_prods_q2, estmethod, data_object$n,
        data_object$p,
        spcov_profiled = FALSE, randcov_profiled = FALSE
      )

      objval <- min(objval_q1, objval_q2)
    } else {

      # compute relevant products
      gll_prods <- gloglik_products(
        spcov_params_val, data_object, estmethod,
        dist_matrix_list, randcov_params_val
      )

      # find -2loglik
      objval <- get_minustwologlik(gll_prods, estmethod, data_object$n,
        data_object$p,
        spcov_profiled = FALSE, randcov_profiled = FALSE
      )
    }
  } else if (estmethod == "sv-wls") {
    # find sv-wls objective function value
    objval <- get_svloss(spcov_params_val,
      esv = esv,
      weights = weights
    )
  } else if (estmethod == "sv-cl") {
    # find cl objective function value
    objval <- get_glogclikloss(spcov_params_val, residual_vector2, dist_vector)
  }
  objval
}

eval_grid_glm <- function(cov_grid_split, data_object, spcov_type,
                          family, estmethod, dist_matrix_list) {

  # convert list structure to a vector
  cov_grid <- unlist(cov_grid_split)

  # find spatial covariance parameter vector
  spcov_grid <- cov_grid[c("de", "ie", "range", "extra", "rotate", "scale")]
  spcov_grid <- spcov_grid[!is.na(spcov_grid)]
  spcov_params_val <- do.call("spcov_params", c(list(spcov_type = spcov_type), as.list(spcov_grid)))

  # dispersion parameter vector
  dispersion_params_val <- dispersion_params(family = data_object$family, dispersion = cov_grid[["dispersion"]])

  # find REML or ML objective function value

  # incorporate random effects if necessary
  if (is.null(data_object$randcov_initial)) {
    randcov_params_val <- NULL
  } else {
    randcov_names <- data_object$randcov_names
    randcov_params_val <- randcov_params(cov_grid[randcov_names], nm = randcov_names)
  }

  # incorporate anisotropy if necessary
  if (data_object$anisotropy) {
    new_coords_list_q1 <- lapply(data_object$obdata_list, transform_anis, data_object$xcoord, data_object$ycoord,
      rotate = spcov_params_val[["rotate"]], scale = spcov_params_val[["scale"]]
    )
    dist_matrix_list_q1 <- lapply(new_coords_list_q1, function(x) spdist(xcoord_val = x$xcoord_val, ycoord_val = x$ycoord_val))

    # compute relevant products
    lapll_prods_q1 <- laploglik_products(
      spcov_params_val, dispersion_params_val, data_object, estmethod,
      dist_matrix_list_q1, randcov_params_val
    )

    # find -2loglik
    objval_q1 <- get_minustwolaploglik(lapll_prods_q1, estmethod, data_object$n, data_object$p, spcov_profiled = FALSE, randcov_profiled = FALSE)


    new_coords_list_q2 <- lapply(data_object$obdata_list, transform_anis, data_object$xcoord, data_object$ycoord,
      rotate = abs(pi - spcov_params_val[["rotate"]]), scale = spcov_params_val[["scale"]]
    )
    dist_matrix_list_q2 <- lapply(new_coords_list_q2, function(x) spdist(xcoord_val = x$xcoord_val, ycoord_val = x$ycoord_val))

    # compute relevant products
    lapll_prods_q2 <- laploglik_products(
      spcov_params_val, dispersion_params_val, data_object, estmethod,
      dist_matrix_list_q2, randcov_params_val
    )

    # find -2loglik
    objval_q2 <- get_minustwolaploglik(lapll_prods_q2, estmethod, data_object$n,
      data_object$p,
      spcov_profiled = FALSE, randcov_profiled = FALSE
    )


    objval <- min(objval_q1, objval_q2)
  } else {


    # compute relevant products
    lapll_prods <- laploglik_products(
      spcov_params_val, dispersion_params_val, data_object, estmethod,
      dist_matrix_list, randcov_params_val
    )

    # find -2loglik
    objval <- get_minustwolaploglik(lapll_prods, estmethod, data_object$n,
      data_object$p,
      spcov_profiled = FALSE, randcov_profiled = FALSE
    )
  }
  objval
}
