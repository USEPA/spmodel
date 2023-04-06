#' @param type The scale (\code{response} or \code{link}) of predictions obtained
#'   using \code{spglm()} or \code{spgautor} objects.
#' @param newdata_size The \code{size} value for each observation in \code{newdata}
#'   used when predicting for the binomial family.
#' @param var_correct A logical indicating whether to return the corrected prediction
#'   variances when predicting via models fit using \code{spglm()} or \code{spgautor()}. The default is
#'   \code{TRUE}.
#' @rdname augment.spmodel
#' @method augment spglm
#' @order 3
#' @export
augment.spglm <- function(x, drop = TRUE, newdata = NULL, type = c("link", "response"), se_fit = FALSE,
                         interval = c("none", "confidence", "prediction"), newdata_size,
                         level = 0.95, local = local, var_correct = TRUE, ...) {

  type <- match.arg(type)
  interval <- match.arg(interval)

  # set data and newdata
  if (is.null(newdata)) {
    if (drop) {
      data <- cbind(model.frame(x), x$obdata[, c(x$xcoord, x$ycoord)])
      # keep_cols <- colnames(model.frame(x))
      # data <- x$obdata[, c(keep_cols, x$xcoord, x$ycoord)]
    } else {
      data <- x$obdata
    }
  } else {
    data <- model.frame(x)
  }

  if (is.null(newdata)) {
    augment_data <- tibble::tibble(.fitted = fitted(x))
    if (se_fit) {
      preds_data <- predict(x, newdata = data, se.fit = se_fit, interval = "confidence", ...)
      augment_data$.se.fit <- preds_data$se.fit
    }
    tibble_out <- tibble::tibble(cbind(data, augment_data, influence(x)))
  } else {
    if (missing(newdata_size)) newdata_size <- NULL
    if (missing(local)) local <- NULL
    preds_newdata <- predict(x, newdata = newdata, type = type, se.fit = se_fit, interval = interval,
                             newdata_size = newdata_size, level = level, local = local, var_correct = var_correct)
    if (se_fit) {
      if (interval %in% c("confidence", "prediction")) {
        augment_newdata <- tibble::tibble(
          .fitted = preds_newdata$fit[, "fit"]
        )
        augment_newdata$.lower <- preds_newdata$fit[, "lwr"]
        augment_newdata$.upper <- preds_newdata$fit[, "upr"]
      } else {
        augment_newdata <- tibble::tibble(
          .fitted = preds_newdata$fit
        )
      }
      augment_newdata$.se.fit <- preds_newdata$se.fit
    } else {
      if (interval %in% c("confidence", "prediction")) {
        augment_newdata <- tibble::tibble(
          .fitted = preds_newdata[, "fit"]
        )
        augment_newdata$.lower <- preds_newdata[, "lwr"]
        augment_newdata$.upper <- preds_newdata[, "upr"]
      } else {
        augment_newdata <- tibble::tibble(
          .fitted = preds_newdata
        )
      }
    }

    # inheritance for sf or sp objects
    attr_sp <- attr(class(newdata), "package")
    if (!is.null(attr_sp) && length(attr_sp) == 1 && attr_sp == "sp") {
      stop("sf objects must be used instead of sp objects. To convert your sp object into an sf object, run sf::st_as_sf().", call. = FALSE)
    }

    if (inherits(newdata, "sf")) {

      newdata <- suppressWarnings(sf::st_centroid(newdata))

      newdata <- sf_to_df(newdata)
      names(newdata)[[which(names(newdata) == ".xcoord")]] <-
        as.character(x$xcoord) # only relevant if newdata is sf data is not
      names(newdata)[[which(names(newdata) == ".ycoord")]] <-
        as.character(x$ycoord) # only relevant if newdata is sf data is not
    }

    tibble_out <- tibble::tibble(cbind(newdata, augment_newdata))
  }


  # if (x$is_sf && requireNamespace("sf", quietly = TRUE)) {
  if (x$is_sf) {
    # sf installed
    if (inherits(newdata, "sf")) {
      tibble_out <- sf::st_as_sf(tibble_out,
                                 sf_column_name = x$sf_column_name,
                                 crs = x$crs
      )
    } else {
      tibble_out <- sf::st_as_sf(tibble_out,
                                 coords = c(x$xcoord, x$ycoord),
                                 crs = x$crs
      ) # may want to double check this for mismatching geometry names
      # i.e. geometry names that are not geometry
    }
    names(tibble_out)[names(tibble_out) == "geometry"] <- x$sf_column_name
    sf::st_geometry(tibble_out) <- x$sf_column_name
  }
  tibble_out
}

#' @rdname augment.spmodel
#' @method augment spgautor
#' @order 4
#' @export
augment.spgautor <- function(x, drop = TRUE, newdata = NULL, type = c("link", "response"), se_fit = FALSE,
                            interval = c("none", "confidence", "prediction"), newdata_size,
                            level = 0.95, local, var_correct = TRUE, ...) {

  type <- match.arg(type)
  interval <- match.arg(interval)

  # set data and newdata
  if (is.null(newdata)) {
    if (drop) {
      if (x$is_sf) {
        data_sf <- x$data[x$observed_index, x$sf_column_name, drop = FALSE]
      }
      data <- model.frame(x)
      # keep_cols <- colnames(model.frame(x))
      # data <- data[, keep_cols, drop = FALSE]
    } else {
      data <- x$data[x$observed_index, , drop = FALSE]
    }
  }


  if (is.null(newdata)) {
    augment_data <- tibble::tibble(.fitted = fitted(x))
    if (se_fit) {
      preds_data <- predict(x, newdata = data, se.fit = se_fit, interval = "confidence", ...)
      augment_data$.se.fit <- preds_data$se.fit
    }
    if (x$is_sf && drop) {
      tibble_out <- tibble::tibble(cbind(data, augment_data, influence(x), data_sf))
    } else {
      tibble_out <- tibble::tibble(cbind(data, augment_data, influence(x)))
    }
  } else {
    if (missing(newdata_size)) newdata_size <- NULL
    if (missing(local)) local <- NULL
    preds_newdata <- predict(x, newdata = newdata, type = type, se.fit = se_fit, interval = interval,
                             newdata_size = newdata_size, level = level, local = local, var_correct = var_correct)
    if (se_fit) {
      if (interval %in% c("confidence", "prediction")) {
        augment_newdata <- tibble::tibble(
          .fitted = preds_newdata$fit[, "fit"]
        )
        augment_newdata$.lower <- preds_newdata$fit[, "lwr"]
        augment_newdata$.upper <- preds_newdata$fit[, "upr"]
      } else {
        augment_newdata <- tibble::tibble(
          .fitted = preds_newdata$fit
        )
      }
      augment_newdata$.se.fit <- preds_newdata$se.fit
    } else {
      if (interval %in% c("confidence", "prediction")) {
        augment_newdata <- tibble::tibble(
          .fitted = preds_newdata[, "fit"]
        )
        augment_newdata$.lower <- preds_newdata[, "lwr"]
        augment_newdata$.upper <- preds_newdata[, "upr"]
      } else {
        augment_newdata <- tibble::tibble(
          .fitted = preds_newdata
        )
      }
    }
    tibble_out <- tibble::tibble(cbind(newdata, augment_newdata))
  }

  if (x$is_sf) {
    # sf installed
    tibble_out <- sf::st_as_sf(tibble_out,
                               sf_column_name = x$sf_column_name,
                               crs = x$crs
    )
    names(tibble_out)[names(tibble_out) == "geometry"] <- x$sf_column_name
    sf::st_geometry(tibble_out) <- x$sf_column_name
  }
  tibble_out
}

