#' Augment data with information from fitted model objects
#'
#' @description Augment accepts a fitted model object and a data set and adds
#'   information about each observation in the data set. New columns always
#'   begin with a \code{.} prefix to avoid overwriting columns in the original
#'   data set.
#'
#'   Augment behaves differently depending on whether the original data or new data
#'   requires augmenting. Typically, when augmenting the original data, only the fitted
#'   model object is specified, and when augmenting new data, the fitted model object
#'   and \code{newdata} is specified. When augmenting the original data, diagnostic
#'   statistics are augmented to each row in the data set. When augmenting new data,
#'   predictions and optional intervals or standard errors are augmented to each
#'   row in the new data set.
#'
#' @param x A fitted model object from [splm()] or [spautor()].
#' @param drop A logical indicating whether to drop extra variables in the
#'   fitted model object \code{x} when augmenting. The default for \code{drop} is \code{TRUE}.
#'   \code{drop} is ignored if augmenting \code{newdata}.
#' @param newdata A data frame or tibble containing observations requiring prediction.
#'   All of the original explanatory variables used to create the fitted model object \code{x}
#'   must be present in \code{newdata}. Defaults to \code{NULL}, which indicates
#'   that nothing has been passed to \code{newdata}.
#' @param se_fit Logical indicating whether or not a \code{.se.fit} column should
#'   be added to augmented output. Passed to \code{predict()} and
#'   defaults to \code{FALSE}.
#' @param interval Character indicating the type of confidence interval columns to
#'   add to the augmented \code{newdata} output. Passed to \code{predict()} and defaults
#'   to \code{"none"}.
#' @param level Tolerance/confidence level. The default is \code{0.95}.
#' @param local A list or logical. If a list, specific list elements described
#'   in [predict.spmodel()] control the big data approximation behavior.
#'   If a logical, \code{TRUE} chooses default list elements for the list version
#'   of \code{local} as specified in [predict.spmodel()]. Defaults to \code{FALSE},
#'   which performs exact computations.
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @details \code{augment()} returns a tibble with the same class as
#'   \code{data}. That is, if \code{data} is
#'   an \code{sf} object, then the augmented object (obtained via \code{augment(x)})
#'   will be an \code{sf} object as well. When augmenting \code{newdata}, the
#'   augmented object has the same class as \code{data}.
#'
#'   Missing response values from the original data can be augmented as if
#'   they were a \code{newdata} object by providing \code{x$newdata} to the
#'   \code{newdata} argument (where \code{x} is the name of the fitted model
#'   object). This is the only way to compute predictions for
#'   [spautor()] and [spgautor()] fitted model objects.
#'
#' @return When augmenting the original data set, a tibble with additional columns
#'   \itemize{
#'     \item{\code{.fitted}}{ Fitted value}
#'     \item{\code{.resid}}{ Response residual (the difference between observed and fitted values)}
#'     \item{\code{.hat}}{ Leverage (diagonal of the hat matrix)}
#'     \item{\code{.cooksd}}{ Cook's distance}
#'     \item{\code{.std.resid}}{ Standardized residuals}
#'     \item{\code{.se.fit}}{ Standard error of the fitted value.}
#'   }
#'
#'   When augmenting a new data set, a tibble with additional columns
#'   \itemize{
#'     \item{\code{.fitted}}{ Predicted (or fitted) value}
#'     \item{\code{.lower}}{ Lower bound on interval}
#'     \item{\code{.upper}}{ Upper bound on interval}
#'     \item{\code{.se.fit}}{ Standard error of the predicted (or fitted) value}
#'   }
#'
#' @name augment.spmodel
#' @method augment splm
#' @order 1
#' @export
#'
#' @seealso [tidy.spmodel()] [glance.spmodel()] [predict.spmodel()]
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' augment(spmod)
#' spmod_sulf <- splm(sulfate ~ 1, data = sulfate, spcov_type = "exponential")
#' augment(spmod_sulf)
#' augment(spmod_sulf, newdata = sulfate_preds)
#' # missingness in original data
#' spmod_seal <- spautor(log_trend ~ 1, data = seal, spcov_type = "car")
#' augment(spmod_seal)
#' augment(spmod_seal, newdata = spmod_seal$newdata)
augment.splm <- function(x, drop = TRUE, newdata = NULL, se_fit = FALSE,
                         interval = c("none", "confidence", "prediction"), level = 0.95,
                         local, ...) {

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
    if (missing(local)) local <- NULL
    preds_newdata <- predict(x, newdata = newdata, se.fit = se_fit, interval = interval,
                             level = level, local = local)
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
#' @method augment spautor
#' @order 2
#' @export
augment.spautor <- function(x, drop = TRUE, newdata = NULL, se_fit = FALSE,
                            interval = c("none", "confidence", "prediction"),
                            level = 0.95, local, ...) {

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
    if (missing(local)) local <- NULL
    preds_newdata <- predict(x, newdata = newdata, se.fit = se_fit, interval = interval,
                             level = level, local = local)
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
