#' Fit random forest residual spatial autoregressive models
#'
#' @param formula A two-sided linear formula describing the fixed effect structure
#'   of the model, with the response to the left of the \code{~} operator and
#'   the terms on the right, separated by \code{+} operators.
#' @param data A data frame or \code{sf} object object that contains
#'   the variables in \code{fixed}, \code{random}, and \code{partition_factor}
#'   as well as geographical information. If an \code{sf} object is
#'   provided with \code{POINT} geometries, the x-coordinates and y-coordinates
#'   are used directly. If an \code{sf} object is
#'   provided with \code{POLYGON} geometries, the x-coordinates and y-coordinates
#'   are taken as the centroids of each polygon.
#' @param ... Additional named arguments to \code{ranger::ranger} or [spautor()].
#'
#' @details 1. Find fitted values from a random forest model. 2. Fit a spatial
#'   autoregressive model to the residuals of the random forest model.
#'
#' @return An \code{spmodRF} object to be used with \code{predict()}. There are
#'   three elements: \code{ranger}, the output from fitting the mean model with
#'   \code{ranger::ranger}; \code{spmod}, the output from fitting the spatial
#'   linear model to the ranger residuals; and \code{newdata}, the \code{newdata}
#'   object, if relevant.
#'
#' @export
#'
#' @examples
#' \donttest{
#' seal$x <- rnorm(NROW(seal)) # add dummy variable
#' sprfmod <- spautorRF(log_trend ~ x, data = seal, spcov_type = "car")
#' predict(sprfmod)
#' }
spautorRF <- function(formula, data, ...) {

  # check to see if ranger installed
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Install the ranger package before using spautorRF", call. = FALSE)
  } else {

    # save calls for later (NSE can be a bit frustrating)
    ranger_call <- call("ranger", formula = substitute(formula), data = substitute(data), quote(...))
    spautor_call <- call("spautor", formula = .ranger_resid ~ 1, data = substitute(data), quote(...))

    # find NA values for newdata if required
    na_index <- is.na(model.response(model.frame(formula, data, na.action = na.pass)))

    # make sure at least one missing value to predict
    if (any(na_index)) {
      newdata <- data[na_index, , drop = FALSE]
      data <- data[!na_index, , drop = FALSE]
    } else {
      stop("No missing data to predict", call. = FALSE)
    }

    # get ... objects
    call_list <- as.list(match.call())[-1]
    call_list <- call_list[!names(call_list) %in% c("formula", "data")]
    penv <- parent.frame()

    # save ranger ... objects
    ranger_names <- names(formals(ranger::ranger))
    ranger_args <- call_list[names(call_list) %in% ranger_names]

    # perform random forest
    ranger_out <- do.call(ranger::ranger, c(list(formula = formula, data = as.data.frame(data)), ranger_args), envir = penv)
    ranger_out$call <- ranger_call

    # get ... objects
    spautor_names <- names(formals(spmodel::spautor))
    spautor_args <-  call_list[names(call_list) %in% spautor_names]
    # find residuals
    data$.ranger_resid <- model.response(model.frame(formula, data)) - ranger_out$predictions
    newdata$.ranger_resid <- NA
    # perform spautor
    # reset newdata
    data <- rbind(data, newdata)
    # putting back in order
    data <- data[order(c(which(!na_index), which(na_index))), , drop = FALSE]
    # perform spautor
    spmod_out <- do.call(spmodel::spautor, c(list(formula = .ranger_resid ~ 1, data = data), spautor_args), envir = penv)
    spmod_out$call <- spautor_call

  }
  # output list with names and class
  structure(list(ranger = ranger_out, spmod = spmod_out, newdata = spmod_out$newdata), class = "spmodRF")
}
