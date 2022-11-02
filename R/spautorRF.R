#' Fit random forest spatial residual models
#'
#' @description Fit random forest residual spatial linear models
#'   for areal data (i.e., spatial autoregressive models) using
#'   random forest to fit the mean and a spatial linear model to fit the residuals.
#'   The spatial linear model fit to the residuals can incorporate
#'   a variety of estimation methods, allowing for random effects,
#'   partition factors, and row standardization.
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
#' @details The random forest residual spatial linear model is described by
#'   Fox Et al. (2020). A random forest model is fit to the mean portion of the
#'   model specified by \code{formula} using \code{ranger::ranger()}. Residuals
#'   are computed and used as the response variable in an intercept-only spatial
#'   linear model fit using [spautor()]. This model object is intended for use with
#'   \code{predict()} for perform prediction, also called random forest
#'   regression Kriging.
#'
#' @return A list with several elements to be used with \code{predict()}. These
#'   elements include the function call (named \code{call}), the random forest object
#'   fit to the mean (named \code{ranger}),
#'   the spatial linear model object fit to the residuals
#'   (named \code{spmod} or \code{spmod_list}), and an object can contain data for
#'   locations at which to predict (called \code{newdata}). The \code{newdata}
#'   object contains the set of
#'   observations in \code{data} whose response variable is \code{NA}.
#'   If \code{spcov_type} or \code{spcov_initial} (which are passed to [spautor()])
#'   are length one, the list has class \code{spmodRF} and the spatial linear
#'   model object fit to the residuals is called \code{spmod}, which has
#'   class \code{spmod}. If
#'   \code{spcov_type} or \code{spcov_initial} are length greater than one, the
#'   list has class \code{spmodRF_list} and the spatial linear model object
#'   fit to the residuals is called \code{spmod_list}, which has class \code{spmod_list}.
#'   and contains several objects, each with class \code{spmod}.
#'
#' @export
#'
#' @references
#' Fox, E.W., Ver Hoef, J. M., & Olsen, A. R. (2020). Comparing spatial
#'   regression to random forests for large environmental data sets.
#'   \emph{PloS one}, 15(3), e0229509.
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
    # ranger_call <- call("ranger", formula = substitute(formula), data = substitute(data), quote(...))
    # spautor_call <- call("spautor", formula = .ranger_resid ~ 1, data = substitute(data), quote(...))

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
    ranger_out$call <- NA

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
    # perform splm
    spmod_out <- do.call(spmodel::spautor, c(list(formula = .ranger_resid ~ 1, data = data), spautor_args), envir = penv)
    if (inherits(spmod_out, "spmod")) {
      spmod_out$call <- NA
      # output list with names and class
      sprf_out <- structure(list(call = match.call(), ranger = ranger_out, spmod = spmod_out, newdata = newdata), class = "spmodRF")
    } else {
      spmod_out <- lapply(spmod_out, function(x) {
        x$call <- NA
        x
      })
      sprf_out <- structure(list(call = match.call(), ranger = ranger_out, spmod_list = spmod_out, newdata = newdata), class = "spmodRF_list")
    }
  }
  # return object
  sprf_out
}
