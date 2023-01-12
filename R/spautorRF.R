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
#'   (named \code{spautor} or \code{spautor_list}), and an object can contain data for
#'   locations at which to predict (called \code{newdata}). The \code{newdata}
#'   object contains the set of
#'   observations in \code{data} whose response variable is \code{NA}.
#'   If \code{spcov_type} or \code{spcov_initial} (which are passed to [spautor()])
#'   are length one, the list has class \code{spautorRF} and the spatial linear
#'   model object fit to the residuals is called \code{spautor}, which has
#'   class \code{spautor}. If
#'   \code{spcov_type} or \code{spcov_initial} are length greater than one, the
#'   list has class \code{spautorRF_list} and the spatial linear model object
#'   fit to the residuals is called \code{spautor_list}, which has class \code{spautor_list}.
#'   and contains several objects, each with class \code{spautor}.
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
#' seal$var <- rnorm(NROW(seal)) # add noise variable
#' sprfmod <- spautorRF(log_trend ~ var, data = seal, spcov_type = "car")
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
    if (inherits(data, "sf")) {
      model_resp <- model.response(model.frame(formula, sf::st_drop_geometry(data), na.action = na.pass))
    } else {
      model_resp <- model.response(model.frame(formula, data, na.action = na.pass))
    }
    na_index <- is.na(model_resp)
    resp <- model_resp[!na_index]

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
    ## ranger needs a data frame and model frame needs non list objects in formula
    if (inherits(data, "sf")) {
      ranger_out <- do.call(ranger::ranger, c(list(formula = formula, data = sf::st_drop_geometry(data)), ranger_args), envir = penv)
    } else {
      ranger_out <- do.call(ranger::ranger, c(list(formula = formula, data = data), ranger_args), envir = penv)
    }
    ranger_out$call <- NA

    # get ... objects
    spautor_names <- names(formals(spmodel::spautor))
    spautor_args <-  call_list[names(call_list) %in% spautor_names]
    # find residuals
    data$.ranger_resid <- resp - ranger_out$predictions
    newdata$.ranger_resid <- NA
    # perform spautor
    # reset newdata
    data <- rbind(data, newdata)
    # putting back in order
    data <- data[order(c(which(!na_index), which(na_index))), , drop = FALSE]
    # perform splm
    spautor_out <- do.call(spmodel::spautor, c(list(formula = .ranger_resid ~ 1, data = data), spautor_args), envir = penv)
    if (inherits(spautor_out, "spautor")) {
      spautor_out$call <- NA
      # output list with names and class
      sprf_out <- structure(list(call = match.call(), ranger = ranger_out, spautor = spautor_out, newdata = newdata), class = "spautorRF")
    } else {
      spautor_out <- lapply(spautor_out, function(x) {
        x$call <- NA
        x
      })
      sprf_out <- structure(list(call = match.call(), ranger = ranger_out, spautor_list = spautor_out, newdata = newdata), class = "spautorRF_list")
    }
  }
  # return object
  sprf_out
}
