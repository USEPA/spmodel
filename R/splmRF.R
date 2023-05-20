#' Fit random forest spatial residual models
#'
#' @description Fit random forest spatial residual models
#'   for point-referenced data (i.e., geostatistical models) using
#'   random forest to fit the mean and a spatial linear model to fit the residuals.
#'   The spatial linear model fit to the residuals can incorporate variety of estimation methods,
#'   allowing for random effects, anisotropy, partition factors, and big data methods.
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
#' @param ... Additional named arguments to [ranger::ranger()] or [splm()].
#'
#' @details The random forest residual spatial linear model is described by
#'   Fox et al. (2020). A random forest model is fit to the mean portion of the
#'   model specified by \code{formula} using \code{ranger::ranger()}. Residuals
#'   are computed and used as the response variable in an intercept-only spatial
#'   linear model fit using [splm()]. This model object is intended for use with
#'   \code{predict()} to perform prediction, also called random forest
#'   regression Kriging.
#'
#' @return A list with several elements to be used with \code{predict()}. These
#'   elements include the function call (named \code{call}), the random forest object
#'   fit to the mean (named \code{ranger}),
#'   the spatial linear model object fit to the residuals
#'   (named \code{splm} or \code{splm_list}), and an object can contain data for
#'   locations at which to predict (called \code{newdata}). The \code{newdata}
#'   object contains the set of
#'   observations in \code{data} whose response variable is \code{NA}.
#'   If \code{spcov_type} or \code{spcov_initial} (which are passed to [splm()])
#'   are length one, the list has class \code{splmRF} and the spatial linear
#'   model object fit to the residuals is called \code{splm}, which has
#'   class \code{splm}. If
#'   \code{spcov_type} or \code{spcov_initial} are length greater than one, the
#'   list has class \code{splmRF_list} and the spatial linear model object
#'   fit to the residuals is called \code{splm_list}, which has class \code{splm_list}.
#'   and contains several objects, each with class \code{splm}.
#'
#'
#' An \code{splmRF} object to be used with \code{predict()}. There are
#'   three elements: \code{ranger}, the output from fitting the mean model with
#'   [ranger::ranger()]; \code{splm}, the output from fitting the spatial
#'   linear model to the ranger residuals; and \code{newdata}, the \code{newdata}
#'   object, if relevant.
#'
#' @note This function does not perform any internal scaling. If optimization is not
#'   stable due to large extremely large variances, scale relevant variables
#'   so they have variance 1 before optimization.
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
#' sulfate$var <- rnorm(NROW(sulfate)) # add noise variable
#' sulfate_preds$var <- rnorm(NROW(sulfate_preds)) # add noise variable
#' sprfmod <- splmRF(sulfate ~ var, data = sulfate, spcov_type = "exponential")
#' predict(sprfmod, sulfate_preds)
#' }
splmRF <- function(formula, data, ...) {

  # check to see if ranger installed
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Install the ranger package before using splmRF", call. = FALSE)
  } else {

    # save calls for later (NSE can be a bit frustrating)
    # ranger_call <- call("ranger", formula = substitute(formula), data = substitute(data), quote(...))
    # splm_call <- call("splm", formula = .ranger_resid ~ 1, data = substitute(data), quote(...))

    # find NA values for newdata if required
    if (inherits(data, "sf")) {
      model_resp <- model.response(model.frame(formula, sf::st_drop_geometry(data), na.action = na.pass))
    } else {
      model_resp <- model.response(model.frame(formula, data, na.action = na.pass))
    }
    na_index <- is.na(model_resp)
    resp <- model_resp[!na_index]

    if (any(na_index)) {
      newdata <- data[na_index, , drop = FALSE]
      data <- data[!na_index, , drop = FALSE]
    } else {
      newdata <- NULL
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
    splm_names <- names(formals(spmodel::splm))
    splm_args <-  call_list[names(call_list) %in% splm_names]
    # find residuals
    data$.ranger_resid <- resp - ranger_out$predictions
    # perform splm
    splm_out <- do.call(spmodel::splm, c(list(formula = .ranger_resid ~ 1, data = data), splm_args), envir = penv)
    if (inherits(splm_out, "splm")) {
      splm_out$call <- NA
      # output list with names and class
      sprf_out <- structure(list(call = match.call(), ranger = ranger_out, splm = splm_out, newdata = newdata), class = "splmRF")
    } else {
      splm_out <- lapply(splm_out, function(x) {
        x$call <- NA
        x
      })
      sprf_out <- structure(list(call = match.call(), ranger = ranger_out, splm_list = splm_out, newdata = newdata), class = "splmRF_list")
    }
  }
  # return object
  sprf_out

}
