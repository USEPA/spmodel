#' @export
predict.spglm <- function(object, newdata, type = c("link", "response"), se.fit = FALSE, interval = c("none", "confidence", "prediction"),
                         newdata_size, level = 0.95, local, ...) {


  # match type argument so the two display
  type <- match.arg(type)

  # match interval argument so the three display
  interval <- match.arg(interval)

  # deal with newdata_size
  if (missing(newdata_size)) newdata_size <- NULL

  # deal with local
  if (missing(local)) {
    local <- NULL
  }

  # error if newdata missing from arguments and object
  if (missing(newdata) && is.null(object$newdata)) {
    stop("No missing data to predict. newdata must be specified in the newdata argument or object$newdata must be non-NULL.", call. = FALSE)
  }

  # rename relevant quantities
  obdata <- object$obdata
  xcoord <- object$xcoord
  ycoord <- object$ycoord

  # write newdata if predicting missing data
  if (missing(newdata)) {
    add_newdata_rows <- TRUE
    newdata <- object$newdata
  } else {
    add_newdata_rows <- FALSE
  }

  # set newdata_size if needed
  if (is.null(newdata_size) && object$family == "binomial") {
    newdata_size <- rep(1, NROW(newdata))
  }

  # deal with local
  if (is.null(local)) {
    if (object$n > 5000 || NROW(newdata) > 5000) {
      local <- TRUE
      message("Because either the sample size of the fitted model object or the number of desired predictions exceeds 5000, we are setting local = TRUE to perform computationally efficient approximations. To override this behavior and compute the exact solution, rerun predict() with local = FALSE. Be aware that setting local = FALSE may result in exceedingly long computational times.")
    } else {
      local <- FALSE
    }
  }

  # save spcov param vector
  spcov_params_val <- coef(object, type = "spcov")

  # save dispersion param vector FIX THIS TO USE COEF LATER
  dispersion_params_val <- object$coefficients$dispersion

  # save randcov param vector
  randcov_params_val <- coef(object, type = "randcov")

}
