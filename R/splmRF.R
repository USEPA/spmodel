splmRF <- function(formula, data, ...) {

  # check to see if ranger installed
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Install the ranger package before using splmRF", call. = FALSE)
  } else {

    # save calls for later (NSE can be a bit frustrating)
    ranger_call <- call("ranger", formula = substitute(formula), data = substitute(data), quote(...))
    splm_call <- call("splm", formula = .ranger_resid ~ 1, data = substitute(data), quote(...))

    # find NA values for newdata if required
    na_index <- is.na(model.response(model.frame(formula, data, na.action = na.pass)))
    if (any(na_index)) {
      newdata <- data[na_index, , drop = FALSE]
      data <- data[!na_index, , drop = FALSE]
    }

    # get ... objects
    dotlist <- as.list(substitute(alist(...)))[-1]
    dotlist_names <- names(dotlist)

    # save ranger ... objects
    ranger_names <- names(formals(ranger::ranger))
    ranger_args <- dotlist[dotlist_names %in% ranger_names]

    # perform random forest
    ranger_out <- do.call(ranger::ranger, c(list(formula = formula, data = as.data.frame(data)), ranger_args))
    ranger_out$call <- ranger_call

    # get ... objects
    splm_args <-  dotlist[!dotlist_names %in% ranger_names]
    # find residuals
    data$.ranger_resid <- model.response(model.frame(formula, data)) - ranger_out$predictions
    # perform splm
    spmod_out <- do.call(spmodel::splm, c(list(formula = .ranger_resid ~ 1, data = data), splm_args))
    spmod_out$call <- splm_call

  }
  # output list with names and class
  structure(list(ranger = ranger_out, spmod = spmod_out, newdata = newdata), class = "spmodRF")
}

#' @export
predict.spmodRF <- function(object, newdata, ...) {

  # check to see if ranger installed
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Install the ranger package before using predict with an splmRF object", call. = FALSE)
  } else {

    # find newdata if required (conditional if spautor? this would force them to use newdata object)
    if ((missing(newdata) && !is.null(object$newdata)) || object$spmod$fn == "spautor") {
      newdata <- object$newdata
    }

    # get ... objects
    dotlist <- as.list(substitute(alist(...)))[-1]
    dotlist_names <- names(dotlist)

    # hardcode ranger names because of predict.ranger export issue
    ranger_names <- c("predict.all", "num.trees", "se.method", "quantiles", "what", "seed", "num.threads", "verbose")
    ranger_args <- dotlist[dotlist_names %in% ranger_names]

    # do random forest prediction
    ranger_pred <- do.call(predict, c(list(object = object$ranger, data = as.data.frame(newdata), type = "response"), ranger_args))

    # find spmod names
    spmod_args <-  dotlist[!dotlist_names %in% ranger_names]
    # do splm prediction
    spmod_pred <- do.call(predict, c(list(object = object$spmod, newdata = newdata), spmod_args))

  }
  # obtain final predictions
  ranger_pred$predictions + spmod_pred
}
