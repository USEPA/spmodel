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
    dotlist <- as.list(substitute(alist(...)))[-1]
    dotlist_names <- names(dotlist)

    # save ranger ... objects
    ranger_names <- names(formals(ranger::ranger))
    ranger_args <- dotlist[dotlist_names %in% ranger_names]

    # perform random forest
    ranger_out <- do.call(ranger::ranger, c(list(formula = formula, data = as.data.frame(data)), ranger_args))
    ranger_out$call <- ranger_call

    # get ... objects
    spautor_args <-  dotlist[!dotlist_names %in% ranger_names]
    # find residuals
    data$.ranger_resid <- model.response(model.frame(formula, data)) - ranger_out$predictions
    newdata$.ranger_resid <- NA
    # perform spautor
    # reset newdata
    data <- rbind(data, newdata)
    # putting back in order
    data <- data[order(c(which(!na_index), which(na_index))), , drop = FALSE]
    # perform spautor
    spmod_out <- do.call(spmodel::spautor, c(list(formula = .ranger_resid ~ 1, data = data), spautor_args))
    spmod_out$call <- spautor_call

  }
  # output list with names and class
  structure(list(ranger = ranger_out, spmod = spmod_out, newdata = spmod_out$newdata), class = "spmodRF")
}
