#' Plot fitted model diagnostics
#'
#' @description Plot fitted model diagnostics such as residuals vs fitted values,
#'   quantile-quantile, scale-location, Cook's distance, residuals vs leverage,
#'   Cook's distance vs leverage, a fitted spatial covariance function, and a
#'   fitted anisotropic level curve of equal correlation.
#'
#' @param x A fitted model object from [splm()] or [spautor()].
#' @param which An integer vector taking on values between 1 and 7, which indicates
#'   the plots to return. Available plots are described in Details. If \code{which}
#'   has length greater than one, additional plots are stepped through in order
#'   using \code{<Return>}. The default for [splm()] fitted model objects is
#'   \code{which = c(1, 2, 7)}. The default for [spautor()] fitted model objects
#'   is \code{which = c(1, 2)}.
#' @param ... Other arguments passed to other methods.
#'
#' @details For [splm()] and [spautor()], the values of \code{which} make the
#'   corresponding plot:
#'   \itemize{
#'     \item{1:}{ Standardized residuals vs fitted values (of the response)}
#'     \item{2:}{ Normal quantile-quantile plot of standardized residuals}
#'     \item{3:}{ Scale-location plot of standardized residuals}
#'     \item{4:}{ Cook's distance}
#'     \item{5:}{ Standardized residuals vs leverage}
#'     \item{6:}{ Cook's distance vs leverage}
#'   }
#'   For [splm()], there are two additional values of \code{which}:
#'   \itemize{
#'     \item{7:}{ Fitted spatial covariance function vs distance}
#'     \item{8:}{ Fitted anisotropic level curve of equal correlation}
#'   }
#'
#' @return No return value. Function called for plotting side effects.
#'
#' @name plot.spmodel
#' @method plot splm
#' @export
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' plot(spmod)
#' plot(spmod, which = c(1, 2, 4, 6))
plot.splm <- function(x, which, ...) {
  if (missing(which)) {
    which <- c(1, 2, 7)
  }

  if (any(!(which %in% 1:8))) {
    stop("Values of which can only take on 1, 2, 3, 4, 5, 6, 7, or 8.", call. = FALSE)
  }

  # setting old graphical parameter value
  oldpar <- par(no.readonly = TRUE)
  # setting exit handler
  on.exit(par(ask = oldpar$ask), add = TRUE)
  # set ask
  if (length(which) > 1) {
    par(ask = TRUE)
  }


  cal <- x$call
  if (!is.na(m.f <- match("formula", names(cal)))) {
    cal <- cal[c(1, m.f)]
    names(cal)[2L] <- ""
  }
  cc <- deparse(cal, 80)
  nc <- nchar(cc[1L], "c")
  abbr <- length(cc) > 1 || nc > 75
  sub.caption <- if (abbr) {
    paste(substr(cc[1L], 1L, min(75L, nc)), "...")
  } else {
    cc[1L]
  }


  # plot 1
  if (1 %in% which) {
    plot(
      x = fitted(x),
      y = rstandard(x),
      main = "Standardized Residuals vs Fitted",
      xlab = "Fitted values",
      ylab = "Standardized residuals",
      ...
    )
    title(sub = sub.caption)
    abline(h = 0, lty = 3, col = "gray")
  }

  # plot 2
  if (2 %in% which) {
    qqnorm(rstandard(x), main = "Normal Q-Q", ylab = "Standardized residuals", ...)
    qqline(rstandard(x), ...)
    title(sub = sub.caption)
  }


  # plot 3
  if (3 %in% which) {
    plot(
      x = fitted(x),
      y = sqrt(abs(rstandard(x))),
      main = "Scale-Location",
      xlab = "Fitted values",
      ylab = as.expression(substitute(sqrt(abs(YL)), list(YL = as.name("Standardized residuals")))),
      ...
    )
    title(sub = sub.caption)
  }




  # plot 4
  if (4 %in% which) {
    plot(
      x = seq_len(x$n),
      y = cooks.distance(x),
      xlab = "Obs. Number",
      ylab = "Cook's distance",
      main = "Cook's Distance",
      type = "h",
      ...
    )
    title(sub = sub.caption)
  }


  # plot 5
  if (5 %in% which) {
    plot(
      x = hatvalues(x),
      y = rstandard(x),
      xlab = "Leverage",
      ylab = "Standardized residuals",
      main = "Standardized residuals vs Leverage",
      ...
    )
    title(sub = sub.caption)
    abline(h = 0, v = 0, lty = 3, col = "gray")
  }

  # plot 6
  if (6 %in% which) {
    plot(
      x = hatvalues(x),
      y = cooks.distance(x),
      xlab = "Leverage",
      ylab = "Cook's distance",
      main = "Cook's dist vs Leverage",
      ...
    )
    title(sub = sub.caption)
  }


  # plot 7
  if (7 %in% which) {

    h <- seq(0, x$max_dist, length.out = 1000)
    spcoef <- coefficients(x, type = "spcov")
    spcov_vector_val <- spcov_vector(spcoef, h)
    plot(
      x = h[-1],
      y = spcov_vector_val[-1],
      xlab = "Distance",
      ylab = paste("Covariance:", class(coef(x, type = "spcov"))),
      main = "Fitted spatial covariance function",
      ylim = c(x = h[1], y = sum(spcoef[c("de", "ie")])),
      type = "l",
      ...
    )
    points(x = h[1], y = sum(spcoef[c("de", "ie")]), ...)
    title(sub = sub.caption)
  }

  if (8 %in% which) {
    r <- 1
    theta_seq <- seq(0, 2 * pi, length.out = 1000)
    x_orig <- r * cos(theta_seq)
    y_orig <- r * sin(theta_seq)
    spcoef <- coefficients(x, type = "spcov")
    rotate <- spcoef[["rotate"]]
    scale <- spcoef[["scale"]]
    if (rotate != 0 || scale != 1) {
      dat <- data.frame(x_orig = x_orig, y_orig = y_orig)
      new_coords <- transform_anis_inv(dat, "x_orig", "y_orig", rotate, scale)
      x_new <- new_coords$xcoord_val
      y_new <- new_coords$ycoord_val
    } else {
      x_new <- x_orig
      y_new <- y_orig
    }
    plot(
      x = x_new,
      y = y_new,
      xlab = "x-distance",
      ylab = "y-distance",
      main = "Anisotropic level curve", # of equal correlation
      type = "l",
      xlim = c(-1, 1),
      ylim = c(-1, 1),
      xaxt = "n", # remove axis information
      yaxt = "n", # remove axis information
      ...
    )
    title(sub = sub.caption)
  }
}

#' @rdname plot.spmodel
#' @method plot spautor
#' @export
plot.spautor <- function(x, which, ...) {
  if (missing(which)) {
    which <- c(1, 2)
  }

  if (any(!(which %in% 1:6))) {
    stop("Values of which can only take on 1, 2, 3, 4, 5, or 6.", call. = FALSE)
  }

  # setting old graphical parameter value
  oldpar <- par(no.readonly = TRUE)
  # setting exit handler
  on.exit(par(ask = oldpar$ask), add = TRUE)
  # set ask
  if (length(which) > 1) {
    par(ask = TRUE)
  }


  cal <- x$call
  if (!is.na(m.f <- match("formula", names(cal)))) {
    cal <- cal[c(1, m.f)]
    names(cal)[2L] <- ""
  }
  cc <- deparse(cal, 80)
  nc <- nchar(cc[1L], "c")
  abbr <- length(cc) > 1 || nc > 75
  sub.caption <- if (abbr) {
    paste(substr(cc[1L], 1L, min(75L, nc)), "...")
  } else {
    cc[1L]
  }


  # plot 1
  if (1 %in% which) {
    plot(
      x = fitted(x),
      y = rstandard(x),
      main = "Standardized Residuals vs Fitted",
      xlab = "Fitted values",
      ylab = "Standardized residuals",
      ...
    )
    title(sub = sub.caption)
    abline(h = 0, lty = 3, col = "gray")
  }

  # plot 2
  if (2 %in% which) {
    qqnorm(rstandard(x), main = "Normal Q-Q", ylab = "Standardized residuals", ...)
    qqline(rstandard(x), ...)
    title(sub = sub.caption)
  }


  # plot 3
  if (3 %in% which) {
    plot(
      x = fitted(x),
      y = sqrt(abs(rstandard(x))),
      main = "Scale-Location",
      xlab = "Fitted values",
      ylab = as.expression(substitute(sqrt(abs(YL)), list(YL = as.name("Standardized residuals")))),
      ...
    )
    title(sub = sub.caption)
  }




  # plot 4
  if (4 %in% which) {
    plot(
      x = seq_len(x$n),
      y = cooks.distance(x),
      xlab = "Obs. Number",
      ylab = "Cook's distance",
      main = "Cook's Distance",
      type = "h",
      ...
    )
    title(sub = sub.caption)
  }


  # plot 5
  if (5 %in% which) {
    plot(
      x = hatvalues(x),
      y = rstandard(x),
      xlab = "Leverage",
      ylab = "Standardized residuals",
      main = "Standardized residuals vs Leverage",
      ...
    )
    title(sub = sub.caption)
    abline(h = 0, v = 0, lty = 3, col = "gray")
  }

  # plot 6
  if (6 %in% which) {
    plot(
      x = hatvalues(x),
      y = cooks.distance(x),
      xlab = "Leverage",
      ylab = "Cook's distance",
      main = "Cook's dist vs Leverage",
      ...
    )
    title(sub = sub.caption)
  }
}
