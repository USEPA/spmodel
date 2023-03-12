#' Tidy a fitted model object
#'
#' @description Tidy a fitted model object into a summarized tibble.
#'
#' @param x A fitted model object from [splm()] or [spautor()]
#' @param conf.int Logical indicating whether or not to include a confidence interval
#'   in the tidied output. The default is \code{FALSE}.
#' @param conf.level The confidence level to use for the confidence interval if
#'   \code{conf.int} is \code{TRUE}. Must be strictly greater than 0 and less than 1.
#'   The default is 0.95, which corresponds to a 95 percent confidence interval.
#' @param effects The type of effects to tidy. Available options are \code{"fixed"}
#'   (fixed effects), \code{"spcov"} (spatial covariance parameters),
#'   and \code{"randcov"} (random effect variances). The default is \code{"fixed"}.
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return A tidy tibble of summary information \code{effects}.
#'
#' @name tidy.spmodel
#' @method tidy splm
#' @order 1
#' @export
#'
#' @seealso [glance.spmodel()] [augment.spmodel()]
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' tidy(spmod)
#' tidy(spmod, effects = "spcov")
tidy.splm <- function(x, conf.int = FALSE,
                       conf.level = 0.95, effects = "fixed", ...) {
  if (effects == "fixed") {
    result <- tibble::as_tibble(summary(x)$coefficients$fixed,
      rownames = "term"
    )
    colnames(result) <- c(
      "term", "estimate", "std.error",
      "statistic", "p.value"
    )

    if (conf.int) {
      ci <- tibble::as_tibble(confint(x,
        level = conf.level,
        type = "fixed"
      ),
      rownames = "term"
      )
      colnames(ci) <- c("term", "conf.low", "conf.high")
      result <- tibble::as_tibble(base::merge(result, ci, by = "term"))
    }
  } else if (effects == "spcov") {
    spcoef <- coefficients(x, type = "spcov")
    result <- tibble::as_tibble(unclass(spcoef),
      rownames = "term"
    )
    colnames(result) <- c("term", "estimate")
    result$is_known <- x$is_known$spcov

    if (!x$anisotropy) {
      which_rotate <- which(result$term == "rotate")
      which_scale <- which(result$term == "scale")
      result <- result[-c(which_rotate, which_scale), , drop = FALSE]
    }
    if (inherits(spcoef, "none")) {
      which_ie <- which(result$term == "ie")
      result <- result[which_ie, , drop = FALSE]
    }

  } else if (effects == "randcov") {
    if (is.null(summary(x)$coefficients$randcov)) {
      result <- NULL
    } else {
      result <- tibble::as_tibble(summary(x)$coefficients$randcov,
        rownames = "term"
      )
      colnames(result) <- c("term", "estimate")
      result$is_known <- x$is_known$randcov
    }
  }

  result
}

#' @rdname tidy.spmodel
#' @method tidy spautor
#' @order 2
#' @export
tidy.spautor <- function(x, conf.int = FALSE,
                      conf.level = 0.95, effects = "fixed", ...) {
  if (effects == "fixed") {
    result <- tibble::as_tibble(summary(x)$coefficients$fixed,
                                rownames = "term"
    )
    colnames(result) <- c(
      "term", "estimate", "std.error",
      "statistic", "p.value"
    )

    if (conf.int) {
      ci <- tibble::as_tibble(confint(x,
                                      level = conf.level,
                                      type = "fixed"
      ),
      rownames = "term"
      )
      colnames(ci) <- c("term", "conf.low", "conf.high")
      result <- tibble::as_tibble(base::merge(result, ci, by = "term"))
    }
  } else if (effects == "spcov") {
    spcoef <- coefficients(x, type = "spcov")
    result <- tibble::as_tibble(unclass(spcoef),
                                rownames = "term"
    )
    colnames(result) <- c("term", "estimate")
    result$is_known <- x$is_known$spcov


    no_ie <- spcoef[["ie"]] == 0 && x$is_known$spcov[["ie"]]
    if (no_ie) {
      which_ie <- which(result$term == "ie")
      result <- result[-c(which_ie), , drop = FALSE]
    }
    no_extra <- spcoef[["extra"]] == 0 && x$is_known$spcov[["extra"]]
    if (no_extra) {
      which_extra <- which(result$term == "extra")
      result <- result[-c(which_extra), , drop = FALSE]
    }

  } else if (effects == "randcov") {
    if (is.null(summary(x)$coefficients$randcov)) {
      result <- NULL
    } else {
      result <- tibble::as_tibble(summary(x)$coefficients$randcov,
                                  rownames = "term"
      )
      colnames(result) <- c("term", "estimate")
      result$is_known <- x$is_known$randcov
    }
  }

  result
}
