#' @rdname tidy.spmodel
#' @method tidy spglm
#' @order 3
#' @export
tidy.spglm <- function(x, conf.int = FALSE,
                       conf.level = 0.95, effects = "fixed", ...) {

  if (conf.int && (conf.level < 0 || conf.level > 1)) {
    stop("conf.level must be between 0 and 1.", call. = FALSE)
  }

  if (effects == "fixed") {
    result <- tibble::as_tibble(summary(x)$coefficients$fixed,
      rownames = "term", .name_repair = "minimal"
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
      rownames = "term", .name_repair = "minimal"
      )
      colnames(ci) <- c("term", "conf.low", "conf.high")
      result <- tibble::as_tibble(base::merge(result, ci, by = "term"), .name_repair = "minimal")
    }
  } else if (effects == "spcov") {
    spcoef <- coefficients(x, type = "spcov")
    result <- tibble::as_tibble(unclass(spcoef),
      rownames = "term", .name_repair = "minimal"
    )
    colnames(result) <- c("term", "estimate")
    result$is_known <- x$is_known$spcov

    if (!x$anisotropy) {
      which_rotate <- which(result$term == "rotate")
      which_scale <- which(result$term == "scale")
      result <- result[-c(which_rotate, which_scale), , drop = FALSE]
    }
    if (inherits(spcoef, c("none", "ie"))) {
      which_ie <- which(result$term == "ie")
      result <- result[which_ie, , drop = FALSE]
    }
  } else if (effects == "dispersion") {
    result <- tibble::as_tibble(unclass(summary(x)$coefficients$dispersion),
      rownames = "term", .name_repair = "minimal"
    )
    colnames(result) <- c("term", "estimate")
    result$is_known <- x$is_known$dispersion
  } else if (effects == "randcov") {
    if (is.null(summary(x)$coefficients$randcov)) {
      result <- NULL
    } else {
      result <- tibble::as_tibble(summary(x)$coefficients$randcov,
        rownames = "term", .name_repair = "minimal"
      )
      colnames(result) <- c("term", "estimate")
      result$is_known <- x$is_known$randcov
    }
  }

  result
}

#' @rdname tidy.spmodel
#' @method tidy spgautor
#' @order 4
#' @export
tidy.spgautor <- function(x, conf.int = FALSE,
                          conf.level = 0.95, effects = "fixed", ...) {

  if (conf.int && (conf.level < 0 || conf.level > 1)) {
    stop("conf.level must be between 0 and 1.", call. = FALSE)
  }

  if (effects == "fixed") {
    result <- tibble::as_tibble(summary(x)$coefficients$fixed,
      rownames = "term", .name_repair = "minimal"
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
      rownames = "term", .name_repair = "minimal"
      )
      colnames(ci) <- c("term", "conf.low", "conf.high")
      result <- tibble::as_tibble(base::merge(result, ci, by = "term"), .name_repair = "minimal")
    }
  } else if (effects == "spcov") {
    spcoef <- coefficients(x, type = "spcov")
    result <- tibble::as_tibble(unclass(spcoef),
      rownames = "term", .name_repair = "minimal"
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
  } else if (effects == "dispersion") {
    result <- tibble::as_tibble(unclass(summary(x)$coefficients$dispersion),
      rownames = "term", .name_repair = "minimal"
    )
    colnames(result) <- c("term", "estimate")
    result$is_known <- x$is_known$dispersion
  } else if (effects == "randcov") {
    if (is.null(summary(x)$coefficients$randcov)) {
      result <- NULL
    } else {
      result <- tibble::as_tibble(summary(x)$coefficients$randcov,
        rownames = "term", .name_repair = "minimal"
      )
      colnames(result) <- c("term", "estimate")
      result$is_known <- x$is_known$randcov
    }
  }

  result
}
