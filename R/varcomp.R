#' Variability component comparison
#'
#' @description Compare the proportion of total variability explained by the fixed effects
#'   and each variance parameter.
#'
#' @param object A fitted model object from [splm()] or [spautor()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return A tibble that partitions the the total variability by the fixed effects
#'   and each variance parameter. The proportion of variability explained by the
#'   fixed effects is the pseudo R-squared obtained by \code{psuedoR2()}. The
#'   remaining proportion is spread accordingly among each variance parameter:
#'   \code{"de"}, \code{"ie"}, and if random effects are used, each named random effect.
#'   If \code{spautor()} objects have unconnected sites, a list is returned with three elements:
#'   \code{"connected"} for a variability comparison among the connected sites;
#'   \code{"unconnected"} for a variability comparison among the unconnected
#'   sites; and \code{"ratio"} for the ratio of the variance of the connected
#'   sites relative to the variance of the unconnected sites.
#'
#' @export
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' varcomp(spmod)
varcomp <- function(object, ...) {
  UseMethod("varcomp", object)
}

#' @rdname varcomp
#' @method varcomp spmod
#' @export
varcomp.spmod <- function(object, ...) {

  PR2 <- pseudoR2(object)
  spcov_coef <- coef(object, type = "spcov")
  de <- spcov_coef[["de"]]
  ie <- spcov_coef[["ie"]]
  randcov_coef <- coef(object, type = "randcov")
  total_var <- sum(de, ie, randcov_coef)
  varcomp_names <- c("Covariates (PR-sq)", "de", "ie", c(names(randcov_coef)))
  varcomp_values <- c(PR2, (1 - PR2) * c(de, ie, randcov_coef) / total_var)
  if (object$fn == "spautor" && spcov_coef[["extra"]] != 0) {
    varcomp_val1 <- tibble::tibble(varcomp = varcomp_names, proportion = varcomp_values)
    extra <- spcov_coef[["extra"]]
    total_var2 <- sum(extra, randcov_coef)
    varcomp_names2 <- c("Covariates (PR-sq)", "extra", c(names(randcov_coef)))
    varcomp_values2 <- c(PR2, (1 - PR2) * c(extra, randcov_coef) / total_var2)
    varcomp_val2 <- tibble::tibble(varcomp = varcomp_names2, proportion = varcomp_values2)
    ratio <- total_var / total_var2
    varcomp_val <- list(connected = varcomp_val1, unconnected = varcomp_val2, ratio = ratio)
  } else {
    varcomp_val <- tibble::tibble(varcomp = varcomp_names, proportion = varcomp_values)
  }
  varcomp_val
}

