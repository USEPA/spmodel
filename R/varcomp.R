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
#' @order 1
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
#' @method varcomp splm
#' @order 2
#' @export
varcomp.splm <- function(object, ...) {

  PR2 <- pseudoR2(object)
  spcov_coef <- coef(object, type = "spcov")
  de <- spcov_coef[["de"]]
  ie <- spcov_coef[["ie"]]
  randcov_coef <- coef(object, type = "randcov")
  total_var <- sum(de, ie, randcov_coef)
  varcomp_names <- c("Covariates (PR-sq)", "de", "ie", c(names(randcov_coef)))
  varcomp_values <- c(PR2, (1 - PR2) * c(de, ie, randcov_coef) / total_var)
  varcomp_val <- tibble::tibble(varcomp = varcomp_names, proportion = varcomp_values)
  varcomp_val
}

#' @rdname varcomp
#' @method varcomp spautor
#' @order 3
#' @export
varcomp.spautor <- function(object, ...) {

  PR2 <- pseudoR2(object)
  spcov_coef <- coef(object, type = "spcov")
  de <- spcov_coef[["de"]]
  ie <- spcov_coef[["ie"]]
  extra <- spcov_coef[["extra"]]
  randcov_coef <- coef(object, type = "randcov")
  total_var_con <- sum(de, ie, randcov_coef)
  varcomp_names_con <- c("Covariates (PR-sq)", "de", "ie", c(names(randcov_coef)))
  varcomp_values_con <- c(PR2, (1 - PR2) * c(de, ie, randcov_coef) / total_var_con)
  varcomp_val <- tibble::tibble(varcomp = varcomp_names_con, proportion = varcomp_values_con)
  if (extra != 0) {
    total_var_uncon <- sum(extra, randcov_coef)
    varcomp_names_uncon <- c("Covariates (PR-sq)", "extra", c(names(randcov_coef)))
    varcomp_values_uncon <- c(PR2, (1 - PR2) * c(extra, randcov_coef) / total_var_uncon)
    varcomp_val_2 <- tibble::tibble(varcomp = varcomp_names_uncon, proportion = varcomp_values_uncon)
    ratio <- total_var_con / total_var_uncon
    varcomp_val <- list(connected = varcomp_val, unconnected = varcomp_val_2, ratio = ratio)
  }
  varcomp_val
}
