#' Package load hook
#'
#' @param libname The library directory
#' @param pkgname The package name
#'
#' @return Nothing; if \code{emmeans} is installed, dynamically registers
#'   \code{splm}, \code{spautor}, \code{spglm}, and \code{spgautor} support
#'   for it (see \code{recover_data.splm()} and \code{emm_basis.splm()})
#'
#' @noRd
.onLoad <- function(libname, pkgname) {
  # emmeans is a Suggests (not Imports) dependency, so emmeans support is
  # only wired up conditionally, and only when emmeans is actually installed
  if (requireNamespace("emmeans", quietly = TRUE)) {
    # suggestion from glmmTMB zzz.R for dynamically loading emmeans
    # can use utils:: because it is part of base R (even though it is not in
    # suggests)
    # emmeans::.emm_register() is only available starting in emmeans 1.4, so
    # older installed versions are rejected with an informative error instead
    # of failing later with an obscure "could not find function" error
    if (utils::packageVersion("emmeans") < "1.4") {
      stop("please install a newer version of emmeans (> 1.4)", call. = FALSE)
    }

    # registers recover_data.*() / emm_basis.*() S3 methods for spmodel's
    # fitted model classes so emmeans::emmeans() etc. work on them without
    # spmodel needing to formally Import/Depend on emmeans
    emmeans::.emm_register(c("splm", "spautor", "spglm", "spgautor"), pkgname)
  }
  ## https://stackoverflow.com/questions/49056642/how-to-make-variable-available-to-namespace-at-loading-time/
}
