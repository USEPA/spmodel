.onLoad <- function(libname, pkgname) {
  if (requireNamespace("emmeans", quietly = TRUE)) {


    # suggestion from glmmTMB zzz.R for dynamically loading emmeans
    # can use utils:: because it is part of base R (even though it is not in
    # suggests)
    if (utils::packageVersion("emmeans") < "1.4") {
      stop("please install a newer version of emmeans (> 1.4)", call. = FALSE)
    }

    emmeans::.emm_register(c("splm", "spautor", "spglm", "spgautor"), pkgname)
  }
  ## https://stackoverflow.com/questions/49056642/how-to-make-variable-available-to-namespace-at-loading-time/
}
