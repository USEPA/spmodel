# recover data
# do NOT export or document (see zzz.R)
recover_data.splm <- function(object, frame = model.frame(object), ...) {
  # check to see if emmeans installed
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop("Install the emmeans package before using", call. = FALSE)
  }
  # recover data (using emmeans code)
  fcall = object$call
  # recognize that lm objects have a $model element that is model.frame(object)
  emmeans::recover_data(fcall, delete.response(terms(object)), frame = frame, na.action = NULL, ...)
}

recover_data.spautor <- recover_data.splm

recover_data.spglm <- recover_data.splm

recover_data.spgautor <- recover_data.splm



# get emm basis
# do NOT export or document (see zzz.R)
emm_basis.splm <- function(object, trms, xlev, grid, ...) {
  # emm_basis
  # check to see if emmeans installed
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop("Install the emmeans package before using", call. = FALSE)
  }
  bhat = coef(object)
  nm = if (is.null(names(bhat)))
    row.names(bhat)
  else names(bhat)
  m = suppressWarnings(model.frame(trms, grid, na.action = na.pass,
                                   xlev = xlev))
  X = model.matrix(trms, m, contrasts.arg = object$contrasts)
  assign = attr(X, "assign")
  X = X[, nm, drop = FALSE]
  bhat = as.numeric(bhat)
  V = emmeans::.my.vcov(object, ...)
  nbasis = estimability::all.estble # returns a 1x1 NA which says all functions estimable
  dfargs = misc = list()
  dffun = function(k, dfargs) Inf
  attr(dffun, "mesg") = "asymptotic"
  mm <- model.matrix(object)
  mm = emmeans::.cmpMM(mm, assign = attr(mm, "assign"))
  if (inherits(object, c("spglm", "spgautor"))) {
    famdat <- data.frame(
      family = c("poisson", "nbinomial", "binomial", "beta", "Gamma", "inverse.gaussian")
    )
    famdat$link <- ifelse(famdat$family %in% c("binomial", "beta"), "logit", "log")
    fam = famdat[match(object$family, famdat$family), ]
    fam = list(family = fam$family, link = fam$link)
    misc = emmeans::.std.link.labels(fam, misc)
  }
  list(X = X, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun,
       dfargs = dfargs, misc = misc, model.matrix = mm)
}

emm_basis.spautor <- emm_basis.splm

emm_basis.spglm <- emm_basis.splm

emm_basis.spgautor <- emm_basis.splm
