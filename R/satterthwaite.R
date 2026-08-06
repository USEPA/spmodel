#' Title
#'
#' @param object filler
#' @param ... filler
#'
#' @returns filler
#' @export
satterthwaite <- function(object, ...) {
  UseMethod("satterthwaite")
}

#' @method satterthwaite splm
#' @order 2
#' @export
satterthwaite.splm <- function(object, method, ...) {

  if (missing(method)) method <- NULL
  method <- get_satterthwaite_method(object, method)

  validate_satterthwaite_scope(object, method)
  context <- get_cov_gradients_context_splm(object)

  vcov_theta <- get_vcov_theta(method, context, object)

  X <- context$X
  p <- ncol(X)
  betahat <- coef(object, type = "fixed")

  L <- diag(p)
  colnames(L) <- names(betahat)
  rownames(L) <- names(betahat)

  ddf <- lapply(seq_len(p), function(i) {
    Li <- L[i, ]
    g <- as.numeric(crossprod(Li, vcov(object)) %*% Li)
    grad_g <- get_grad_g(Li, method, context, object)
    satterthwaite_df <- get_satterthwaite_df(g, grad_g, vcov_theta)
  })

  ddf <- unlist(ddf)
  names(ddf) <- names(betahat)
  ddf


}

get_satterthwaite_method <- function(object, method) {
  if (is.null(method)) {
    if (inherits(coef(object, type = "spcov"), "exponential") && !isTRUE(object$anisotropy)) {
      method <- "closed"
    } else {
      method <- "numeric"
    }
  }
  if (method == "closed" && (!inherits(coef(object, type = "spcov"), "exponential") || isTRUE(object$anisotropy))) {
    method <- "numeric"
    warning("Closed form not available for this model. Using method = \"numeric\".", call. = FALSE)
  }
  method
}

get_satterthwaite_df <- function(g, grad_g, vcov_theta) {
  as.numeric(2 * g^2 / (crossprod(grad_g, vcov_theta) %*% grad_g))
}