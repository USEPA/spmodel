satterthwaite_anova <- function(object, ..., test = TRUE, Terms, L, method) {

  if (missing(Terms)) Terms <- NULL
  if (missing(L)) L <- NULL
  if (missing(method)) method <- NULL

  method <- get_satterthwaite_method(object, method)

  validate_satterthwaite_scope(object, method)

  L <- get_L(L, Terms, object)

  context <- get_cov_gradients_context_splm(object)
  vcov_theta <- get_vcov_theta(method, context, object)
  betahat <- coef(object, type = "fixed")

  DenDF <- unlist(lapply(L, fai_cornelius, betahat = betahat,
                 vcov_theta = vcov_theta, context = context,
                 method = method, object))

  # call anova() function and appropriate scale quantities to be F instead of X2
  anova_val <- do.call(rbind, lapply(L, get_marginal_Chi2, object))
  anova_val$`F value`<- anova_val$Chi2/anova_val$Df
  anova_val$NumDF <- anova_val$Df
  anova_val$DenDF <- DenDF
  anova_val <- anova_val[, c("NumDF", "DenDF", "F value")]
  anova_val$`Pr(>F)` <- pf(anova_val$`F value`, anova_val$NumDF, anova_val$DenDF, lower.tail = FALSE)
  anova_val
}

fai_cornelius <- function(Lsub, betahat, vcov_theta, context, method, object) {

  if (!is.matrix(Lsub)) Lsub <- matrix(Lsub, nrow = 1)
  q <- nrow(Lsub)
  C <- Lsub %*% tcrossprod(vcov(object), Lsub)
  eig <- eigen(C)
  U <- eig$vectors
  eta <- crossprod(U, Lsub)

  nu <- lapply(seq_len(q), function(i) {
    eta_i <- eta[i, ]
    g <- as.numeric(crossprod(eta_i, vcov(object)) %*% eta_i)
    grad_g <- get_grad_g(eta_i, method, context, object)
    satterthwaite_df <- get_satterthwaite_df(g, grad_g, vcov_theta)
  })

  nu <- unlist(nu)

  if (any(!is.finite(nu)) || any(nu <= 2)) {
    DenDF <- NA
  } else {
    E <- sum(nu / (nu - 2))
    DenDF <- if (E > q) 2 * E / (E - q) else NA
  }
  DenDF
}
