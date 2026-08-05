satterthwaite_anova <- function(object, ..., test = TRUE, Terms, L, method) {

  if (missing(method)) method <- NULL
  method <- get_satterthwaite_method(object, method)

  if (method == "numeric" && !requireNamespace("numDeriv", quietly = TRUE)) {
    stop("Install the numDeriv package before using satterthwaite(method = \"numeric\").", call. = FALSE)
  }
  if (!method %in% c("closed", "numeric")) stop("method must be \"closed\" or \"numeric\".", call. = FALSE)

  validate_satterthwaite_scope(object)

  # build a hypothesis matrix L (or list of them) and run a general linear
  # hypothesis test (GLHT) L*beta = 0 for each set of terms
  if (missing(L)) {
    # "assign" attribute maps each column of the model matrix to the model
    # term that generated it (0 = intercept), used to group coefficients
    # belonging to the same term (e.g. all dummy columns of a factor)
    assign_indices <- attr(model.matrix(object), "assign") + 1
    # attr(model.matrix(object), "assign") if centering at zero
    if (missing(Terms)) {
      # default: test each term separately (type III / marginal tests)
      assign_index <- unique(assign_indices)
      L <- lapply(assign_index, get_L_list, assign_indices)
      label <- labels(object)
      if (attr(terms(object), "intercept") == 1) {
        label <- c("(Intercept)", label)
      }
      names(L) <- label
    } else {
      # Terms specified: build one L testing the listed terms jointly
      if (is.character(Terms)) {
        Terms <- which(c("(Intercept)", labels(object)) %in% Terms) # - 1 if centering at zero
      }
      L <- list(do.call(rbind, lapply(Terms, get_L_list, assign_indices)))
      label <- c("(Intercept)", labels(object))
      label <- label[Terms] # label[Terms + 1] if centering at zero
      names(L) <- paste(label, collapse = ", ")
    }
  } else {
    # user supplied custom contrast matrix/matrices directly
    if (!is.list(L)) {
      L <- list(L)
    }
    names(L) <- paste("contrast", seq_along(L), sep = "")
  }

  context <- get_satterthwaite_context_splm(object)
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
    satterthwaite_ddf <- as.numeric(2 * g^2 / (crossprod(grad_g, vcov_theta) %*% grad_g))
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
