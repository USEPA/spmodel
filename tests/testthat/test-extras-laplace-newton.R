skip_on_cran()
skip_if_not(identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"))

# Frozen spmodel 0.15.0 functions keep the comparison independent of the new helpers.
v0.15.0_reference <- new.env(parent = asNamespace("spmodel"))
v0.15.0_reference$get_w_and_H_spglm <- 
function (data_object, dispersion, SigInv_list, SigInv_X, cov_betahat, cov_betahat_Inv, estmethod, ret_mHInv = FALSE) 
{
    family <- data_object$family
    SigInv <- Matrix::bdiag(SigInv_list)
    Ptheta <- SigInv - SigInv_X %*% tcrossprod(cov_betahat, SigInv_X)
    y <- as.vector(do.call("rbind", data_object$y_list))
    size <- data_object$size
    w <- get_w_init(family, y, dispersion)
    wdiffmax <- Inf
    iter <- 0
    offset <- if (is.null(data_object$offset)) 
        0
    else as.vector(data_object$offset)
    if (length(SigInv_list) == 1) {
        while (iter < 50 && wdiffmax > 1e-04) {
            iter <- iter + 1
            d <- get_d(family, w + offset, y, size, dispersion)
            g <- d - Ptheta %*% w
            D <- get_D(family, w + offset, y, size, dispersion)
            H <- D - Ptheta
            solveHg <- solve(H, g)
            wnew <- w - solveHg
            dnew <- get_d(family, wnew + offset, y, size, dispersion)
            gnew <- dnew - Ptheta %*% wnew
            if (any(is.na(gnew) | is.infinite(gnew))) 
                stop("Convergence problem. Try using a different family, removing extreme observations, rescaling the response variable (if continuous), fixing ie at a known, non-zero value (via spcov_initial), or fixing dispersion at one (via dispersion_initial).", 
                  call. = FALSE)
            if (max(abs(gnew)) > max(abs(g))) 
                wnew <- w - 0.1 * solveHg
            wdiffmax <- max(abs(wnew - w))
            w <- wnew
        }
        mHldet <- as.numeric(determinant(if (estmethod == "ml") SigInv - D else -H, logarithm = TRUE)$modulus)
        w_and_H_list <- list(w = w, H = NULL, mHldet = mHldet)
        if (ret_mHInv) {
            HInv <- solve(H)
            w_and_H_list$mHInv <- -HInv
        }
    }
    else {
        diag(cov_betahat_Inv) <- diag(cov_betahat_Inv) + data_object$diagtol
        while (iter < 50 && wdiffmax > 1e-04) {
            iter <- iter + 1
            d <- get_d(family, w + offset, y, size, dispersion)
            g <- d - Ptheta %*% w
            D <- get_D(family, w + offset, y, size, dispersion)
            D_diag <- diag(D)
            D_list <- lapply(split(D_diag, sort(data_object$local_index)), function(x) Diagonal(x = x))
            if (data_object$parallel) {
                cluster_list <- lapply(seq_along(D_list), function(l) {
                  cluster_list_element <- list(D = D_list[[l]], S = SigInv_list[[l]])
                })
                DSigInv_list <- parallel::parLapply(data_object$cl, cluster_list, get_DSigInv_parallel)
                names(DSigInv_list) <- names(D_list)
            }
            else {
                DSigInv_list <- mapply(D = D_list, S = SigInv_list, function(D, S) get_DSigInv(D, S), SIMPLIFY = FALSE)
            }
            if (data_object$parallel) {
                cluster_list <- DSigInv_list
                DSigInv_Inv_list <- parallel::parLapply(data_object$cl, cluster_list, solve)
                names(DSigInv_Inv_list) <- names(D_list)
            }
            else {
                DSigInv_Inv_list <- lapply(DSigInv_list, function(x) solve(x))
            }
            DSigInv_Inv <- Matrix::bdiag(DSigInv_Inv_list)
            HInv <- smw_HInv(AInv = DSigInv_Inv, U = SigInv_X, CInv = cov_betahat_Inv)
            solveHg <- HInv %*% g
            wnew <- w - solveHg
            dnew <- get_d(family, wnew + offset, y, size, dispersion)
            gnew <- dnew - Ptheta %*% wnew
            if (any(is.na(gnew) | is.infinite(gnew))) 
                stop("Convergence problem. Try using a different family, removing extreme observations, rescaling the response variable (if continuous), fixing ie at a known, non-zero value (via spcov_initial), or fixing dispersion at one (via dispersion_initial).", 
                  call. = FALSE)
            if (max(abs(gnew)) > max(abs(g))) 
                wnew <- w - 0.1 * solveHg
            wdiffmax <- max(abs(wnew - w))
            w <- wnew
        }
        if (estmethod == "ml") {
            mHldet <- sum(vapply(DSigInv_list, function(x) {
                as.numeric(determinant(x, logarithm = TRUE)$modulus)
            }, numeric(1)))
        }
        else {
            mHldet <- smw_mHldet(A_list = DSigInv_list, AInv = DSigInv_Inv, U = SigInv_X, C = cov_betahat, CInv = cov_betahat_Inv)
        }
        w_and_H_list <- list(w = w, H = NULL, mHldet = mHldet)
        if (ret_mHInv) {
            w_and_H_list$mHInv <- -HInv
        }
    }
    w_and_H_list
}
environment(v0.15.0_reference$get_w_and_H_spglm) <- v0.15.0_reference

v0.15.0_reference$smw_HInv <- 
function (AInv, U, CInv) 
{
    mid <- CInv + t(U) %*% AInv %*% U
    AInv - (AInv %*% U) %*% solve(mid) %*% (t(U) %*% AInv)
}
environment(v0.15.0_reference$smw_HInv) <- v0.15.0_reference

v0.15.0_reference$smw_mHldet <- 
function (A_list, AInv, U, C, CInv) 
{
    Aldet <- sum(unlist(lapply(A_list, function(x) determinant(x, logarithm = TRUE)$modulus)))
    Cldet <- 2 * sum(log(diag(t(chol(C)))))
    mid <- CInv + t(U) %*% AInv %*% U
    midldet <- determinant(mid, logarithm = TRUE)$modulus
    as.numeric(Aldet + Cldet + midldet)
}
environment(v0.15.0_reference$smw_mHldet) <- v0.15.0_reference

v0.15.0_reference$get_DSigInv <- 
function (D, SigInv) 
{
    D - SigInv
}
environment(v0.15.0_reference$get_DSigInv) <- v0.15.0_reference

v0.15.0_reference$get_DSigInv_parallel <- 
function (cluster_list) 
{
    D <- cluster_list$D
    S <- cluster_list$S
    get_DSigInv(D, S)
}
environment(v0.15.0_reference$get_DSigInv_parallel) <- v0.15.0_reference

v0.15.0_reference$get_laplace_correction <- 
function (mHInv, family, eta, y, size, dispersion) 
{
    v <- Matrix::diag(mHInv)
    t3 <- get_D3(family, eta, y, size, dispersion)
    delta <- 0.5 * as.numeric(mHInv %*% (t3 * v))
    list(delta = delta)
}
environment(v0.15.0_reference$get_laplace_correction) <- v0.15.0_reference


v0.15.0_data <- function(n = 48L) {
  set.seed(2)
  d <- data.frame(cx = runif(n), cy = runif(n), x = rnorm(n), x2 = rnorm(n),
    x3 = 1e3 * rnorm(n), off = runif(n, -.15, .15),
    part = factor(rep(c('z', 'b'), length.out = n)), g = factor(rep(1:8, length.out = n)))
  eta <- .2 + .25 * d$x + rnorm(n, sd = .25)
  mu <- exp(eta); pr <- plogis(eta)
  d$poisson <- rpois(n, mu)
  d$nbinomial <- rnbinom(n, mu = mu, size = 4)
  d$binomial <- rbinom(n, 5, pr)
  d$bernoulli <- rbinom(n, 1, pr)
  d$Gamma <- rgamma(n, shape = 4, rate = 4 / mu)
  d$inverse.gaussian <- mu * exp(rnorm(n, sd = .25))
  d$beta <- .2 + .6 * pr
  d
}

v0.15.0_cases <- function() {
  d <- v0.15.0_data()
  out <- list()
  for (fam in c('poisson', 'binomial', 'bernoulli', 'nbinomial', 'Gamma', 'inverse.gaussian', 'beta')) {
    for (layout in c('single', 'local', 'partition')) for (method in c('ml', 'reml')) {
      out[[paste(fam, layout, method, sep = '-')]] <- v0.15.0_args(d, fam, layout, method)
    }
  }
  for (method in c('ml', 'reml')) {
    a <- v0.15.0_args(d, 'poisson', 'local', method)
    a$formula <- poisson ~ x + x2 + x3 + offset(off)
    a$local <- list(index = sample(rep(c(7, 2, 9, 4), c(1, 7, 15, 25))))
    out[[paste0('unequal-singleton-scaled-', method)]] <- a
    a <- v0.15.0_args(d, 'poisson', 'local', method)
    a$formula <- poisson ~ x + offset(off)
    a$random <- ~g
    a$randcov_initial <- randcov_initial(g = .1, known = 'given')
    a$partition_factor <- ~part
    out[[paste0('random-partition-offset-', method)]] <- a
    a <- v0.15.0_args(v0.15.0_data(240L), 'bernoulli', 'local', method)
    a$local <- list(index = rep(c(7, 2, 9, 4), each = 60))
    a$spcov_initial <- spcov_initial('exponential', de = .8, ie = 1e-4, range = .3, known = 'given')
    out[[paste0('bernoulli-floor-', method)]] <- a
  }
  out
}

v0.15.0_args <- function(d, fam, layout, method) {
  response <- if (fam == 'binomial') 'cbind(binomial, 5 - binomial)' else fam
  f <- as.formula(paste(response, '~ x'))
  a <- list(formula = f, data = d, family = if (fam == 'bernoulli') 'binomial' else fam,
    xcoord = 'cx', ycoord = 'cy', estmethod = method,
    spcov_initial = spcov_initial('exponential', de = .8, ie = .2, range = .3, known = 'given'),
    local = if (layout == 'local') list(index = sample(rep(c(7, 2, 9, 4), length.out = nrow(d)))) else FALSE)
  if (layout == 'partition') a$partition_factor <- ~part
  if (!fam %in% c('poisson', 'binomial', 'bernoulli'))
    a$dispersion_initial <- dispersion_initial(as.character(fam), 4, known = 'given')
  a
}

v0.15.0_capture <- function(args, solver = get_w_and_H_spglm) {
  force(solver)
  captured <- NULL
  fit <- testthat::with_mocked_bindings(do.call(spglm, args),
    get_w_and_H_spglm = function(data_object, dispersion, SigInv_list, SigInv_X,
      cov_betahat, cov_betahat_Inv, estmethod, ret_mHInv = FALSE) {
      a <- as.list(environment())
      a$captured <- a$solver <- a$args <- a$fit <- NULL
      if (ret_mHInv) captured <<- a
      do.call(solver, a)
    }, .package = 'spmodel')
  list(fit = fit, args = captured)
}

v0.15.0_outputs <- function(fit, d) {
  nd <- d[1:3, , drop = FALSE]
  list(coefficients = fit$coefficients, logLik = as.numeric(logLik(fit)), AIC = AIC(fit),
    deviance = deviance(fit), vcov = vcov(fit),
    link = fitted(fit, type = 'link'), response = fitted(fit, type = 'response'),
    laplace = fit$laplace, prediction = predict(fit, nd, local = FALSE, se.fit = TRUE))
}

v0.15.0_scaled_error <- function(new, ref) max(abs(as.numeric(new) - as.numeric(ref))) / max(1, max(abs(as.numeric(ref))))

v0.15.0_solve <- function(solver, a) {
  force(solver)
  calls <- 0L
  last_D <- NULL
  original_d <- get_d
  original_D <- get_D
  value <- testthat::with_mocked_bindings(do.call(solver, a),
    get_d = function(...) { calls <<- calls + 1L; original_d(...) },
    get_D = function(...) { last_D <<- original_D(...); last_D }, .package = 'spmodel')
  list(value = value, iterations = calls / 2L, D = last_D)
}

v0.15.0_check_solver <- function(a, solver = get_w_and_H_spglm) {
  old <- v0.15.0_solve(v0.15.0_reference$get_w_and_H_spglm, a)
  new <- v0.15.0_solve(solver, a)
  expect_identical(new$iterations, old$iterations)
  expect_lte(v0.15.0_scaled_error(new$value$w, old$value$w), 1e-10)
  expect_lte(v0.15.0_scaled_error(new$value$mHldet, old$value$mHldet), 1e-8)
  if (isTRUE(a$ret_mHInv)) {
    M <- new$value$mHInv
    n <- length(new$value$w)
    expect_lte(v0.15.0_scaled_error(mhinv_apply(M, diag(n)), old$value$mHInv), 1e-10)
    expect_lte(v0.15.0_scaled_error(mhinv_diag(M), Matrix::diag(old$value$mHInv)), 1e-10)
    V <- cbind(seq_along(new$value$w) / length(new$value$w), 1)
    expect_lte(v0.15.0_scaled_error(mhinv_apply(M, V), old$value$mHInv %*% V), 1e-10)
    expect_lte(v0.15.0_scaled_error(mhinv_apply(M, V[, 1]), old$value$mHInv %*% V[, 1]), 1e-10)
    eta <- as.vector(new$value$w) + if (is.null(a$data_object$offset)) 0 else a$data_object$offset
    y <- as.vector(do.call('rbind', a$data_object$y_list))
    delta <- get_laplace_correction(M, a$data_object$family, eta, y, a$data_object$size, a$dispersion)
    ref <- v0.15.0_reference$get_laplace_correction(old$value$mHInv, a$data_object$family, eta, y,
      a$data_object$size, a$dispersion)
    expect_lte(v0.15.0_scaled_error(delta$delta, ref$delta), 1e-10)
    W <- as.matrix(tcrossprod(a$cov_betahat, a$SigInv_X))
    expect_lte(v0.15.0_scaled_error(W %*% mhinv_apply(M, t(W)), W %*% old$value$mHInv %*% t(W)), 1e-10)
    if (length(a$SigInv_list) > 1L) {
      expect_s3_class(M, 'laplace_mhinv')
      expect_named(M, c('F', 'idx', 'AU', 'mid_inv'))
      expect_equal(dim(M$AU), c(n, ncol(a$SigInv_X)))
    }
  }
  invisible(list(old = old, new = new))
}

test_that('block Newton solves preserve all families, methods, and model terms', {
  set.seed(2)
  for (args in v0.15.0_cases()) {
    captured <- v0.15.0_capture(args)
    v0.15.0_check_solver(captured$args)
    a <- captured$args
    ends <- cumsum(vapply(a$SigInv_list, nrow, 1L))
    idx <- Map(seq.int, c(1L, head(ends, -1L) + 1L), ends)
    r <- seq_len(nrow(a$SigInv_X)) / 10
    block <- sum(mapply(function(S, i) as.numeric(crossprod(r[i], S %*% r[i])), a$SigInv_list, idx))
    dense <- as.numeric(crossprod(r, Matrix::bdiag(a$SigInv_list)) %*% r)
    expect_lte(v0.15.0_scaled_error(block, dense), 1e-8)
  }
})

test_that('block Newton factorization preserves parallel fitting', {
  set.seed(2)
  args <- v0.15.0_args(v0.15.0_data(), 'poisson', 'local', 'reml')
  serial <- v0.15.0_capture(args)
  cl <- parallel::makePSOCKcluster(2L)
  on.exit(parallel::stopCluster(cl))
  parallel::clusterEvalQ(cl, library(spmodel))
  a <- serial$args
  a$data_object$parallel <- TRUE
  a$data_object$cl <- cl
  checked <- v0.15.0_check_solver(a)
  args$local$parallel <- TRUE
  args$local$ncores <- 2L
  fit <- do.call(spglm, args)
  expect_lte(v0.15.0_scaled_error(coef(fit), coef(serial$fit)), 1e-10)
  expect_lte(v0.15.0_scaled_error(vcov(fit), vcov(serial$fit)), 1e-10)
})

test_that('block factors preserve singleton and indefinite beta curvature', {
  set.seed(2)
  candidates <- expand.grid(w = seq(-8, 8, length.out = 65), y = c(.01, .1, .5, .9, .99))
  D <- diag(get_D('beta', candidates$w, candidates$y, NULL, 4))
  k <- which.max(D)
  expect_gt(D[k], 0)
  S <- list(matrix(D[k] / 2, 1, 1), matrix(c(2, .2, .2, 2), 2))
  dv <- c(D[k], -1, -2)
  idx <- list(1L, 2:3)
  F <- factor_laplace_blocks(S, dv, idx, FALSE, NULL)
  expect_identical(vapply(F, `[[`, '', 'type'), c('lu', 'chol'))
  A <- Matrix::bdiag(Map(function(s, i) v0.15.0_reference$get_DSigInv(Matrix::Diagonal(x = dv[i]), s), S, idx))
  V <- cbind(c(.5, -1, 2), 1)
  expect_lte(v0.15.0_scaled_error(apply_Ainv(F, V, idx), solve(A, V)), 1e-10)
  expect_lte(v0.15.0_scaled_error(apply_Ainv(F, V[, 1], idx), solve(A, V[, 1])), 1e-10)
  expect_lte(v0.15.0_scaled_error(logdet_A(F), as.numeric(determinant(A)$modulus)), 1e-8)
  U <- matrix(c(.1, .2, -.1), 3, 1)
  AU <- apply_Ainv(F, U, idx)
  mid <- matrix(1, 1, 1) + crossprod(U, AU)
  M <- structure(list(F = F, idx = idx, AU = AU, mid_inv = solve(mid)), class = 'laplace_mhinv')
  ref <- -as.matrix(v0.15.0_reference$smw_HInv(solve(A), U, matrix(1, 1, 1)))
  expect_lte(v0.15.0_scaled_error(mhinv_apply(M, V), ref %*% V), 1e-10)
  expect_lte(v0.15.0_scaled_error(mhinv_apply(M, V[, 1]), ref %*% V[, 1]), 1e-10)
  expect_lte(v0.15.0_scaled_error(mhinv_diag(M), diag(ref)), 1e-10)
  bad <- factor_laplace_blocks(list(matrix(1, 1, 1)), 1, list(1L), FALSE, NULL)
  expect_error(apply_Ainv(bad, matrix(1, 1, 1), list(1L)))
  expect_error(solve(v0.15.0_reference$get_DSigInv(Matrix::Diagonal(x = 1), matrix(1, 1, 1))))
})

test_that('LU block solves preserve the complete Newton iteration', {
  args <- v0.15.0_args(v0.15.0_data(), 'beta', 'local', 'reml')
  a <- v0.15.0_capture(args)$args
  count <- 0L
  with_mocked_bindings(v0.15.0_check_solver(a), factor_laplace_block = function(block) {
    count <<- count + 1L
    m <- block$S
    diag(m) <- diag(m) - block$D
    list(type = 'lu', A = -m)
  })
  expect_gt(count, 0L)
})

test_that('larger Bernoulli blocks preserve the solve at the nugget floor', {
  args <- v0.15.0_args(v0.15.0_data(1000L), 'bernoulli', 'local', 'reml')
  args$local <- list(index = rep(c(7, 2, 9, 4), each = 250))
  args$spcov_initial <- spcov_initial('exponential', de = .8, ie = 1e-4, range = .3, known = 'given')
  a <- v0.15.0_capture(args, solver = v0.15.0_reference$get_w_and_H_spglm)$args
  v0.15.0_check_solver(a)
})

test_that('inverse helpers preserve ordinary and Matrix inputs', {
  V <- cbind(c(.5, -1, 2), 1)
  dense <- crossprod(matrix(seq_len(9) / 10, 3)) + diag(3)
  for (M in list(dense, Matrix::Matrix(dense, sparse = FALSE), Matrix::Matrix(dense, sparse = TRUE))) {
    expect_identical(mhinv_diag(M), Matrix::diag(M))
    expect_identical(mhinv_apply(M, V), as.matrix(M %*% V))
    expect_identical(mhinv_apply(M, V[, 1]), as.matrix(M %*% V[, 1]))
    expect_identical(get_laplace_correction(M, 'poisson', rep(.2, 3), 1:3, NULL, 1),
      v0.15.0_reference$get_laplace_correction(M, 'poisson', rep(.2, 3), 1:3, NULL, 1))
  }
})

test_that('single-block beta solves preserve the general-solve fallback', {
  a <- v0.15.0_capture(v0.15.0_args(v0.15.0_data(), 'beta', 'single', 'reml'))$args
  calls <- 0L
  fallback <- get_w_and_H_spglm
  environment(fallback) <- list2env(list(chol = function(...) {
    calls <<- calls + 1L
    stop('Forced Cholesky failure')
  }), parent = environment(get_w_and_H_spglm))
  for (method in c('ml', 'reml')) {
    a$estmethod <- method
    v0.15.0_check_solver(a, fallback)
  }
  expect_gt(calls, 0L)
})
