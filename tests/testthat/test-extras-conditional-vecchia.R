skip_on_cran()
skip_if_not(identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"))

set.seed(1)

# Small fixed-covariance fixtures keep the extras checks focused on simulation.
revised_vecchia_fixture <- function(family = "poisson", random = FALSE, anisotropy = FALSE) {
  dat <- data.frame(cx = seq(0, 1, length.out = 18), cy = rep(c(0, .3, .1), 6),
    x = rep(c(-1, 0, 1), 6), off = seq(-.3, .2, length.out = 18),
    group = factor(rep(1:3, 6)), y = c(2,1,4,3,6,2,1,4,5,2,3,6,4,2,5,3,1,7))
  if (family == "binomial") dat$y <- dat$y %% 5
  new <- dat[c(2,6,11,15), ]; new$cx <- new$cx + .025
  initial <- if (anisotropy) spcov_initial("exponential",de=.4,ie=.2,range=.3,rotate=.5,scale=.4,known="given") else
    spcov_initial("exponential",de=.4,ie=.2,range=.3,known="given")
  fit <- do.call(spglm,list(formula=if(family=="binomial")cbind(y,5-y)~x+offset(off) else y~x+offset(off),
    family=family,data=dat,xcoord="cx",ycoord="cy",spcov_initial=initial,
    random=if(random)~(1|group) else NULL,randcov_initial=if(random)randcov_initial(group=.2,known="given") else NULL))
  list(fit=fit,newdata=new)
}

test_that("revised Gaussian preparation and conditioning match dense algebra", {
  S <- exp(-abs(outer(1:5,1:5,"-")))+diag(.2,5)
  X <- cbind(1,1:5); w <- seq(-.2,.4,length.out=5); D <- c(-1,0,.1,-.4,-2)
  V <- solve(solve(S)-diag(D)); M <- V%*%solve(S,X)
  a <- prepare_conditional_local_gaussian(S,X,w,D)
  expect_equal(a$V,V,tolerance=1e-11); expect_equal(a$M,M,tolerance=1e-11)
  for(N in list(integer(),c(2,4))) {
    z <- prepare_conditional_latent_operator(a,1,N,"test")
    g <- if(length(N))as.numeric(V[1,N,drop=FALSE]%*%solve(V[N,N])) else numeric()
    expect_equal(z$weights,g)
    expect_equal(z$intercept,w[1]-sum(g*w[N]))
    expect_equal(z$variance,V[1,1]-sum(g*V[N,1]))
    expect_equal(z$coefficient,as.numeric(M[1,]-colSums(g*M[N,,drop=FALSE])))
  }
  expect_error(prepare_conditional_local_gaussian(S,X,w,rep(100,5)),"latent precision.*positive definite")
  # The full Gaussian must be proper even when the target-only precision is positive.
  expect_error(prepare_conditional_local_gaussian(diag(2),matrix(1,2,1),c(0,0),c(-1,2)),
    "latent precision.*positive definite")
  # Follow chol's upper-triangle convention even if unused lower entries differ.
  bad <- S; bad[2,1] <- bad[2,1]+.01
  expect_equal(conditional_local_factor(bad,"test")$U,conditional_local_factor(S,"test")$U)
  rho <- 1-1e-12; near <- matrix(c(1,rho,rho,1),2)
  z <- prepare_conditional_local_weights(near,1,2,"near")
  expect_lt(abs(z$variance-(1-rho)*(1+rho)),1e-22)
  # GLM Step 2 can build the QR factor from its saved precision Cholesky instead.
  precision <- conditional_local_factor(base::chol2inv(base::chol(near)),"near precision")
  z_precision <- prepare_conditional_local_weights(near,1,2,"near",precision)
  expect_equal(z_precision$variance/((1-rho)*(1+rho)),1,tolerance=1e-3)
  scales <- c(1e-6,1,1e6,1e-2,1e2)
  scaled <- S*outer(scales,scales)
  ans <- conditional_local_solve(conditional_local_factor(scaled,"scaled"),matrix(1,5,1))
  expect_equal(as.numeric(scaled%*%ans),rep(1,5),tolerance=1e-3)
})

test_that("revised defaults, latent output, and factor preparation are coherent", {
  z <- revised_vecchia_fixture(); fit <- z$fit
  settings <- list(approximation="vecchia",size=5,ordering="none",chunk_size=2)
  expect_equal(get_local_list_conditional(list(approximation="vecchia"),fit,z$newdata)$size,60)
  lm <- fit; class(lm) <- "splm"
  expect_equal(get_local_list_conditional(list(approximation="vecchia"),lm,z$newdata)$size,30)
  expect_error(get_local_list_conditional(c(settings,list(size_response=3)),fit,z$newdata),"not supported")
  prepare <- prepare_conditional_local_gaussian; n_calls <- 0L; sizes <- integer()
  pred_prepare <- prepare_conditional_vecchia; prediction_calls <- 0L
  scores <- conditional_vecchia_scores; score_calls <- 0L
  local_mocked_bindings(get_d=function(...)stop("gradient must not be called"),
    get_conditional_glm_joint=function(...)stop("dense preparation"),
    conditional_vecchia_scores=function(...) { score_calls <<- score_calls+1L; scores(...) },
    prepare_conditional_vecchia=function(...) { prediction_calls <<- prediction_calls+1L; pred_prepare(...) },
    prepare_conditional_local_gaussian=function(Sigma,...) {
      n_calls <<- n_calls+1L; sizes <<- c(sizes,NROW(Sigma)); prepare(Sigma,...)
    })
  for(chunk in c(1,20)) {
    settings$chunk_size <- chunk; n_calls <- prediction_calls <- score_calls <- 0L
    ans <- conditional(fit,z$newdata,output=c("newdata","latent","beta"),samples=7,local=settings)
    expect_equal(n_calls,fit$n); expect_equal(prediction_calls,1L); expect_lte(max(sizes),10)
    expect_equal(score_calls,fit$n+nrow(z$newdata))
    expect_equal(dim(ans$latent),c(18,7)); expect_true(all(is.finite(ans$newdata)))
  }
  set.seed(1); all <- conditional(fit,z$newdata,output=c("newdata","beta","latent"),samples=10,local=settings)
  set.seed(1); part <- conditional(fit,z$newdata,output=c("newdata","beta"),samples=10,local=settings)
  expect_identical(all$newdata,part$newdata); expect_identical(all$beta,part$beta)
  prediction_calls <- 0L
  only <- conditional(fit,output="latent",samples=2,local=settings)
  expect_equal(prediction_calls,0L)
  expect_equal(dim(only),c(18,2))
  expect_equal(rownames(only),rownames(model.matrix(fit)))
})

test_that("Vecchia ranks full covariance and restricts partition eligibility", {
  fit <- revised_vecchia_fixture(random=TRUE)$fit
  dat <- fit$obdata
  # Row 2 is nearer, but row 4 shares the target's random-effect level.
  dat$cx <- seq_len(nrow(dat)); dat$cy <- 0
  dat$cx[1:4] <- c(0,.1,3,.11)
  ctx <- get_conditional_vecchia_covariance(fit,dat)
  expect_equal(conditional_vecchia_neighbors(ctx,1,c(2L,4L),1,"distance"),2L)
  expect_equal(conditional_vecchia_neighbors(ctx,1,c(2L,4L),1,"covariance"),4L)
  full_object <- fit; full_object$obdata <- dat
  full <- covmatrix(full_object)
  expect_equal(-conditional_vecchia_scores(ctx,1,"covariance")[-1],
    abs(as.numeric(full[1,-1])),tolerance=1e-12)

  dat$part <- factor(rep(c("a","b"),length.out=nrow(dat)))
  full_object$partition_factor <- ~part; full_object$obdata <- dat
  ctx <- get_conditional_vecchia_covariance(full_object,dat)
  for(method in c("distance","covariance")) {
    selected <- conditional_vecchia_neighbors(ctx,1,2:18,60,method)
    expect_equal(selected,seq(3L,17L,by=2L))
    expect_length(conditional_vecchia_neighbors(ctx,1,c(2L,4L),60,method),0)
  }
  expect_equal(conditional_vecchia_neighbors(ctx,1,1:18,60,"all"),2:18)
  # Both model classes use the same prediction preparation and partition rules.
  new <- dat[1,,drop=FALSE]; new$part <- factor("new",levels=c("a","b","new"))
  for(model_class in c("spglm","splm")) {
    class(full_object) <- model_class
    op <- prepare_conditional_vecchia(full_object,new,list(size=3,method="covariance",order=1L))$operators[[1]]
    expect_length(op$neighbors,0)
    expect_true(is.finite(op$variance) && op$variance>0)
  }
})

test_that("full observed neighborhoods reproduce exact fitted-centered moments", {
  z <- revised_vecchia_fixture(); fit <- z$fit
  local <- get_local_list_conditional(list(approximation="vecchia",size=100,ordering="none"),fit,z$newdata)
  context <- get_conditional_vecchia_glm_context(fit,z$newdata,model.matrix(~x,z$newdata),local)
  X <- model.matrix(fit); eta <- fitted(fit,type="link"); w <- eta-model.offset(model.frame(fit))
  V <- solve(solve(covmatrix(fit))-get_D(fit$family,eta,fit$y,fit$size,numeric()))
  V <- as.matrix(V); M <- V%*%solve(covmatrix(fit),X)
  F <- matrix(0,18,18); T <- matrix(0,18,2); a <- v <- numeric(18)
  for(op in context$operators) {
    F[op$i,op$neighbors] <- op$weights; T[op$i,] <- op$coefficient
    a[op$i] <- op$intercept; v[op$i] <- op$variance
  }
  H <- solve(diag(18)-F)
  expect_equal(as.numeric(H%*%a),as.numeric(w),tolerance=1e-10)
  expect_equal(H%*%T,unname(M),ignore_attr=TRUE,tolerance=1e-10)
  expect_equal(H%*%diag(v)%*%t(H),V,ignore_attr=TRUE,tolerance=1e-10)
  set.seed(1); exact <- conditional(fit,z$newdata,output=c("newdata","latent"),samples=4,local=FALSE)
  set.seed(1); all <- conditional(fit,z$newdata,output=c("newdata","latent"),samples=4,
    local=list(approximation="vecchia",method="all",ordering="none"))
  expect_identical(exact,all)
})

test_that("Vecchia defaults and all-neighbor exact dispatch apply to both model classes", {
  z <- revised_vecchia_fixture()
  lm <- splm(y ~ x + offset(off), z$fit$obdata, xcoord = cx, ycoord = cy,
    ddf = "asymptotic", spcov_initial = spcov_initial("exponential",
      de = .4, ie = .2, range = .3, known = "given"))
  for (fit in list(lm, z$fit)) {
    for (setting in list(TRUE, list())) {
      local <- get_local_list_conditional(setting, fit, z$newdata)
      expect_identical(local$approximation, "vecchia")
      expect_equal(local$size, if (inherits(fit, "spglm")) 60 else 30)
      expect_false(local$exact)
      set.seed(1)
      draws <- conditional(fit, z$newdata, samples = 3, local = setting)
      set.seed(1)
      explicit <- conditional(fit, z$newdata, samples = 3,
        local = list(approximation = "vecchia"))
      expect_identical(draws, explicit)
    }
    large <- fit; large$n <- 5001
    expect_message(local <- get_local_list_conditional(NULL, large, z$newdata), "local = TRUE")
    expect_identical(local$approximation, "vecchia")
    expect_false(local$exact)
    expect_error(get_local_list_conditional(list(size_base = 5), fit, z$newdata), "approximation.*low-rank")
  }
  local_mocked_bindings(
    conditional_vecchia_order = function(...) stop("exact calls must not order sites"),
    get_conditional_vecchia = function(...) stop("exact calls must not use Vecchia"),
    get_conditional_vecchia_glm = function(...) stop("exact calls must not use Vecchia"))
  for (fit in list(lm, z$fit)) {
    output <- c("all", if (inherits(fit, "spglm")) "latent")
    set.seed(1)
    exact <- conditional(fit, z$newdata, samples = 4, output = output, local = FALSE)
    for (ordering in c("none", "random", "grts")) {
      set.seed(1)
      all <- conditional(fit, z$newdata, samples = 4, output = output,
        local = list(method = "all", ordering = ordering))
      expect_identical(all, exact)
    }
    large <- fit; large$n <- 5001
    expect_warning(local <- get_local_list_conditional(list(method = "all"), large, z$newdata), NA)
    expect_true(local$exact)
  }
})

test_that("the requested ordering is applied to both observed and prediction sites", {
  z <- revised_vecchia_fixture()
  original <- conditional_vecchia_order
  calls <- list()
  local_mocked_bindings(conditional_vecchia_order = function(coords, ordering = "maxmin") {
    value <- original(coords, ordering)
    calls[[length(calls) + 1L]] <<- list(size = NROW(coords), method = ordering, order = value)
    value
  })
  for (ordering in c("none", "random", "coordinate", "middleout", "outsidein", "maxmin")) {
    calls <- list()
    set.seed(1)
    draws <- conditional(z$fit, z$newdata, samples = 3, output = c("newdata", "latent"),
      local = list(size = 5, ordering = ordering))
    expect_identical(vapply(calls, `[[`, "", "method"), rep(ordering, 2))
    expect_equal(vapply(calls, `[[`, 0L, "size"), c(4, 18))
    expect_identical(sort(calls[[2]]$order), seq_len(18))
    if (ordering == "none") expect_identical(calls[[2]]$order, seq_len(18))
    expect_identical(rownames(draws$latent), rownames(model.matrix(z$fit)))
    expect_true(all(is.finite(draws$newdata)))
  }
})

test_that("GRTS orders observed latent sites as well as prediction sites", {
  skip_if_not_installed("spsurvey")
  z <- revised_vecchia_fixture()
  original <- conditional_vecchia_order
  calls <- character()
  local_mocked_bindings(conditional_vecchia_order = function(coords, ordering = "maxmin") {
    calls <<- c(calls, ordering)
    original(coords, ordering)
  })
  set.seed(1)
  draws <- conditional(z$fit, z$newdata, samples = 3, output = c("newdata", "latent"),
    local = list(size = 5, ordering = "grts"))
  expect_identical(calls, c("grts", "grts"))
  expect_equal(dim(draws$latent), c(18, 3))
  expect_true(all(is.finite(draws$latent)))
})

test_that("small covariance blocks preserve anisotropy, random effects and diagtol", {
  for(random in c(FALSE,TRUE)) {
    z <- revised_vecchia_fixture(random=random,anisotropy=TRUE); fit <- z$fit
    fit$diagtol <- .5
    ctx <- get_conditional_vecchia_covariance(fit,fit$obdata)
    expect_equal(conditional_vecchia_covariance(ctx,c(1,4,7)),covmatrix(fit)[c(1,4,7),c(1,4,7)],ignore_attr=TRUE)
    ans <- conditional(fit,z$newdata,samples=3,local=list(approximation="vecchia",size=3,ordering="none"))
    expect_true(all(is.finite(ans)))
  }
})

test_that("latent output offsets and low-rank limitations are explicit", {
  z <- revised_vecchia_fixture("binomial"); fit <- z$fit
  set.seed(1); out <- conditional(fit,z$newdata,output=c("latent","beta"),samples=5,local=FALSE,type="new")
  set.seed(1); joint <- draw_conditional_glm_joint(get_conditional_glm_joint(fit),5,residual=TRUE)
  target <- sweep(joint$residual+model.matrix(fit)%*%joint$beta,1,model.offset(model.frame(fit)),"+")
  expect_equal(out$latent,target,ignore_attr=TRUE)
  expect_error(conditional(fit,z$newdata,output="latent",local=list(approximation = "low-rank", size_base=5,reorder_base="none")),"reduced low-rank")
  expect_equal(dim(conditional(fit,output="latent",samples=2,local=list(approximation = "low-rank", method_base="all",method_new="all"))),c(18,2))
  for(size in list(0,NA,Inf,1.5,c(1,2))) expect_error(conditional(fit,z$newdata,local=list(approximation="vecchia",size=size)),"size must")
})

test_that("coincident observed sites remain separate latent draws", {
  xy <- rbind(c(0,0),c(0,0),c(1,1),c(1,1),c(2,1))
  ord <- conditional_vecchia_order(xy)
  expect_equal(sort(ord),1:5)
  expect_equal(conditional_vecchia_order(matrix(0,3,2)),1:3)
})

test_that("areal latent output includes offsets without changing prediction draws", {
  dat <- data.frame(x=seq(-1,1,length.out=12),off=seq(-.3,.4,length.out=12),
    y=c(1,2,3,1,4,2,3,5,2,3,4,1))
  W <- 1*(abs(outer(1:12,1:12,"-"))==1)
  for(missing in list(c(3,8),integer())) {
    data <- dat; data$y[missing] <- NA
    fit <- spgautor(y~x+offset(off),data=data,W=W,family="poisson",
      spcov_initial=spcov_initial("car",de=.4,ie=.2,range=.2,known="given"))
    set.seed(1); output <- conditional(fit,output=c("latent","beta"),samples=3)
    set.seed(1); draw <- draw_conditional_glm_joint(get_conditional_glm_joint(fit),3)
    expect_equal(output$latent,sweep(draw$w,1,model.offset(model.frame(fit)),"+"),ignore_attr=TRUE)
    expect_equal(output$beta,draw$beta)
    if(length(missing)) {
      set.seed(1); all <- conditional(fit,output=c("all","latent"),samples=3)
      set.seed(1); usual <- conditional(fit,output="all",samples=3)
      expect_identical(all$newdata,usual$newdata)
      expect_equal(nrow(all$latent),10)
      expect_equal(rownames(all$latent),as.character(fit$observed_index))
    }
  }
})
