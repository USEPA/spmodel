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

test_that("latent output preserves offsets and supports calls without newdata", {
  z <- revised_vecchia_fixture("binomial"); fit <- z$fit
  set.seed(1); out <- conditional(fit,z$newdata,output=c("latent","beta"),samples=5,type="new")
  set.seed(1); joint <- draw_conditional_glm_joint(get_conditional_glm_joint(fit),5,residual=TRUE)
  target <- sweep(joint$residual+model.matrix(fit)%*%joint$beta,1,model.offset(model.frame(fit)),"+")
  expect_equal(out$latent,target,ignore_attr=TRUE)
  expect_equal(dim(conditional(fit,output="latent",samples=2)),c(18,2))
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
