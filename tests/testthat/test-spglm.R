# SPMODEL PACKAGE NEEDS TO BE INSTALLED VIA DEVTOOLS::INSTALL() BEFORE RUNNING TESTS IF THOSE TESTS HAVE PARALLELIZATION

test_that("generics work spglm point data", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))

  spmod1 <- spglm(abs(y) ~ x, "Gamma", exdata, spcov_type = "exponential", xcoord = "xcoord", ycoord = "ycoord", estmethod = "reml")
  spmod2 <- spglm(abs(y) ~ x, Gamma, exdata, spcov_type = "none", xcoord = xcoord, ycoord = ycoord, estmethod = "reml")

  # AIC
  expect_vector(AIC(spmod1))
  expect_s3_class(AIC(spmod1, spmod2), "data.frame") # turn reml fixed effects warning off

  # anova
  expect_s3_class(anova(spmod1), "data.frame")
  expect_s3_class(anova(spmod1), "anova.spglm")
  expect_s3_class(tidy(anova(spmod1)), "data.frame")
  expect_s3_class(anova(spmod1, spmod2), "data.frame")
  expect_s3_class(anova(spmod1, spmod2), "anova.spglm")
  expect_s3_class(tidy(anova(spmod1, spmod2)), "data.frame")

  # augment
  expect_s3_class(augment(spmod1), "data.frame")
  expect_s3_class(augment(spmod1, newdata = newexdata), "data.frame")

  # coef
  expect_vector(coef(spmod1))
  expect_s3_class(coef(spmod1, type = "spcov"), "exponential")
  expect_null(coef(spmod1, type = "randcov"))
  expect_vector(coefficients(spmod1))
  expect_s3_class(coefficients(spmod1, type = "spcov"), "exponential")
  expect_null(coefficients(spmod1, type = "randcov"))

  # confint
  expect_true(inherits(confint(spmod1), "matrix"))
  expect_true(inherits(confint(spmod1, parm = c("x"), level = 0.9), "matrix"))

  # cooks.distance
  expect_vector(cooks.distance(spmod1))

  # covmatrix
  expect_true(inherits(covmatrix(spmod1), "matrix"))
  expect_true(inherits(covmatrix(spmod1, newexdata), "matrix"))

  # deviance
  expect_vector(deviance(spmod1))

  # esv
  expect_s3_class(esv(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord), "data.frame")

  # fitted
  expect_vector(fitted(spmod1))
  expect_vector(fitted(spmod1, type = "link"))
  expect_type(fitted(spmod1, type = "spcov"), "list")
  expect_null(fitted(spmod1, type = "randcov"))
  expect_vector(fitted.values(spmod1))
  expect_vector(fitted.values(spmod1, type = "link"))
  expect_type(fitted.values(spmod1, type = "spcov"), "list")
  expect_null(fitted.values(spmod1, type = "randcov"))

  # formula
  expect_type(formula(spmod1), "language")

  # getCall
  expect_type(getCall(spmod1), "language")

  # glance
  expect_s3_class(glance(spmod1), "data.frame")

  # glances
  expect_s3_class(glances(spmod1), "data.frame")
  expect_s3_class(glances(spmod1, spmod2), "data.frame")

  # hatvalues
  expect_vector(hatvalues(spmod1))

  # influence
  expect_s3_class(influence(spmod1), "data.frame")

  # labels
  expect_type(labels(spmod1), "character")

  # logLik
  expect_vector(logLik(spmod1))

  # loocv
  expect_vector(loocv(spmod1))
  expect_type(loocv(spmod1, cv_predict = TRUE, se.fit = TRUE, local = TRUE), "list")

  # model.frame
  expect_s3_class(model.frame(spmod1), "data.frame")

  # model.matrix
  expect_true(inherits(model.matrix(spmod1), "matrix"))

  # model.offset
  expect_null(model.offset(model.frame(spmod1)))

  # model.response
  expect_vector(model.response(model.frame(spmod1)))

  # plot
  expect_invisible(plot(spmod1, which = 1))
  expect_invisible(plot(spmod1, which = 2))
  expect_invisible(plot(spmod1, which = 7))

  # predict
  expect_vector(predict(spmod1, newdata = newexdata))
  expect_type(predict(spmod1, newdata = newexdata, interval = "prediction", se.fit = TRUE, local = TRUE), "list")
  expect_true(inherits(predict(spmod1, newdata = newexdata, interval = "confidence", level = 0.9), "matrix"))
  expect_vector(predict(spmod1, newdata = newexdata, type = "response"))
  expect_type(predict(spmod1, newdata = newexdata, type = "response", interval = "prediction", se.fit = TRUE, local = TRUE, var_correct = FALSE), "list")
  expect_true(inherits(predict(spmod1, newdata = newexdata, type = "response", interval = "confidence", level = 0.9), "matrix"))

  # print
  expect_output(print(spmod1))
  expect_output(print(summary(spmod1)))
  expect_output(print(anova(spmod1)))

  # pseudoR2
  expect_vector(pseudoR2(spmod1))

  # residuals
  expect_vector(residuals(spmod1))
  expect_vector(residuals(spmod1, type = "response"))
  expect_vector(residuals(spmod1, type = "pearson"))
  expect_vector(residuals(spmod1, type = "standardized"))
  expect_vector(resid(spmod1))
  expect_vector(resid(spmod1, type = "response"))
  expect_vector(resid(spmod1, type = "pearson"))
  expect_vector(resid(spmod1, type = "standardized"))
  expect_vector(rstandard(spmod1))

  # summary
  expect_type(summary(spmod1), "list")

  # terms
  expect_type(terms(spmod1), "language")

  # tidy
  expect_s3_class(tidy(spmod1), "data.frame")

  # update
  expect_s3_class(update(spmod2), "spglm")

  # varcomp
  expect_s3_class(varcomp(spmod1), "data.frame")

  # vcov
  expect_true(inherits(vcov(spmod1), "matrix"))
  expect_true(inherits(vcov(spmod1, var_correct = FALSE), "matrix"))
})

test_that("generics work spglm point data with missing", {
  load(file = system.file("extdata", "exdata_M.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))

  spmod1 <- spglm(abs(y) ~ x, family = "Gamma", data = exdata_M, spcov_type = "exponential", xcoord = "xcoord", ycoord = "ycoord", estmethod = "reml")
  spmod2 <- spglm(abs(y) ~ x, family = Gamma, data = exdata_M, spcov_type = "none", xcoord = xcoord, ycoord = ycoord, estmethod = "reml")

  # AIC
  expect_vector(AIC(spmod1))
  expect_s3_class(AIC(spmod1, spmod2), "data.frame") # turn reml fixed effects warning off

  # anova
  expect_s3_class(anova(spmod1), "data.frame")
  expect_s3_class(anova(spmod1), "anova.spglm")
  expect_s3_class(tidy(anova(spmod1)), "data.frame")
  expect_s3_class(anova(spmod1, spmod2), "data.frame")
  expect_s3_class(anova(spmod1, spmod2), "anova.spglm")
  expect_s3_class(tidy(anova(spmod1, spmod2)), "data.frame")

  # augment
  expect_s3_class(augment(spmod1), "data.frame")
  expect_s3_class(augment(spmod1, newdata = spmod1$newdata), "data.frame")

  # coef
  expect_vector(coef(spmod1))
  expect_s3_class(coef(spmod1, type = "spcov"), "exponential")
  expect_null(coef(spmod1, type = "randcov"))
  expect_vector(coefficients(spmod1))
  expect_s3_class(coefficients(spmod1, type = "spcov"), "exponential")
  expect_null(coefficients(spmod1, type = "randcov"))

  # confint
  expect_true(inherits(confint(spmod1), "matrix"))
  expect_true(inherits(confint(spmod1, parm = c("x"), level = 0.9), "matrix"))

  # cooks.distance
  expect_vector(cooks.distance(spmod1))

  # covmatrix
  expect_true(inherits(covmatrix(spmod1), "matrix"))
  expect_true(inherits(covmatrix(spmod1, newdata = spmod1$newdata), "matrix"))

  # deviance
  expect_vector(deviance(spmod1))

  # esv
  expect_s3_class(esv(y ~ x, exdata_M, xcoord = xcoord, ycoord = ycoord), "data.frame")

  # fitted
  expect_vector(fitted(spmod1))
  expect_vector(fitted(spmod1, type = "link"))
  expect_type(fitted(spmod1, type = "spcov"), "list")
  expect_null(fitted(spmod1, type = "randcov"))
  expect_vector(fitted.values(spmod1))
  expect_vector(fitted.values(spmod1, type = "link"))
  expect_type(fitted.values(spmod1, type = "spcov"), "list")
  expect_null(fitted.values(spmod1, type = "randcov"))

  # formula
  expect_type(formula(spmod1), "language")

  # getCall
  expect_type(getCall(spmod1), "language")

  # glance
  expect_s3_class(glance(spmod1), "data.frame")

  # glances
  expect_s3_class(glances(spmod1), "data.frame")
  expect_s3_class(glances(spmod1, spmod2), "data.frame")

  # hatvalues
  expect_vector(hatvalues(spmod1))

  # influence
  expect_s3_class(influence(spmod1), "data.frame")

  # labels
  expect_type(labels(spmod1), "character")

  # logLik
  expect_vector(logLik(spmod1))

  # loocv
  expect_vector(loocv(spmod1))
  expect_type(loocv(spmod1, cv_predict = TRUE, se.fit = TRUE, local = TRUE), "list")

  # model.frame
  expect_s3_class(model.frame(spmod1), "data.frame")

  # model.matrix
  expect_true(inherits(model.matrix(spmod1), "matrix"))

  # model.offset
  expect_null(model.offset(model.frame(spmod1)))

  # model.response
  expect_vector(model.response(model.frame(spmod1)))

  # plot
  expect_invisible(plot(spmod1, which = 1))
  expect_invisible(plot(spmod1, which = 2))
  expect_invisible(plot(spmod1, which = 7))

  # predict
  expect_vector(predict(spmod1, newdata = newexdata))
  expect_type(predict(spmod1, newdata = newexdata, interval = "prediction", se.fit = TRUE, local = TRUE), "list")
  expect_true(inherits(predict(spmod1, newdata = newexdata, interval = "confidence", level = 0.9), "matrix"))
  expect_vector(predict(spmod1, newdata = newexdata, type = "response"))
  expect_type(predict(spmod1, newdata = newexdata, type = "response", interval = "prediction", se.fit = TRUE, local = TRUE, var_correct = FALSE), "list")
  expect_true(inherits(predict(spmod1, newdata = newexdata, type = "response", interval = "confidence", level = 0.9), "matrix"))

  # print
  expect_output(print(spmod1))
  expect_output(print(summary(spmod1)))
  expect_output(print(anova(spmod1)))

  # pseudoR2
  expect_vector(pseudoR2(spmod1))

  # residuals
  expect_vector(residuals(spmod1))
  expect_vector(residuals(spmod1, type = "response"))
  expect_vector(residuals(spmod1, type = "pearson"))
  expect_vector(residuals(spmod1, type = "standardized"))
  expect_vector(resid(spmod1))
  expect_vector(resid(spmod1, type = "response"))
  expect_vector(resid(spmod1, type = "pearson"))
  expect_vector(resid(spmod1, type = "standardized"))
  expect_vector(rstandard(spmod1))

  # summary
  expect_type(summary(spmod1), "list")

  # terms
  expect_type(terms(spmod1), "language")

  # tidy
  expect_s3_class(tidy(spmod1), "data.frame")

  # update
  expect_s3_class(update(spmod2), "spglm")

  # varcomp
  expect_s3_class(varcomp(spmod1), "data.frame")

  # vcov
  expect_true(inherits(vcov(spmod1), "matrix"))
  expect_true(inherits(vcov(spmod1, var_correct = FALSE), "matrix"))
})

test_that("generics work spglm polygon data with missing", {
  load(file = system.file("extdata", "exdata_Mpoly.rda", package = "spmodel"))

  spmod1 <- spglm(abs(y) ~ x, "Gamma", exdata_Mpoly, spcov_type = "exponential", xcoord = "xcoord", ycoord = "ycoord", estmethod = "reml")
  spmod2 <- spglm(abs(y) ~ x, Gamma, exdata_Mpoly, spcov_type = "none", xcoord = xcoord, ycoord = ycoord, estmethod = "reml")

  # AIC
  expect_vector(AIC(spmod1))
  expect_s3_class(AIC(spmod1, spmod2), "data.frame") # turn reml fixed effects warning off

  # anova
  expect_s3_class(anova(spmod1), "data.frame")
  expect_s3_class(anova(spmod1), "anova.spglm")
  expect_s3_class(tidy(anova(spmod1)), "data.frame")
  expect_s3_class(anova(spmod1, spmod2), "data.frame")
  expect_s3_class(anova(spmod1, spmod2), "anova.spglm")
  expect_s3_class(tidy(anova(spmod1, spmod2)), "data.frame")

  # augment
  expect_s3_class(augment(spmod1), "data.frame")
  expect_s3_class(augment(spmod1, newdata = spmod1$newdata), "data.frame")

  # coef
  expect_vector(coef(spmod1))
  expect_s3_class(coef(spmod1, type = "spcov"), "exponential")
  expect_null(coef(spmod1, type = "randcov"))
  expect_vector(coefficients(spmod1))
  expect_s3_class(coefficients(spmod1, type = "spcov"), "exponential")
  expect_null(coefficients(spmod1, type = "randcov"))

  # confint
  expect_true(inherits(confint(spmod1), "matrix"))
  expect_true(inherits(confint(spmod1, parm = c("x"), level = 0.9), "matrix"))

  # cooks.distance
  expect_vector(cooks.distance(spmod1))

  # covmatrix
  expect_true(inherits(covmatrix(spmod1), "matrix"))
  expect_true(inherits(covmatrix(spmod1, newdata = spmod1$newdata), "matrix"))

  # deviance
  expect_vector(deviance(spmod1))

  # esv
  expect_s3_class(esv(y ~ x, exdata_Mpoly, xcoord = xcoord, ycoord = ycoord), "data.frame")

  # fitted
  expect_vector(fitted(spmod1))
  expect_vector(fitted(spmod1, type = "link"))
  expect_type(fitted(spmod1, type = "spcov"), "list")
  expect_null(fitted(spmod1, type = "randcov"))
  expect_vector(fitted.values(spmod1))
  expect_vector(fitted.values(spmod1, type = "link"))
  expect_type(fitted.values(spmod1, type = "spcov"), "list")
  expect_null(fitted.values(spmod1, type = "randcov"))

  # formula
  expect_type(formula(spmod1), "language")

  # getCall
  expect_type(getCall(spmod1), "language")

  # glance
  expect_s3_class(glance(spmod1), "data.frame")

  # glances
  expect_s3_class(glances(spmod1), "data.frame")
  expect_s3_class(glances(spmod1, spmod2), "data.frame")

  # hatvalues
  expect_vector(hatvalues(spmod1))

  # influence
  expect_s3_class(influence(spmod1), "data.frame")

  # labels
  expect_type(labels(spmod1), "character")

  # logLik
  expect_vector(logLik(spmod1))

  # loocv
  expect_vector(loocv(spmod1))
  expect_type(loocv(spmod1, cv_predict = TRUE, se.fit = TRUE, local = TRUE), "list")

  # model.frame
  expect_s3_class(model.frame(spmod1), "data.frame")

  # model.matrix
  expect_true(inherits(model.matrix(spmod1), "matrix"))

  # model.offset
  expect_null(model.offset(model.frame(spmod1)))

  # model.response
  expect_vector(model.response(model.frame(spmod1)))

  # plot
  expect_invisible(plot(spmod1, which = 1))
  expect_invisible(plot(spmod1, which = 2))
  expect_invisible(plot(spmod1, which = 7))

  # predict
  expect_vector(predict(spmod1))
  expect_type(predict(spmod1, interval = "prediction", se.fit = TRUE, local = TRUE), "list")
  expect_true(inherits(predict(spmod1, interval = "confidence", level = 0.9), "matrix"))
  expect_vector(predict(spmod1, type = "response"))
  expect_type(predict(spmod1, type = "response", interval = "prediction", se.fit = TRUE, local = TRUE, var_correct = FALSE), "list")
  expect_true(inherits(predict(spmod1, type = "response", interval = "confidence", level = 0.9), "matrix"))

  # print
  expect_output(print(spmod1))
  expect_output(print(summary(spmod1)))
  expect_output(print(anova(spmod1)))

  # pseudoR2
  expect_vector(pseudoR2(spmod1))

  # residuals
  expect_vector(residuals(spmod1))
  expect_vector(residuals(spmod1, type = "response"))
  expect_vector(residuals(spmod1, type = "pearson"))
  expect_vector(residuals(spmod1, type = "standardized"))
  expect_vector(resid(spmod1))
  expect_vector(resid(spmod1, type = "response"))
  expect_vector(resid(spmod1, type = "pearson"))
  expect_vector(resid(spmod1, type = "standardized"))
  expect_vector(rstandard(spmod1))

  # summary
  expect_type(summary(spmod1), "list")

  # terms
  expect_type(terms(spmod1), "language")

  # tidy
  expect_s3_class(tidy(spmod1), "data.frame")

  # update
  expect_s3_class(update(spmod2), "spglm")

  # varcomp
  expect_s3_class(varcomp(spmod1), "data.frame")

  # vcov
  expect_true(inherits(vcov(spmod1), "matrix"))
  expect_true(inherits(vcov(spmod1, var_correct = FALSE), "matrix"))
})
