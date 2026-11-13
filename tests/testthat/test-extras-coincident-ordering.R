skip_on_cran()
skip_if_not(identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"))
skip_if_not_installed("GPvecchia")
skip_if_not_installed("spsurvey")

coincident_fixture <- function(glm = FALSE, observed = TRUE, prediction = TRUE, random = FALSE) {
  set.seed(1)
  d <- data.frame(cx = if (observed) rep(0:5, each = 2) else seq(0, 5, length.out = 12),
    cy = 0, x = rep(c(-1, 1), 6), g = factor(rep(1:3, 4)),
    y = c(2, 1, 3, 4, 2, 5, 1, 2, 4, 3, 1, 2))
  new <- d[c(1, 2, 7), ]
  new$cx <- if (prediction) c(.1, .1, 3.1) else c(.1, 1.1, 3.1)
  args <- list(formula = y ~ x, data = d, xcoord = "cx", ycoord = "cy",
    spcov_initial = spcov_initial("exponential", de = .4, ie = .2, range = .15, known = "given"))
  if (glm) args$family <- "poisson"
  if (random) {
    args$random <- ~ (1 | g)
    args$randcov_initial <- randcov_initial(g = .05, known = "given")
  }
  list(fit = do.call(if (glm) spglm else splm, args), new = new, data = d)
}


test_that("decorrelation and block prediction reject explicit coincident maxmin", {
  f <- coincident_fixture()
  params <- coef(f$fit, type = "spcov")
  expect_error(decorrelate_data(y ~ x, f$data, params, xcoord = cx, ycoord = cy,
    ordering = "maxmin"), "coincident")
  expect_error(decorrelate(y ~ x, f$data, spcov_params = params, xcoord = cx, ycoord = cy,
    ordering = "maxmin"), "coincident")
  expect_error(predict(f$fit, f$new, block = TRUE,
    local = list(method_new = "subset", size_new = 2, ordering = "maxmin")), "coincident")
  set.seed(1)
  default <- decorrelate_data(y ~ x, f$data, params, xcoord = cx, ycoord = cy)
  set.seed(1)
  explicit <- decorrelate_data(y ~ x, f$data, params, xcoord = cx, ycoord = cy, ordering = "grts")
  expect_identical(default, explicit)
})
