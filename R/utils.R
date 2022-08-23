#' @export
# use generics to export tidy
generics::tidy

#' @export
# use generics to export glance
generics::glance

#' @export
# use generics to export augment
generics::augment

# logit function
logit <- function(x) {
  if (x < 0 | x > 1) {
    stop("logit argument must be between zero and one", call. = FALSE)
  }
  log(x / (1 - x))
}

# expit function
expit <- function(x) {
  1 / (1 + exp(-x))
}

# CRAN release questions
release_questions <- function() {
  c(
    "Have you turned off local tests in test-AIC, test-anova, test-augment, test-coef, test-fitted, test-generics, test-glances, test-loocv, test-predict, test-spautor, test-splm, and test-tidy?",
    "Have you changed version numbers in DESCRIPTION, CITATION, and README?"
  )
}
