#' @export
# use generics to export tidy
generics::tidy

#' @export
# use generics to export glance
generics::glance

#' @export
# use generics to export augment
generics::augment

# CRAN release questions
release_questions <- function() {
  c(
    "Have you turned off local tests in test-loocv, test-spautor, and test-splm?",
    "Have you changed version numbers in DESCRIPTION, CITATION, and README?"
  )
}
