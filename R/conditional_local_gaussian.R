# factor conditional matrix
conditional_local_factor <- function(A, label) {
  A <- as.matrix(A)
  n <- NROW(A)
  fail <- function(reason) stop(paste0(label, ": ", reason), call. = FALSE)
  if (!n || NCOL(A) != n || any(!is.finite(A))) fail("matrix must be finite and square.")
  if (any(diag(A) <= 0)) fail("matrix is not positive definite.")
  U <- tryCatch(base::chol(A), error = function(e) fail("matrix is not positive definite."))
  list(U = U)
}
