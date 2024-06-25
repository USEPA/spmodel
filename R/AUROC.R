#' Area Under Receiver Operating Characteristic Curve
#'
#' @description Compare area under the receiver operating characteristic curve (AUROC) for binary (e.g.,
#'   logistic) models. The area under the ROC curve provides a measure of the
#'   model's classification accuracy averaged over all possible threshold values.
#'
#' @param object A fitted model object from [spglm()] or [spgautor()]) where \code{family = "binomial"}
#'   and the response values are binary, representing a single success or failure for each datum.
#' @param ... Additional arguments to [pROC::auc()].
#'
#' @return The area under the receiver operating characteristic curve.
#'
#' @order 1
#' @export
#'
#' @examples
#' spgmod <- spglm(presence ~ elev,
#'   family = "binomial", data = moose,
#'   spcov_type = "exponential"
#' )
#' AUROC(spgmod)
#' @references
#' Robin, X., Turck, N., Hainard, A., Tiberti, N., Lisacek, F., Sanchez, J. C., & Müller, M.
#' (2011). pROC: an open-source package for R and S+ to analyze and compare ROC curves.
#' \emph{BMC bioinformatics}, 12, 1-8.
AUROC <- function(object, ...) {
  UseMethod("AUROC", object)
}

#' @rdname AUROC
#' @method AUROC spglm
#' @order 2
#' @export
AUROC.spglm <- function(object, ...) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("Install the pROC package before using AUROC().", call. = FALSE)
  } else {

    if (object$family != "binomial") {
      stop("AUROC() only available when family is \"binomial\".", call. = FALSE)
    }

    if (any(object$size != 1)) {
      stop("AUROC() only available for binary models (i.e., models whose response indicates a single success or failure).", call. = FALSE)
    }
    dotlist <- list(...)
    if (!("quiet" %in% names(dotlist))) {
      dotlist$quiet <- TRUE
    }
    as.numeric(do.call(pROC::auc, c(list(response = as.vector(object$y), predictor = fitted(object)), dotlist)))
  }
}

#' @rdname AUROC
#' @method AUROC spgautor
#' @order 3
#' @export
AUROC.spgautor <- AUROC.spglm
