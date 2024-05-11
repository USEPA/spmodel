#' Area Under Receiver Operating Characteristic (ROC) Curve
#'
#' @description Compare area under the receiver operating characteristic (ROC) curve for binary (e.g.,
#'   logistic) models. The area under the ROC curve (AUC) provided a measure of the
#'   model's classification acuracy averaged over all possible threshold values.
#'
#' @param object A fitted model object from [spglm()] or [spgautor()]).
#' @param ... Additional arguments to [pROC::auc()].
#'
#' @return The area under the receiver operator curve.
#'
#' @order 1
#' @export
#'
#' @examples
#' spgmod <- spglm(presence ~ elev,
#'   family = "binomial", data = moose,
#'   spcov_type = "exponential"
#' )
#' AUC(spgmod)
#' @references
#' Robin, X., Turck, N., Hainard, A., Tiberti, N., Lisacek, F., Sanchez, J. C., & Müller, M.
#' (2011). pROC: an open-source package for R and S+ to analyze and compare ROC curves.
#' \emph{BMC bioinformatics}, 12, 1-8.
AUC <- function(object, ...) {
  UseMethod("AUC", object)
}

#' @rdname AUC
#' @method AUC spglm
#' @order 2
#' @export
AUC.spglm <- function(object, ...) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("Install the pROC package before using AUC().", call. = FALSE)
  } else {

    if (object$family != "binomial") {
      stop("AUC() only available when family is \"binomial\".", call. = FALSE)
    }

    if (any(object$size != 1)) {
      stop("AUC() only available for binary models (i.e., models whose response indicates a single success or failure).", call. = FALSE)
    }
    dotlist <- list(...)
    if (!("quiet" %in% names(dotlist))) {
      dotlist$quiet <- TRUE
    }
    as.numeric(do.call(pROC::auc, c(list(response = as.vector(object$y), predictor = fitted(object)), dotlist)))
  }
}

#' @rdname AUC
#' @method AUC spgautor
#' @order 3
#' @export
AUC.spgautor <- AUC.spglm
