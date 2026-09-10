#' @rdname glances
#' @method glances spglm
#' @order 6
#' @export
glances.spglm <- function(object, ..., sort_by = "AICc", decreasing = FALSE, warning = TRUE) {
  # mirrors glances.splm() but restricted to GLM-type model classes
  model_list <- c(list(object), list(...))
  if (any(!(vapply(model_list, function(x) class(x), character(1)) %in% c("spglm", "spgautor")))) {
    stop("All models must be of class spglm or spgautor", call. = FALSE)
  }
  # recover the argument names/expressions as passed by the caller so the
  # output table can label each row by model name, not just position
  model_list_names <- c(as.character(as.list(substitute(list(object)))[-1]), as.character(as.list(substitute(list(...)))[-1]))
  if (warning && length(model_list) > 1) {
    check_likstat_use(model_list)
  }
  model_glance <- lapply(model_list, function(x) glance(x))
  model_bind <- do.call(rbind, model_glance)
  model_bind <- cbind(data.frame(model = model_list_names), model_bind)
  if (sort_by == "order") {
    model_bind <- model_bind[order(seq_len(NROW(model_bind)), decreasing = decreasing), , drop = FALSE]
  } else {
    model_bind <- model_bind[order(model_bind[[substitute(sort_by)]], decreasing = decreasing), , drop = FALSE]
  }
  tibble::as_tibble(model_bind)
}

#' @rdname glances
#' @method glances spgautor
#' @order 7
#' @export
glances.spgautor <- glances.spglm

#' @rdname glances
#' @method glances spglm_list
#' @order 8
#' @export
glances.spglm_list <- glances.splm_list

#' @rdname glances
#' @method glances spgautor_list
#' @order 9
#' @export
glances.spgautor_list <- glances.spautor_list
