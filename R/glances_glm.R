#' @rdname glances.spmodel
#' @method glances spglm
#' @export
glances.spglm <- function(object, ..., sort_by = "AICc", decreasing = FALSE) {
  model_list <- c(list(object), list(...))
  if (any(!(vapply(model_list, function(x) class(x), character(1)) %in% c("spglm", "spgautor")))) {
    stop("All models must be of class spglm or spgautor", call. = FALSE)
  }
  model_list_names <- c(as.character(as.list(substitute(list(object)))[-1]), as.character(as.list(substitute(list(...)))[-1]))
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

#' @rdname glances.spmodel
#' @method glances spgautor
#' @export
glances.spgautor <- glance.spglm
