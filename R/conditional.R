conditional <- function(object, ...) {
  UseMethod("conditional", object)
}

conditional.splm <- function(object, newdata, samples = 1, local) {

  if (missing(local)) {
    local <- NULL
  }

  local_list <- get_local_list_conditional(local, object, newdata)

}
