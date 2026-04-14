#' @title Internal Utility Function for Parameter Selection
#' @description Resolves and validates the \code{bin} and \code{sp_order} parameters.
#' If \code{bin} is not provided or valid, both \code{bin} and \code{sp_order} are generated automatically.
#' If \code{bin} is valid but \code{sp_order} is not provided or valid, only \code{sp_order} is generated.
#' This function is intended for internal use only and should not be accessed by package users.
#' @keywords internal
#'
#' @param n Integer. Length of the input vector.
#' @param bin Number or \code{NULL}. Number of bins.
#' @param sp_order Number or \code{NULL}. Spline order.
#'
#' @return A list containing the validated or automatically generated values
#' \code{bin} and \code{sp_order}.
resolve_params <- function(n, bin = NULL, sp_order = NULL) {
  if (!is.numeric(n) || length(n) != 1 || is.na(n) || n <= 0) {
    stop("Parameter 'n' must be a positive integer.")
  }

  n <- as.integer(n)

  #bin <- as.integer(bin)
  #sp_order <- as.integer(sp_order)

  if (length(bin)!=1 || is.null(bin) || bin <= 0 || is.na(bin)) {
    bin <- max(4L, round(n/10))
    bin <- min(bin, 12)
    sp_order <- round(bin/6+bin/8)
  } else {
    bin <- as.integer(bin)

    if (length(sp_order)!=1 || is.null(sp_order) || sp_order <= 0 || is.na(sp_order) ) {
      sp_order <- round(bin/6+bin/8)
    } else {
      sp_order <- as.integer(sp_order)
    }
  }

  if (sp_order >= bin) {
    stop("Spline order 'sp_order' must be less than 'bin'.")
  }

  list(bin = bin, sp_order = sp_order)
}
