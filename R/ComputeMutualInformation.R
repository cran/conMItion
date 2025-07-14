#' @title Calculate Mutual Information Between Two Vectors
#' @description Computes the mutual information (MI) between two numeric vectors
#' using the specified number of bins and spline order.
#'
#' @param x_1 A numeric vector representing the first variable.
#' @param x_2 A numeric vector representing the second variable. Must be the same length as `x_1`.
#' @param bin An integer specifying the number of bins. Default is 6.
#' @param sp_order An integer specifying the spline order. Must be less than `bin`. Default is 2.
#' @return A numeric value representing the mutual information (MI).
#'
#' @examples
#' x_1 <- rnorm(100)
#' x_2 <- rnorm(100)
#' getMI(x_1, x_2)
#'
#' @export
getMI <- function(x_1, x_2, bin = 6, sp_order = 2) {
  # Validate inputs as non-empty numeric vectors
  bin <- as.integer(bin)
  sp_order <- as.integer(sp_order)

  if (!is.numeric(x_1) || length(x_1) == 0) {
    stop("Input x_1 must be a non-empty numeric vector.")
  }

  if (!is.numeric(x_2) || length(x_2) == 0) {
    stop("Input x_2 must be a non-empty numeric vector.")
  }

  # Check that x_1 and x_2 are of the same length
  if (length(x_1) != length(x_2)) {
    stop("Vectors x_1 and x_2 must be of the same length.")
  }

  # Validate bin value
  if (!is.integer(bin) || bin <= 0) {
    stop("Parameter 'bin' must be a positive integer.")
  }

  # Validate sp_order value
  if (!is.integer(sp_order) || sp_order <= 0) {
    stop("Parameter 'sp_order' must be a positive integer.")
  }

  # Ensure spline order is less than number of bins
  if (sp_order >= bin) {
    stop("Spline order 'sp_order' must be less than 'bin'.")
  }

  # Get the length of the vectors
  n <- length(x_1)

  # Call the C functions to compute the mutual information
  out <- .C("compute_mutual_information",
            x_1 = as.double(x_1),
            x_2 = as.double(x_2),
            n = as.integer(n),
            bin = bin,
            sp_order = sp_order,
            eOut = as.double(0))

  # Extract the mutual information value
  mi <- out$eOut

  return(mi)
}


#' @title Calculate Joint Mutual Information
#' @description Computes the joint mutual information I(x_1, x_2; x_3)
#' using the specified number of bins and spline order.
#'
#' @param x_1 A numeric vector for the first variable.
#' @param x_2 A numeric vector for the second variable. Must match the length of `x_1`.
#' @param x_3 A numeric vector for the third variable. Must match the length of `x_1`.
#' @param bin An integer specifying the number of bins. Default is 6.
#' @param sp_order An integer specifying the spline order. Must be less than `bin`. Default is 2.
#' @return A numeric value representing joint mutual information (MI).
#'
#' @examples
#' x_1 <- rnorm(100)
#' x_2 <- rnorm(100)
#' x_3 <- rnorm(100)
#' getMIBi(x_1, x_2, x_3)
#'
#' @export
getMIBi <- function(x_1, x_2, x_3, bin = 6, sp_order = 2) {
  # Validate inputs as non-empty numeric vectors
  bin <- as.integer(bin)
  sp_order <- as.integer(sp_order)

  if (!is.numeric(x_1) || length(x_1) == 0) {
    stop("Input x_1 must be a non-empty numeric vector.")
  }

  if (!is.numeric(x_2) || length(x_2) == 0) {
    stop("Input x_2 must be a non-empty numeric vector.")
  }

  if (!is.numeric(x_3) || length(x_3) == 0) {
    stop("Input x_3 must be a non-empty numeric vector.")
  }

  # Check that x_1, x_2 and x_3 are of the same length
  if (length(x_1) != length(x_2)) {
    stop("Vectors x_1, x_2, and x_3 must be of the same length.")
  }

  if (length(x_2) != length(x_3)) {
    stop("Vectors x_1, x_2, and x_3 must be of the same length.")
  }

  # Validate bin value
  if (!is.integer(bin) || bin <= 0) {
    stop("Parameter 'bin' must be a positive integer.")
  }

  # Validate sp_order value
  if (!is.integer(sp_order) || sp_order <= 0) {
    stop("Parameter 'sp_order' must be a positive integer.")
  }

  # Ensure spline order is less than number of bins
  if (sp_order >= bin) {
    stop("Spline order 'sp_order' must be less than 'bin'.")
  }

  # Get the length of the vectors
  n <- length(x_1)

  # Call the C functions to compute the mutual information
  out <- .C("compute_mutual_information_12_to_3",
            x_1 = as.double(x_1),
            x_2 = as.double(x_2),
            x_3 = as.double(x_3),
            n = as.integer(n),
            bin = bin,
            sp_order = sp_order,
            eOut = as.double(0))

  # Extract the mutual information value
  mi <- out$eOut

  return(mi)
}


#' @title Calculate Conditional Mutual Information
#' @description Computes the conditional mutual information I(x_1; x_2 | x_3)
#' using the specified number of bins and spline order.
#'
#' @param x_1 A numeric vector for the first variable.
#' @param x_2 A numeric vector for the second variable. Must match `x_1` length.
#' @param x_3 A numeric vector for the condition variable. Must match `x_1` length.
#' @param bin An integer specifying the number of bins. Default is 6.
#' @param sp_order An integer specifying the spline order. Must be less than `bin`. Default is 2.
#' @return A numeric value representing the conditional mutual information (CMI).
#'
#' @examples
#' x_1 <- rnorm(100)
#' x_2 <- rnorm(100)
#' x_3 <- rnorm(100)
#' getCMI(x_1, x_2, x_3)
#'
#' @export
getCMI <- function(x_1, x_2, x_3, bin = 6, sp_order = 2) {
  # Validate inputs as non-empty numeric vectors
  bin <- as.integer(bin)
  sp_order <- as.integer(sp_order)

  if (!is.numeric(x_1) || length(x_1) == 0) {
    stop("Input x_1 must be a non-empty numeric vector.")
  }

  if (!is.numeric(x_2) || length(x_2) == 0) {
    stop("Input x_2 must be a non-empty numeric vector.")
  }

  if (!is.numeric(x_3) || length(x_3) == 0) {
    stop("Input x_3 must be a non-empty numeric vector.")
  }

  # Check that x_1, x_2 and x_3 are of the same length
  if (length(x_1) != length(x_2)) {
    stop("Vectors x_1, x_2, and x_3 must be of the same length.")
  }

  if (length(x_2) != length(x_3)) {
    stop("Vectors x_1, x_2, and x_3 must be of the same length.")
  }

  # Validate bin value
  if (!is.integer(bin) || bin <= 0) {
    stop("Parameter 'bin' must be a positive integer.")
  }

  # Validate sp_order value
  if (!is.integer(sp_order) || sp_order <= 0) {
    stop("Parameter 'sp_order' must be a positive integer.")
  }

  # Ensure spline order is less than number of bins
  if (sp_order >= bin) {
    stop("Spline order 'sp_order' must be less than 'bin'.")
  }

  # Get the length of the vectors
  n <- length(x_1)

  # Call the C functions to compute the conditional mutual information
  out <- .C("compute_conditional_mutual_information_12_to_3",
            x_1 = as.double(x_1),
            x_2 = as.double(x_2),
            x_3 = as.double(x_3),
            n = as.integer(n),
            bin = bin,
            sp_order = sp_order,
            eOut = as.double(0))

  # Extract the conditional mutual information value
  mi <- out$eOut

  return(mi)
}


#' @title Calculate Bivariate Conditional Mutual Information
#' @description Computes conditional mutual information I(x_1; x_2 | x_3, x_4)
#' using the specified number of bins and spline order.
#'
#' @param x_1 A numeric vector for the first variable.
#' @param x_2 A numeric vector for the second variable. Must match `x_1` length.
#' @param x_3 A numeric vector for the first condition variable. Must match `x_1` length.
#' @param x_4 A numeric vector for the second condition variable. Must match `x_1` length.
#' @param bin An integer specifying the number of bins. Default is 6.
#' @param sp_order An integer specifying the spline order. Must be less than `bin`. Default is 2.
#' @return A numeric value representing the bivariate conditional mutual information (CMI).
#'
#' @examples
#' x_1 <- rnorm(100)
#' x_2 <- rnorm(100)
#' x_3 <- rnorm(100)
#' x_4 <- rnorm(100)
#' getCMIBiCondi(x_1, x_2, x_3, x_4)
#'
#' @export
getCMIBiCondi <- function(x_1, x_2, x_3, x_4, bin = 6, sp_order = 2) {
  # Validate inputs as non-empty numeric vectors
  bin <- as.integer(bin)
  sp_order <- as.integer(sp_order)

  if (!is.numeric(x_1) || length(x_1) == 0) {
    stop("Input x_1 must be a non-empty numeric vector.")
  }

  if (!is.numeric(x_2) || length(x_2) == 0) {
    stop("Input x_2 must be a non-empty numeric vector.")
  }

  if (!is.numeric(x_3) || length(x_3) == 0) {
    stop("Input x_3 must be a non-empty numeric vector.")
  }

  if (!is.numeric(x_4) || length(x_4) == 0) {
    stop("Input x_4 must be a non-empty numeric vector.")
  }

  # Check that x_1, x_2, x_3 and x_4 are of the same length
  if (length(x_1) != length(x_2)) {
    stop("Vectors x_1, x_2, x_3 and x_4 must be of the same length.")
  }

  if (length(x_2) != length(x_3)) {
    stop("Vectors x_1, x_2, x_3 and x_4 must be of the same length.")
  }

  if (length(x_3) != length(x_4)) {
    stop("Vectors x_1, x_2, x_3 and x_4 must be of the same length.")
  }

  # Validate bin value
  if (!is.integer(bin) || bin <= 0) {
    stop("Parameter 'bin' must be a positive integer.")
  }

  # Validate sp_order value
  if (!is.integer(sp_order) || sp_order <= 0) {
    stop("Parameter 'sp_order' must be a positive integer.")
  }

  # Ensure spline order is less than number of bins
  if (sp_order >= bin) {
    stop("Spline order 'sp_order' must be less than 'bin'.")
  }

  # Get the length of the vectors
  n <- length(x_1)

  # Call the C functions to compute the conditional mutual information
  out <- .C("compute_conditional_mutual_information_12_to_34",
            x_1 = as.double(x_1),
            x_2 = as.double(x_2),
            x_3 = as.double(x_3),
            x_3 = as.double(x_4),
            n = as.integer(n),
            bin = bin,
            sp_order = sp_order,
            eOut = as.double(0))

  # Extract the conditional mutual information value
  mi <- out$eOut

  return(mi)
}
