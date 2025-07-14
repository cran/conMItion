#' @title Calculate Univariate Entropy
#' @description This function calculates the univariate entropy of a numeric vector using the specified
#' number of bins and spline order.
#'
#' @param x_1 A numeric vector for the only variable.
#' @param bin An integer specifying the number of bins. Default is 6.
#' @param sp_order An integer specifying the spline order. Must be less than `bin`. Default is 2.
#' @return A numeric value representing the entropy of the vector.
#'
#' @examples
#' x_1 <- rnorm(100)
#' getEntropy(x_1)
#'
#' @export
getEntropy <- function(x_1, bin = 6, sp_order = 2) {
  # Check inputs
  bin <- as.integer(bin)
  sp_order <- as.integer(sp_order)

  if (!is.numeric(x_1) || length(x_1) == 0) {
    stop("Input x_1 must be a non-empty numeric vector.")
  }

  if (!is.integer(bin) || bin <= 0) {
    stop("Parameter 'bin' must be a positive integer.")
  }

  if (!is.integer(sp_order) || sp_order <= 0) {
    stop("Parameter 'sp_order' must be a positive integer.")
  }

  if (sp_order >= bin) {
    stop("Spline order 'sp_order' must be less than 'bin'")
  }

  # Get the length of x_1
  n <- length(x_1)

  # Call the C function for entropy calculation
  out <- .C("compute_entropy_univar",
            x_1 = as.double(x_1),
            n = as.integer(n),
            bin = bin,
            sp_order = sp_order,
            eOut = as.double(0))

  # Extract the entropy value
  e <- out$eOut
  return(e)
}

#' @title Calculate Joint Entropy for Two Variables
#' @description Computes the joint entropy of two numeric vectors using the specified number of bins and spline order.
#'
#' @param x_1 A numeric vector for the first variable.
#' @param x_2 A numeric vector for the second variable. Must be the same length as `x_1`.
#' @param bin An integer specifying the number of bins. Default is 6.
#' @param sp_order An integer specifying the spline order. Must be less than `bin`. Default is 2.
#' @return A numeric value representing the joint entropy of the two vectors.
#'
#' @examples
#' x_1 <- rnorm(100)
#' x_2 <- rnorm(100)
#' getEntropyBi(x_1, x_2)
#'
#' @export
getEntropyBi <- function(x_1, x_2, bin = 6, sp_order = 2) {
  bin <- as.integer(bin)
  sp_order <- as.integer(sp_order)
  # Validate inputs are numeric vectors
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

  # Define the length of the vectors
  n <- length(x_1)

  # Call the C function for joint entropy calculation
  out <- .C("compute_entropy_bivar",
            x_1 = as.double(x_1),
            x_2 = as.double(x_2),
            n = as.integer(n),
            bin = bin,
            sp_order = sp_order,
            eOut = as.double(0))

  # Extract the joint entropy value
  joint_entropy <- out$eOut
  return(joint_entropy)
}


#' @title Calculate Joint Entropy for Three Variables
#' @description Computes the joint entropy of three numeric vectors using the specified number of bins and spline order.
#'
#' @param x_1 A numeric vector for the first variable.
#' @param x_2 A numeric vector for the second variable. Must be the same length as `x_1`.
#' @param x_3 A numeric vector for the third variable. Must be the same length as `x_1`.
#' @param bin An integer specifying the number of bins. Default is 6.
#' @param sp_order An integer specifying the spline order. Must be less than `bin`. Default is 2.
#' @return A numeric value representing the joint entropy of the three vectors.
#'
#' @examples
#' x_1 <- rnorm(100)
#' x_2 <- rnorm(100)
#' x_3 <- rnorm(100)
#' getEntropyTri(x_1, x_2, x_3)
#'
#' @export
getEntropyTri <- function(x_1, x_2, x_3, bin = 6, sp_order = 2) {
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

  # Call the C functions to compute the entropy
  out <- .C("compute_entropy_trivar",
            x_1 = as.double(x_1),
            x_2 = as.double(x_2),
            x_3 = as.double(x_3),
            n = as.integer(n),
            bin = bin,
            sp_order = sp_order,
            eOut = as.double(0))

  # Extract the joint entropy value
  joint_entropy <- out$eOut

  return(joint_entropy)
}


#' @title Calculate Joint Entropy for Four Variables
#' @description Computes the joint entropy of four numeric vectors using the specified number of bins and spline order.
#'
#' @param x_1 A numeric vector for the first variable.
#' @param x_2 A numeric vector for the second variable. Must be the same length as `x_1`.
#' @param x_3 A numeric vector for the third variable. Must be the same length as `x_1`.
#' @param x_4 A numeric vector for the fourth variable. Must be the same length as `x_1`.
#' @param bin An integer specifying the number of bins. Default is 6.
#' @param sp_order An integer specifying the spline order. Must be less than `bin`. Default is 2.
#' @return A numeric value representing the joint entropy of the four vectors.
#'
#' @examples
#' x_1 <- rnorm(100)
#' x_2 <- rnorm(100)
#' x_3 <- rnorm(100)
#' x_4 <- rnorm(100)
#' getEntropyQuadri(x_1, x_2, x_3, x_4)
#'
#' @export
getEntropyQuadri <- function(x_1, x_2, x_3, x_4, bin = 6, sp_order = 2) {
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

  # Check that x_1, x_2, x_3, x_4 are of the same length
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

  # Call the C functions to compute the entropy
  out <- .C("compute_entropy_quadrivar",
            x_1 = as.double(x_1),
            x_2 = as.double(x_2),
            x_3 = as.double(x_3),
            x_4 = as.double(x_4),
            n = as.integer(n),
            bin = bin,
            sp_order = sp_order,
            eOut = as.double(0))

  # Extract the joint entropy value
  joint_entropy <- out$eOut

  return(joint_entropy)
}
