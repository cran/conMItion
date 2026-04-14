#' @title Normalized Mutual Information Between Matrix and Vector
#'
#' @description Computes the normalized mutual information (MI) between each row of a matrix and
#' a numeric vector normalized by the mutual information of the vector with itself using the specified number of bins and spline order.
#'
#' @param mat A numeric matrix. For example, each row represents a gene and each column represents a sample.
#' @param vec A numeric vector, with length equal to the number of columns in `mat`.
#' @param bin An integer specifying the number of bins. Default is NULL.
#' @param sp_order An integer specifying the spline order. Must be less than `bin`. Default is NULL.
#'
#' @return A numeric vector representing the normalized mutual information (MI) between each row of `mat` and `vec`.
#'
#' @examples
#' mat <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' vec <- rnorm(100)
#' MImat2vec(mat, vec)
#'
#' @export
MImat2vec <- function(mat, vec, bin = NULL, sp_order = NULL) {
  # Validate inputs as non-empty numeric vectors

  if (!is.numeric(vec) || length(vec) == 0) {
    stop("Input vec must be a non-empty numeric vector.")
  }

  if (!is.matrix(mat) || !is.numeric(mat) || is.na(sum(mat)) ) {
    stop("Input mat must be a non-NA numeric matrix.")
  }

  # Check that mat and vec are of the same length
  if (ncol(mat) != length(vec)) {
    stop("The column of mat and the length of vector vec must be the same.")
  }

  max_MI <- getMI(vec, vec, bin = bin, sp_order = sp_order)
  # Ensure vec contains information
  if (max_MI == 0) {
    stop("Vector 'vec' must contain information.")
  }

  mi_vector <- sapply(1:nrow(mat), function(x)
    getMI(mat[x,], vec, bin = bin, sp_order = sp_order))

  normalized_mi_vector <- mi_vector / max_MI
  names(normalized_mi_vector) <- row.names(mat)

  return(normalized_mi_vector)
}


#' @title Normalized Conditional Mutual Information Between Matrix and Vector
#' @description Computes the normalized conditional mutual information (CMI) between each row of a matrix and a vector,
#' given a third condition vector, normalized by the mutual information of the vector with itself using the specified bins and spline order.
#'
#' @param mat A numeric matrix. For example, each row represents a gene and each column represents a sample.
#' @param vec A numeric vector, with length equal to the number of columns in `mat`.
#' @param condi A numeric condition vector, matching the number of columns in `mat`.
#' @param bin An integer specifying the number of bins. Default is NULL.
#' @param sp_order An integer specifying the spline order. Must be less than `bin`. Default is NULL.
#' @return A numeric vector representing the normalized conditional mutual information (CMI) between each row of `mat` and `vec`, given `condi`.
#'
#' @examples
#' mat <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' vec <- rnorm(100)
#' condi <- rnorm(100)
#' CMImat2vec(mat, vec, condi)
#'
#' @export
CMImat2vec <- function(mat, vec, condi, bin = NULL, sp_order = NULL) {
  # Validate inputs as non-empty numeric vectors

  if (!is.numeric(vec) || length(vec) == 0) {
    stop("Input vec must be a non-empty numeric vector.")
  }

  if (!is.numeric(condi) || length(condi) == 0) {
    stop("Input condi must be a non-empty numeric vector.")
  }

  if (!is.matrix(mat) || !is.numeric(mat) || is.na(sum(mat)) ) {
    stop("Input mat must be a non-NA numeric matrix.")
  }

  # Check that mat, vec and condi are of the same length
  if (ncol(mat) != length(vec)) {
    stop("The column of mat and the length of vec must be the same.")
  }

  if (ncol(mat) != length(condi)) {
    stop("The column of mat and the length of vector condi must be the same.")
  }

  max_MI <- getMI(vec, vec, bin = bin, sp_order = sp_order)
  # Ensure vec contains information
  if (max_MI == 0) {
    stop("Vector 'vec' must contain information.")
  }

  condi_mi_vector <- sapply(1:nrow(mat), function(x)
    getCMI(mat[x,], vec, condi, bin = bin, sp_order = sp_order))

  normalized_condi_mi_vector <- condi_mi_vector / max_MI
  names(normalized_condi_mi_vector) <- row.names(mat)

  return(normalized_condi_mi_vector)
}


#' @title Normalized Mutual Information Between Two Matrices
#' @description Computes the normalized mutual information (MI) between corresponding rows of two matrices
#' normalized by their individual information content, using the specified number of bins and spline order.
#'
#' @param mat1 A numeric matrix. For example, each row represents a gene and each column represents a sample.
#' @param mat2 Another numeric matrix to compare against. Must have the same dimensions as `mat1`.
#' @param bin An integer specifying the number of bins. Default is NULL.
#' @param sp_order An integer specifying the spline order. Must be less than `bin`. Default is NULL.
#' @return A numeric vector where each element corresponds to the normalized mutual information (MI) between
#' respective rows of `mat1` and `mat2`.
#'
#' @examples
#' mat1 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' mat2 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' MImat2mat(mat1, mat2)
#'
#' @export
MImat2mat <- function(mat1, mat2, bin = NULL, sp_order = NULL) {
  # Validate inputs as non-empty numeric vectors

  if (!is.matrix(mat1) || !is.numeric(mat1) || is.na(sum(mat1)) ) {
    stop("Input mat must be a non-NA numeric matrix.")
  }

  if (!is.matrix(mat2) || !is.numeric(mat2) || is.na(sum(mat2)) ) {
    stop("Input mat must be a non-NA numeric matrix.")
  }

  # Check that 2 mats are of the same size
  if (ncol(mat1) != ncol(mat2) | nrow(mat1) != nrow(mat2) ) {
    stop("The size of mat1 and mat2 must be the same.")
  }

  max_MI <- sapply(1:nrow(mat1), function(x)
    max(getMI(mat1[x,], mat1[x,], bin = bin, sp_order = sp_order),
        getMI(mat2[x,], mat2[x,], bin = bin, sp_order = sp_order))
  )

  mi_vector <- sapply(1:nrow(mat1), function(x)
    getMI(mat1[x,], mat2[x,], bin = bin, sp_order = sp_order))

  normalized_mi_vector <- ifelse(max_MI == 0, 0, mi_vector / max_MI)
  names(normalized_mi_vector) <- row.names(mat1)
  normalized_mi_vector[normalized_mi_vector<0] = 0
  return(normalized_mi_vector)
}


#' @title Normalized Conditional Mutual Information Between Two Matrices
#' @description Computes the normalized conditional mutual information (CMI) between corresponding rows of two matrices,
#' given a condition variable, normalized by their individual information content. CMI is calculated using the specified number of bins and spline order.
#'
#' @param mat1 A numeric matrix. For example, each row represents a gene and each column represents a sample.
#' @param mat2 Another numeric matrix to compare against. Must have the same dimensions as `mat1`.
#' @param condi A numeric condition vector, matching the number of columns in `mat1`.
#' @param bin An integer specifying the number of bins. Default is NULL.
#' @param sp_order An integer specifying the spline order. Must be less than `bin`. Default is NULL.
#' @return A numeric vector representing the normalized conditional mutual information (CMI) between pairs of rows
#' from `mat1` and `mat2`, conditioned on `condi`.
#'
#' @examples
#' mat1 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' mat2 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' condi <- rnorm(100)
#' CMImat2mat(mat1, mat2, condi)
#'
#' @export
CMImat2mat <- function(mat1, mat2, condi, bin = NULL, sp_order = NULL) {
  # Validate inputs as non-empty numeric vectors

  if (!is.numeric(condi) || length(condi) == 0) {
    stop("Input condi must be a non-empty numeric vector.")
  }

  if (!is.matrix(mat1) || !is.numeric(mat1) || is.na(sum(mat1)) ) {
    stop("Input mat must be a non-NA numeric matrix.")
  }

  if (!is.matrix(mat2) || !is.numeric(mat2) || is.na(sum(mat2)) ) {
    stop("Input mat must be a non-NA numeric matrix.")
  }

  # Check that 2 mats are of the same size
  if (ncol(mat1) != ncol(mat2) | nrow(mat1) != nrow(mat2) ) {
    stop("The size of mat1 and mat2 must be the same.")
  }

  if (ncol(mat1) != length(condi)) {
    stop("The column of mat1 and the length of vector condi must be the same.")
  }

  max_MI <- sapply(1:nrow(mat1), function(x)
    max(getMI(mat1[x,], mat1[x,], bin = bin, sp_order = sp_order),
        getMI(mat2[x,], mat2[x,], bin = bin, sp_order = sp_order))
  )

  mi_vector <- sapply(1:nrow(mat1), function(x)
    getCMI(mat1[x,], mat2[x,], condi, bin = bin, sp_order = sp_order))

  normalized_mi_vector <- ifelse(max_MI == 0, 0, mi_vector / max_MI)
  names(normalized_mi_vector) <- row.names(mat1)
  normalized_mi_vector[normalized_mi_vector<0] = 0
  return(normalized_mi_vector)
}

#' @title Normalized Conditional Mutual Information Between Matrix and Vector Given Two Conditions
#' @description Computes the normalized conditional mutual information (CMI) between each row of a matrix and a vector,
#' given two condition vectors, normalized by the mutual information of the vector with itself using the specified bins and spline order.
#'
#' @param mat A numeric matrix. For example, each row represents a gene and each column represents a sample.
#' @param vec A numeric vector, with length equal to the number of columns in `mat`.
#' @param condi1 A numeric condition vector, matching the number of columns in `mat`.
#' @param condi2 Another numeric condition vector, matching the number of columns in `mat`.
#' @param bin An integer specifying the number of bins. Default is NULL.
#' @param sp_order An integer specifying the spline order. Must be less than `bin`. Default is NULL.
#' @return A numeric vector representing the normalized conditional mutual information (CMI) between each row of `mat` and `vec`, given `condi1` and `condi2`.
#'
#' @examples
#' mat <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' vec <- rnorm(100)
#' condi1 <- rnorm(100)
#' condi2 <- rnorm(100)
#' CMIBiCondimat2vec(mat, vec, condi1, condi2)
#'
#' @export
CMIBiCondimat2vec <- function(mat, vec, condi1, condi2, bin = NULL, sp_order = NULL) {
  # Validate inputs as non-empty numeric vectors

  if (!is.numeric(vec) || length(vec) == 0) {
    stop("Input vec must be a non-empty numeric vector.")
  }

  if (!is.numeric(condi1) || length(condi1) == 0) {
    stop("Input condi1 must be a non-empty numeric vector.")
  }

  if (!is.numeric(condi2) || length(condi2) == 0) {
    stop("Input condi2 must be a non-empty numeric vector.")
  }

  if (!is.matrix(mat) || !is.numeric(mat) || is.na(sum(mat)) ) {
    stop("Input mat must be a non-NA numeric matrix.")
  }

  # Check that mat, vec and condi1 are of the same length
  if (ncol(mat) != length(vec)) {
    stop("The column of mat and the length of vec must be the same.")
  }

  if (ncol(mat) != length(condi1)) {
    stop("The column of mat and the length of vector condi1 must be the same.")
  }

  if (ncol(mat) != length(condi2)) {
    stop("The column of mat and the length of vector condi2 must be the same.")
  }

  max_MI <- getMI(vec, vec, bin = bin, sp_order = sp_order)
  # Ensure vec contains information
  if (max_MI == 0) {
    stop("Vector 'vec' must contain information.")
  }

  condi_mi_vector <- sapply(1:nrow(mat), function(x)
    getCMIBiCondi(mat[x,], vec, condi1, condi2, bin = bin, sp_order = sp_order))

  normalized_condi_mi_vector <- condi_mi_vector / max_MI
  names(normalized_condi_mi_vector) <- row.names(mat)

  return(normalized_condi_mi_vector)
}


#' @title Normalized Conditional Mutual Information Between Two Matrices Given Two Conditions
#' @description Computes the normalized conditional mutual information (CMI) between corresponding rows of two matrices,
#' given two condition variables, normalized by their individual information content. CMI is calculated using the specified number of bins and spline order.
#'
#' @param mat1 A numeric matrix. For example, each row represents a gene and each column represents a sample.
#' @param mat2 Another numeric matrix to compare against. Must have the same dimensions as `mat1`.
#' @param condi1 A numeric condition vector, matching the number of columns in `mat1`.
#' @param condi2 Another numeric condition vector, matching the number of columns in `mat1`.
#' @param bin An integer specifying the number of bins. Default is NULL.
#' @param sp_order An integer specifying the spline order. Must be less than `bin`. Default is NULL.
#' @return A numeric vector representing the normalized conditional mutual information (CMI) between pairs of rows
#' from `mat1` and `mat2`, conditioned on `condi1` and `condi2`.
#'
#' @examples
#' mat1 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' mat2 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' condi1 <- rnorm(100)
#' condi2 <- rnorm(100)
#' CMIBiCondimat2mat(mat1, mat2, condi1, condi2)
#'
#' @export
CMIBiCondimat2mat <- function(mat1, mat2, condi1, condi2, bin = NULL, sp_order = NULL) {
  # Validate inputs as non-empty numeric vectors

  if (!is.numeric(condi1) || length(condi1) == 0) {
    stop("Input condi1 must be a non-empty numeric vector.")
  }

  if (!is.numeric(condi2) || length(condi2) == 0) {
    stop("Input condi2 must be a non-empty numeric vector.")
  }

  if (!is.matrix(mat1) || !is.numeric(mat1) || is.na(sum(mat1)) ) {
    stop("Input mat1 must be a non-NA numeric matrix.")
  }

  if (!is.matrix(mat2) || !is.numeric(mat2) || is.na(sum(mat2)) ) {
    stop("Input mat2 must be a non-NA numeric matrix.")
  }

  # Check that 2 mats are of the same size
  if (ncol(mat1) != ncol(mat2) | nrow(mat1) != nrow(mat2) ) {
    stop("The size of mat1 and mat2 must be the same.")
  }

  if (ncol(mat1) != length(condi1)) {
    stop("The column of mat1 and the length of vector condi1 must be the same.")
  }

  if (ncol(mat1) != length(condi2)) {
    stop("The column of mat1 and the length of vector condi2 must be the same.")
  }

  max_MI <- sapply(1:nrow(mat1), function(x)
    max(getMI(mat1[x,], mat1[x,], bin = bin, sp_order = sp_order),
        getMI(mat2[x,], mat2[x,], bin = bin, sp_order = sp_order))
  )

  mi_vector <- sapply(1:nrow(mat1), function(x)
    getCMIBiCondi(mat1[x,], mat2[x,], condi1, condi2, bin = bin, sp_order = sp_order))

  normalized_mi_vector <- ifelse(max_MI == 0, 0, mi_vector / max_MI)
  names(normalized_mi_vector) <- row.names(mat1)
  normalized_mi_vector[normalized_mi_vector<0] = 0
  return(normalized_mi_vector)
}
