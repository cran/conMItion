#' @title Internal Utility Function for Random Sampling Loops
#' @description Executes nested loops to perform random sampling operations.
#' This function is intended for internal use only and should not be accessed by package users.
#' @keywords internal
#'
#' @return No return value. This function is designed to allow large permutation tests distributed to multiple nodes.
additionalRunIF <- function(x, y, z, p){
  for (i in 1:x) {
    for (j in 1:y) {
      for (k in 1:z) {
        runif(p)
      }
    }
  }
}

#' @title Internal Utility Function for Within-Bin Permutation
#' @description Permutes the values of one vector within bins defined by another vector.
#' This function is intended for internal use only and should not be accessed by package users.
#'
#' The binning vector is split into evenly spaced bins between its minimum and maximum
#' values. Values of the target vector are then randomly permuted only within each bin.
#' This preserves the approximate relationship between the target vector and the binning
#' variable while breaking the original sample-wise ordering within local strata.
#'
#' @param x A numeric vector used to define bins.
#' @param y A vector with the same length as `x`. Values in `y` are permuted within bins of `x`.
#' @param bin An integer specifying the number of evenly spaced bins. Default is 8.
#'
#' @return A vector with the same length as `y`, where values have been permuted within bins of `x`.
#'
#' @keywords internal
permute_within_bins <- function(x, y, bin = 8) {

  if (!is.numeric(x) || length(x) == 0) {
    stop("Input x must be a non-empty numeric vector.")
  }

  if (length(x) != length(y)) {
    stop("Input x and y must have the same length.")
  }

  if (!is.numeric(bin) || length(bin) != 1 || is.na(bin) || bin < 1) {
    stop("Input bins must be a positive integer.")
  }

  bin <- as.integer(bin)

  if (all(is.na(x))) {
    stop("Input x cannot contain only NA values.")
  }

  y_perm <- y

  ## If x has only one unique non-NA value, binning is impossible.
  ## In this case, all values belong to one stratum, so ordinary permutation is used.
  if (min(x, na.rm = TRUE) == max(x, na.rm = TRUE)) {
    return(sample(y, length(y), replace = FALSE))
  }

  bins <- cut(
    x,
    breaks = seq(
      min(x, na.rm = TRUE),
      max(x, na.rm = TRUE),
      length.out = bin + 1
    ),
    include.lowest = TRUE,
    labels = FALSE
  )

  for (b in unique(stats::na.omit(bins))) {
    idx <- which(bins == b)
    y_perm[idx] <- sample(y[idx], length(idx), replace = FALSE)
  }

  return(y_perm)
}


#' @title Internal Utility Function for Within-Bin Permutation With Two Conditions
#' @description Permutes the values of one vector within joint bins defined by two condition vectors.
#' This function is intended for internal use only and should not be accessed by package users.
#'
#' The two condition vectors are each split into quantile-based bins. Joint bins are then
#' defined by the combination of the two bin assignments. Values of the target vector are
#' randomly permuted only within each joint bin. This preserves the approximate relationship
#' between the target vector and both condition variables while breaking the original
#' sample-wise ordering within local strata.
#'
#' @param y A vector to be permuted.
#' @param condi1 A numeric vector used to define the first set of bins.
#' @param condi2 A numeric vector used to define the second set of bins.
#' @param bin An integer specifying the number of bins for each condition vector. Default is 8.
#'
#' @return A vector with the same length as `y`, where values have been permuted within
#' joint bins defined by `condi1` and `condi2`.
#'
#' @keywords internal
permute_within_two_condition_bins <- function(y, condi1, condi2, bin = 8) {
  if (length(y) != length(condi1) || length(y) != length(condi2)) {
    stop("y, condi1, and condi2 must have the same length.")
  }

  if (is.null(bin)) {
    stop("bin must be specified for within-bin permutation.")
  }

  if (!is.numeric(condi1) || !is.numeric(condi2)) {
    stop("condi1 and condi2 must be numeric vectors.")
  }

  if (anyNA(y) || anyNA(condi1) || anyNA(condi2)) {
    stop("y, condi1, and condi2 must not contain NA values.")
  }

  bin1 <- cut(condi1,
              breaks = seq(
                min(condi1, na.rm = TRUE),
                max(condi1, na.rm = TRUE),
                length.out = bin + 1),
    include.lowest = TRUE,
    labels = FALSE
  )

  bin2 <- cut(condi2,
    breaks = seq(
      min(condi2, na.rm = TRUE),
      max(condi2, na.rm = TRUE),
      length.out = bin + 1),
    include.lowest = TRUE,
    labels = FALSE
  )

  joint_bins <- interaction(bin1, bin2, drop = TRUE)

  y_perm <- y

  for (b in unique(joint_bins)) {
    idx <- which(joint_bins == b)

    if (length(idx) > 1) {
      y_perm[idx] <- sample(y[idx], length(idx), replace = FALSE)
    }
  }

  return(y_perm)
}


#' @title Permuted Correlation Between Matrix and Vector
#' @description Computes the correlation between a randomly sampled vector from a matrix and a given vector.
#' The sampling is done multiple times to generate a distribution.
#'
#' @param mat A numeric matrix. For example, each row represents a gene and each column represents a sample.
#' @param vec A numeric vector, with length equal to the number of columns in `mat`.
#' @param cor_type Type of correlation to calculate: "Pearson", "Kendall", or "Spearman". Default is "Pearson".
#' @param bulkIdx Index to divide the task when processing many permutations. Default is 0.
#' @param permutationTimes Number of permutations for sampling. Default is 1000.
#' @param seedNum Seed for random number generation. Default is 99999999.
#' @return A numeric vector of correlation values for each permutation.
#'
#' @examples
#' mat <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' vec <- rnorm(100)
#' CORmat2vecPermu(mat, vec)
#'
#' @export
CORmat2vecPermu <- function(mat, vec, cor_type = "pearson", bulkIdx = 0, permutationTimes = 1000, seedNum = 99999999){
  set.seed(seedNum)

  if (!is.numeric(vec) || length(vec) == 0) {
    stop("Input vec must be a non-empty numeric vector.")
  }

  if (!is.matrix(mat) || !is.numeric(mat) || is.na(sum(mat)) ) {
    stop("Input mat must be a non-NA numeric matrix.")
  }

  # Check that mat and vec are of the same length
  if (ncol(mat) != length(vec)) {
    stop("The column of mat and the length of vec must be the same.")
  }

  if (bulkIdx!=0) {
    additionalRunIF(bulkIdx,ncol(mat),2,permutationTimes)
  }

  permutationCorVector <- rep(0,permutationTimes)
  for (p in 1:permutationTimes) {
    randRowIdx <- sample(1:nrow(mat),size = ncol(mat),replace = (ncol(mat)>nrow(mat)) )
    randColIdx <- sample(1:ncol(mat),replace = F)
    idx <- cbind(randRowIdx, randColIdx)
    randomExp1 <- mat[idx]
    permutationCorVector[p] <- cor(randomExp1, vec, method = cor_type)
  }

  return(permutationCorVector)
}

#' @title Permuted Normalized Mutual Information Between Matrix and Vector
#' @description Computes the mutual information (MI) between a random vector sampled from a matrix and a vector,
#' normalized by the mutual information of the vector with itself. The sampling is done multiple times to generate a distribution.
#'
#' @param mat A numeric matrix. For example, each row represents a gene and each column represents a sample.
#' @param vec A numeric vector, with length equal to the number of columns in `mat`.
#' @param bin An integer specifying the number of bins. Default is NULL.
#' @param sp_order An integer specifying the spline order. Must be less than `bin`. Default is NULL.
#' @param bulkIdx Index to divide the task when processing many permutations. Default is 0.
#' @param permutationTimes Number of permutations for sampling. Default is 1000.
#' @param seedNum Seed for random number generation. Default is 99999999.
#' @return A numeric vector of normalized mutual information (MI) values for each permutation.
#'
#' @examples
#' mat <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' vec <- rnorm(100)
#' MImat2vecPermu(mat, vec)
#'
#' @export
MImat2vecPermu <- function(mat, vec, bin = NULL, sp_order = NULL, bulkIdx = 0, permutationTimes = 1000, seedNum = 99999999){
  set.seed(seedNum)

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

  max_MI = getMI(vec, vec, bin = bin, sp_order = sp_order)
  # Ensure vec contains information
  if (max_MI == 0) {
    stop("Vector 'vec' must contain information.")
  }

  if (bulkIdx!=0) {
    additionalRunIF(bulkIdx,ncol(mat),2,permutationTimes)
  }

  permutationMIVector <- rep(0,permutationTimes)
  for (p in 1:permutationTimes) {
    randRowIdx <- sample(1:nrow(mat),size = ncol(mat),replace = (ncol(mat)>nrow(mat)) )
    randColIdx <- sample(1:ncol(mat),replace = F)
    idx <- cbind(randRowIdx, randColIdx)
    randomExp1 <- mat[idx]

    permutationMIVector[p] <- getMI(randomExp1, vec, bin = bin, sp_order = sp_order)
  }
  permutationMIVector <- permutationMIVector/max_MI
  permutationMIVector[permutationMIVector<0] = 0
  return(permutationMIVector)
}

#' @title Permuted Normalized Conditional Mutual Information Between Matrix and Vector
#' @description Computes the conditional mutual information (CMI) between a random vector sampled from a matrix and a vector,
#' conditioned on a third vector, normalized by the mutual information of the vector with itself.
#' The sampling is done multiple times to generate a distribution.
#'
#' @param mat A numeric matrix. For example, each row represents a gene and each column represents a sample.
#' @param vec A numeric vector, with length equal to the number of columns in `mat`.
#' @param condi A numeric condition vector, matching the number of columns in `mat`.
#' @param bin An integer specifying the number of bins. Default is NULL.
#' @param sp_order An integer specifying the spline order. Must be less than `bin`. Default is NULL.
#' @param bulkIdx Index to divide the task when processing many permutations. Default is 0.
#' @param permutationTimes Number of permutations for sampling. Default is 1000.
#' @param seedNum Seed for random number generation. Default is 99999999.
#' @param permutation_method A character string specifying the permutation scheme. Use `"random"` for the original random row-column sampling procedure, or `"within_bins"` for conditional permutation within bins of `condi`.
#' @return A numeric vector of normalized conditional mutual information (CMI) values for each permutation.
#'
#' @examples
#' mat <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' vec <- rnorm(100)
#' condi <- rnorm(100)
#' CMImat2vecPermu(mat, vec, condi)
#'
#' @export
CMImat2vecPermu <- function(mat, vec, condi, bin = NULL, sp_order = NULL, bulkIdx = 0, permutationTimes = 1000, seedNum = 99999999, permutation_method = c("random", "within_bins")){
  set.seed(seedNum)

  permutation_method <- match.arg(permutation_method)

  if (!is.numeric(vec) || length(vec) == 0) {
    stop("Input vec must be a non-empty numeric vector.")
  }

  if (!is.numeric(condi) || length(condi) == 0) {
    stop("Input condi must be a non-empty numeric vector.")
  }

  if (!is.matrix(mat) || !is.numeric(mat) || is.na(sum(mat)) ) {
    stop("Input mat must be a non-NA numeric matrix.")
  }

  # Check that mat and vec are of the same length
  if (ncol(mat) != length(vec)) {
    stop("The column of mat and the length of vec must be the same.")
  }

  if (ncol(mat) != length(condi)) {
    stop("The column of mat and the length of vector condi must be the same.")
  }

  n <- length(condi)
  # Check parameters and generate automatically if incorrect
  params <- resolve_params(n, bin, sp_order)

  bin=as.integer(params$bin)
  sp_order=as.integer(params$sp_order)

  max_MI = getMI(vec, vec, bin = bin, sp_order = sp_order)
  if (max_MI == 0) {
    stop("Vector 'vec' must contain information.")
  }

  if (bulkIdx!=0) {
    additionalRunIF(bulkIdx,ncol(mat),2,permutationTimes)
  }

  permutationCondiMIVector <- rep(0,permutationTimes)
  for (p in 1:permutationTimes) {
    randRowIdx <- sample(1:nrow(mat),size = ncol(mat),replace = (ncol(mat)>nrow(mat)) )
    if (permutation_method == "random") {
      randColIdx <- sample(1:ncol(mat),replace = F)
      idx <- cbind(randRowIdx, randColIdx)
      randomExp1 <- mat[idx]

    } else if (permutation_method == "within_bins") {
      colIdx <- 1:ncol(mat)
      idx <- cbind(randRowIdx, colIdx)
      randomExp1_raw <- mat[idx]

      randomExp1 <- permute_within_bins(x = condi,y = randomExp1_raw, bin = bin)
    }
    permutationCondiMIVector[p] <- getCMI(randomExp1, vec, condi, bin = bin, sp_order = sp_order)
  }
  permutationCondiMIVector <- permutationCondiMIVector/max_MI
  permutationCondiMIVector[permutationCondiMIVector<0] = 0
  return(permutationCondiMIVector)
}

#' @title Permuted Mutual Information Between Two Matrices
#' @description Computes the normalized mutual information (MI) between vectors sampled from two matrices normalized by the individual information content.
#' The sampling is done multiple times to generate a distribution.
#'
#' @param mat1 A numeric matrix. For example, each row represents a gene and each column represents a sample.
#' @param mat2 Another numeric matrix to compare against. Must have the same dimensions as `mat1`.
#' @param bin An integer specifying the number of bins. Default is NULL.
#' @param sp_order An integer specifying the spline order. Must be less than `bin`. Default is NULL.
#' @param bulkIdx Index to divide the task when processing many permutations. Default is 0.
#' @param permutationTimes Number of permutations for sampling. Default is 1000.
#' @param seedNum Seed for random number generation. Default is 99999999.
#' @return A numeric vector of normalized mutual information (MI) values for each permutation.
#'
#' @examples
#' mat1 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' mat2 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' MImat2matPermu(mat1, mat2)
#'
#' @export
MImat2matPermu <- function(mat1, mat2, bin = NULL, sp_order = NULL, bulkIdx = 0, permutationTimes = 1000, seedNum = 99999999){
  set.seed(seedNum)

  # Validate inputs as non-empty numeric vectors

  if (!is.matrix(mat1) || !is.numeric(mat1) || is.na(sum(mat1)) ) {
    stop("Input mat1 must be a non-NA numeric matrix.")
  }

  if (!is.matrix(mat2) || !is.numeric(mat2) || is.na(sum(mat2)) ) {
    stop("Input mat2 must be a non-NA numeric matrix.")
  }

  # Check that 2 mats are of the same size
  if (ncol(mat1) != ncol(mat2) ) {
    stop("The column size of mat1 and mat2 must be the same.")
  }

  if (bulkIdx!=0) {
    additionalRunIF(bulkIdx,ncol(mat1),4,permutationTimes)
  }

  permutationMIVector <- rep(0,permutationTimes)

  for (p in 1:permutationTimes) {
    randRowIdx <- sample(1:nrow(mat1),size = ncol(mat1),replace = (ncol(mat1)>nrow(mat1)) )
    randColIdx <- sample(1:ncol(mat1),replace = F)
    idx <- cbind(randRowIdx, randColIdx)
    randomExp1 <- mat1[idx]
    randRowIdx <- sample(1:nrow(mat2),size = ncol(mat2),replace = (ncol(mat2)>nrow(mat2)) )
    randColIdx <- sample(1:ncol(mat2),replace = F)
    idx <- cbind(randRowIdx, randColIdx)
    randomExp2 <- mat2[idx]
    permutationMIVector[p] <- getMI(randomExp1,randomExp2, bin = bin, sp_order = sp_order) /
      max(getMI(randomExp1,randomExp1, bin = bin, sp_order = sp_order),
          getMI(randomExp2,randomExp2, bin = bin, sp_order = sp_order))

  }
  permutationMIVector[permutationMIVector<0] = 0
  return(permutationMIVector)
}

#' @title Permuted Conditional Mutual Information Between Two Matrices
#' @description Computes the normalized conditional mutual information (CMI) between vectors sampled from two matrices,
#' conditioned on another vector, normalized by the individual information content. The sampling is done multiple times to generate a distribution.
#'
#' @param mat1 A numeric matrix. For example, each row represents a gene and each column represents a sample.
#' @param mat2 Another numeric matrix to compare against. Must have the same columns as `mat1`.
#' @param condi A numeric condition vector, matching the number of columns in `mat1`.
#' @param bin An integer specifying the number of bins. Default is NULL.
#' @param sp_order An integer specifying the spline order. Must be less than `bin`. Default is NULL.
#' @param bulkIdx Index to divide the task when processing many permutations. Default is 0.
#' @param permutationTimes Number of permutations for sampling. Default is 1000.
#' @param seedNum Seed for random number generation. Default is 99999999.
#' @param permutation_method A character string specifying the permutation scheme. Use `"random"` for the original random row-column sampling procedure, or `"within_bins"` for conditional permutation within bins of `condi`.
#' @return A numeric vector of normalized conditional mutual information (CMI) values for each permutation.
#'
#' @examples
#' mat1 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' mat2 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' condi <- rnorm(100)
#' CMImat2matPermu(mat1, mat2, condi)
#'
#' @export
CMImat2matPermu <- function(mat1, mat2, condi, bin = NULL, sp_order = NULL, bulkIdx = 0, permutationTimes = 1000, seedNum = 99999999, permutation_method = c("random", "within_bins")){
  set.seed(seedNum)

  permutation_method <- match.arg(permutation_method)
  # Validate inputs as non-empty numeric vectors

  if (!is.numeric(condi) || length(condi) == 0) {
    stop("Input condi must be a non-empty numeric vector.")
  }

  if (!is.matrix(mat1) || !is.numeric(mat1) || is.na(sum(mat1)) ) {
    stop("Input mat1 must be a non-NA numeric matrix.")
  }

  if (!is.matrix(mat2) || !is.numeric(mat2) || is.na(sum(mat2)) ) {
    stop("Input mat2 must be a non-NA numeric matrix.")
  }

  # Check that 2 mats are of the same size
  if (ncol(mat1) != ncol(mat2) ) {
    stop("The column size of mat1 and mat2 must be the same.")
  }

  if (ncol(mat1) != length(condi)) {
    stop("The column of mat1 and the length of vector condi must be the same.")
  }

  n <- length(condi)
  # Check parameters and generate automatically if incorrect
  params <- resolve_params(n, bin, sp_order)

  bin=as.integer(params$bin)
  sp_order=as.integer(params$sp_order)

  if (bulkIdx!=0) {
    additionalRunIF(bulkIdx,ncol(mat1),4,permutationTimes)
  }

  permutationMIVector <- rep(0,permutationTimes)

  for (p in 1:permutationTimes) {
    if (permutation_method == "random") {
      randRowIdx <- sample(1:nrow(mat1),size = ncol(mat1),replace = (ncol(mat1)>nrow(mat1)) )
      randColIdx <- sample(1:ncol(mat1),replace = F)
      idx <- cbind(randRowIdx, randColIdx)
      randomExp1 <- mat1[idx]

      randRowIdx <- sample(1:nrow(mat2),size = ncol(mat2),replace = (ncol(mat2)>nrow(mat2)) )
      randColIdx <- sample(1:ncol(mat2),replace = F)
      idx <- cbind(randRowIdx, randColIdx)
      randomExp2 <- mat2[idx]
    } else if (permutation_method == "within_bins") {
      randRowIdx <- sample(1:nrow(mat1),size = ncol(mat1),replace = (ncol(mat1)>nrow(mat1)) )
      colIdx <- 1:ncol(mat1)
      idx <- cbind(randRowIdx, colIdx)
      randomExp1 <- mat1[idx]

      randRowIdx <- sample(1:nrow(mat2),size = ncol(mat2),replace = (ncol(mat2)>nrow(mat2)) )
      colIdx <- 1:ncol(mat2)
      idx <- cbind(randRowIdx, colIdx)
      randomExp2_raw <- mat2[idx]
      randomExp2 <- permute_within_bins(x = condi,y = randomExp2_raw, bin = bin)
    }

    permutationMIVector[p] <- getCMI(randomExp1,randomExp2, condi, bin = bin, sp_order = sp_order) /
      max(getMI(randomExp1,randomExp1, bin = bin, sp_order = sp_order),
          getMI(randomExp2,randomExp2, bin = bin, sp_order = sp_order))

  }
  permutationMIVector[permutationMIVector<0] = 0
  return(permutationMIVector)
}


#' @title Permuted Normalized Conditional Mutual Information Between Matrix and Vector Given Two Conditions
#' @description Computes the conditional mutual information (CMI) between a random vector sampled from a matrix and a vector,
#' conditioned on two vectors, normalized by the mutual information of the vector with itself.
#' The sampling is done multiple times to generate a distribution.
#'
#' @param mat A numeric matrix. For example, each row represents a gene and each column represents a sample.
#' @param vec A numeric vector, with length equal to the number of columns in `mat`.
#' @param condi1 A numeric condition vector, matching the number of columns in `mat`.
#' @param condi2 Another numeric condition vector, matching the number of columns in `mat`.
#' @param bin An integer specifying the number of bins. Default is NULL.
#' @param sp_order An integer specifying the spline order. Must be less than `bin`. Default is NULL.
#' @param bulkIdx Index to divide the task when processing many permutations. Default is 0.
#' @param permutationTimes Number of permutations for sampling. Default is 1000.
#' @param seedNum Seed for random number generation. Default is 99999999.
#' @param permutation_method A character string specifying the permutation scheme. Use `"random"` for the original random row-column sampling procedure, or `"within_bins"` for conditional permutation within joint bins of `condi1` and `condi2`.
#' @return A numeric vector of normalized conditional mutual information (CMI) values for each permutation.
#'
#' @examples
#' mat <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' vec <- rnorm(100)
#' condi1 <- rnorm(100)
#' condi2 <- rnorm(100)
#' CMIBiCondimat2vecPermu(mat, vec, condi1, condi2)
#'
#' @export
CMIBiCondimat2vecPermu <- function(mat, vec, condi1, condi2, bin = NULL, sp_order = NULL, bulkIdx = 0, permutationTimes = 1000, seedNum = 99999999,permutation_method = c("random", "within_bins")){
  set.seed(seedNum)

  permutation_method <- match.arg(permutation_method)

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

  # Check that mat and vec are of the same length
  if (ncol(mat) != length(vec)) {
    stop("The column of mat and the length of vec must be the same.")
  }

  if (ncol(mat) != length(condi1)) {
    stop("The column of mat and the length of vector condi1 must be the same.")
  }

  if (ncol(mat) != length(condi2)) {
    stop("The column of mat and the length of vector condi2 must be the same.")
  }

  n <- length(condi1)
  # Check parameters and generate automatically if incorrect
  params <- resolve_params(n, bin, sp_order)

  bin=as.integer(params$bin)
  sp_order=as.integer(params$sp_order)

  max_MI = getMI(vec, vec, bin = bin, sp_order = sp_order)
  if (max_MI == 0) {
    stop("Vector 'vec' must contain information.")
  }

  if (bulkIdx!=0) {
    additionalRunIF(bulkIdx,ncol(mat),2,permutationTimes)
  }

  permutationCondiMIVector <- rep(0,permutationTimes)
  for (p in 1:permutationTimes) {
    randRowIdx <- sample(1:nrow(mat),size = ncol(mat),replace = (ncol(mat)>nrow(mat)) )
    if (permutation_method == "random") {
      randColIdx <- sample(1:ncol(mat), replace = FALSE)
      idx <- cbind(randRowIdx, randColIdx)
      randomExp1 <- mat[idx]
    } else if (permutation_method == "within_bins") {
      colIdx <- 1:ncol(mat)
      idx <- cbind(randRowIdx, colIdx)

      randomExp1_raw <- mat[idx]
      randomExp1 <- permute_within_two_condition_bins(
        y = randomExp1_raw,
        condi1 = condi1,
        condi2 = condi2,
        bin = bin
      )
    }
    permutationCondiMIVector[p] <- getCMIBiCondi(randomExp1, vec, condi1, condi2, bin = bin, sp_order = sp_order)
  }
  permutationCondiMIVector <- permutationCondiMIVector/max_MI
  permutationCondiMIVector[permutationCondiMIVector<0] = 0
  return(permutationCondiMIVector)
}


#' @title Permuted Conditional Mutual Information Between Two Matrices Given Two Conditions
#' @description Computes the normalized conditional mutual information (CMI) between vectors sampled from two matrices,
#' conditioned on two vectors, normalized by the individual information content.
#' The sampling is done multiple times to generate a distribution.
#'
#' @param mat1 A numeric matrix. For example, each row represents a gene and each column represents a sample.
#' @param mat2 Another numeric matrix to compare against. Must have the same dimensions as `mat1`.
#' @param condi1 A numeric condition vector, matching the number of columns in `mat1`.
#' @param condi2 Another numeric condition vector, matching the number of columns in `mat`.
#' @param bin An integer specifying the number of bins. Default is NULL.
#' @param sp_order An integer specifying the spline order. Must be less than `bin`. Default is NULL.
#' @param bulkIdx Index to divide the task when processing many permutations. Default is 0.
#' @param permutationTimes Number of permutations for sampling. Default is 1000.
#' @param seedNum Seed for random number generation. Default is 99999999.
#' @param permutation_method A character string specifying the permutation scheme. Use `"random"` for the original random row-column sampling procedure, or `"within_bins"` for conditional permutation within joint bins of `condi1` and `condi2`.
#' @return A numeric vector of normalized conditional mutual information (CMI) values for each permutation.
#'
#' @examples
#' mat1 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' mat2 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' condi1 <- rnorm(100)
#' condi2 <- rnorm(100)
#' CMIBiCondimat2matPermu(mat1, mat2, condi1, condi2)
#'
#' @export
CMIBiCondimat2matPermu <- function(mat1, mat2, condi1, condi2, bin = NULL, sp_order = NULL, bulkIdx = 0, permutationTimes = 1000, seedNum = 99999999,permutation_method = c("random", "within_bins")){
  set.seed(seedNum)

  permutation_method <- match.arg(permutation_method)

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
  if (ncol(mat1) != ncol(mat2) ) {
    stop("The column size of mat1 and mat2 must be the same.")
  }

  if (ncol(mat1) != length(condi1)) {
    stop("The column of mat1 and the length of vector condi1 must be the same.")
  }

  if (ncol(mat1) != length(condi2)) {
    stop("The column of mat1 and the length of vector condi2 must be the same.")
  }

  n <- length(condi1)
  # Check parameters and generate automatically if incorrect
  params <- resolve_params(n, bin, sp_order)

  bin=as.integer(params$bin)
  sp_order=as.integer(params$sp_order)

  if (bulkIdx!=0) {
    additionalRunIF(bulkIdx,ncol(mat1),4,permutationTimes)
  }

  permutationMIVector <- rep(0,permutationTimes)

  for (p in 1:permutationTimes) {

    if (permutation_method == "random") {
      randRowIdx <- sample(1:nrow(mat1),size = ncol(mat1),replace = (ncol(mat1)>nrow(mat1)) )
      randColIdx <- sample(1:ncol(mat1),replace = F)
      idx <- cbind(randRowIdx, randColIdx)
      randomExp1 <- mat1[idx]
      randRowIdx <- sample(1:nrow(mat2),size = ncol(mat2),replace = (ncol(mat2)>nrow(mat2)) )
      randColIdx <- sample(1:ncol(mat2),replace = F)
      idx <- cbind(randRowIdx, randColIdx)
      randomExp2 <- mat2[idx]
    } else if (permutation_method == "within_bins") {
      randRowIdx <- sample(1:nrow(mat1),size = ncol(mat1),replace = (ncol(mat1) > nrow(mat1)))
      colIdx <- seq_len(ncol(mat1))
      idx <- cbind(randRowIdx, colIdx)
      randomExp1 <- mat1[idx]

      randRowIdx <- sample(1:nrow(mat2),size = ncol(mat2),replace = (ncol(mat2) > nrow(mat2)))
      colIdx <- seq_len(ncol(mat2))
      idx <- cbind(randRowIdx, colIdx)
      randomExp2_raw <- mat2[idx]

      randomExp2 <- permute_within_two_condition_bins(y = randomExp2_raw,
                                                      condi1 = condi1,condi2 = condi2,
                                                      bin = bin)
    }

    permutationMIVector[p] <- getCMIBiCondi(randomExp1,randomExp2, condi1, condi2, bin = bin, sp_order = sp_order) /
      max(getMI(randomExp1,randomExp1, bin = bin, sp_order = sp_order),
          getMI(randomExp2,randomExp2, bin = bin, sp_order = sp_order))

  }
  permutationMIVector[permutationMIVector<0] = 0
  return(permutationMIVector)
}
