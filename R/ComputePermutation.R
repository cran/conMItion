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
#' @param bin An integer specifying the number of bins. Default is 6.
#' @param sp_order An integer specifying the spline order. Must be less than `bin`. Default is 2.
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
MImat2vecPermu <- function(mat, vec, bin = 6, sp_order = 2, bulkIdx = 0, permutationTimes = 1000, seedNum = 99999999){
  set.seed(seedNum)

  bin <- as.integer(bin)
  sp_order <- as.integer(sp_order)

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
#' @param bin An integer specifying the number of bins. Default is 6.
#' @param sp_order An integer specifying the spline order. Must be less than `bin`. Default is 2.
#' @param bulkIdx Index to divide the task when processing many permutations. Default is 0.
#' @param permutationTimes Number of permutations for sampling. Default is 1000.
#' @param seedNum Seed for random number generation. Default is 99999999.
#' @return A numeric vector of normalized conditional mutual information (CMI) values for each permutation.
#'
#' @examples
#' mat <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' vec <- rnorm(100)
#' condi <- rnorm(100)
#' CMImat2vecPermu(mat, vec, condi)
#'
#' @export
CMImat2vecPermu <- function(mat, vec, condi, bin = 6, sp_order = 2, bulkIdx = 0, permutationTimes = 1000, seedNum = 99999999){
  set.seed(seedNum)

  bin <- as.integer(bin)
  sp_order <- as.integer(sp_order)

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
    randColIdx <- sample(1:ncol(mat),replace = F)
    idx <- cbind(randRowIdx, randColIdx)
    randomExp1 <- mat[idx]

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
#' @param bin An integer specifying the number of bins. Default is 6.
#' @param sp_order An integer specifying the spline order. Must be less than `bin`. Default is 2.
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
MImat2matPermu <- function(mat1, mat2, bin = 6, sp_order = 2, bulkIdx = 0, permutationTimes = 1000, seedNum = 99999999){
  set.seed(seedNum)

  # Validate inputs as non-empty numeric vectors
  bin <- as.integer(bin)
  sp_order <- as.integer(sp_order)

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
#' @param mat2 Another numeric matrix to compare against. Must have the same dimensions as `mat1`.
#' @param condi A numeric condition vector, matching the number of columns in `mat1`.
#' @param bin An integer specifying the number of bins. Default is 6.
#' @param sp_order An integer specifying the spline order. Must be less than `bin`. Default is 2.
#' @param bulkIdx Index to divide the task when processing many permutations. Default is 0.
#' @param permutationTimes Number of permutations for sampling. Default is 1000.
#' @param seedNum Seed for random number generation. Default is 99999999.
#' @return A numeric vector of normalized conditional mutual information (CMI) values for each permutation.
#'
#' @examples
#' mat1 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' mat2 <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' condi <- rnorm(100)
#' CMImat2matPermu(mat1, mat2, condi)
#'
#' @export
CMImat2matPermu <- function(mat1, mat2, condi, bin = 6, sp_order = 2, bulkIdx = 0, permutationTimes = 1000, seedNum = 99999999){
  set.seed(seedNum)

  # Validate inputs as non-empty numeric vectors
  bin <- as.integer(bin)
  sp_order <- as.integer(sp_order)

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
#' @param bin An integer specifying the number of bins. Default is 6.
#' @param sp_order An integer specifying the spline order. Must be less than `bin`. Default is 2.
#' @param bulkIdx Index to divide the task when processing many permutations. Default is 0.
#' @param permutationTimes Number of permutations for sampling. Default is 1000.
#' @param seedNum Seed for random number generation. Default is 99999999.
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
CMIBiCondimat2vecPermu <- function(mat, vec, condi1, condi2, bin = 6, sp_order = 2, bulkIdx = 0, permutationTimes = 1000, seedNum = 99999999){
  set.seed(seedNum)

  bin <- as.integer(bin)
  sp_order <- as.integer(sp_order)

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
    randColIdx <- sample(1:ncol(mat),replace = F)
    idx <- cbind(randRowIdx, randColIdx)
    randomExp1 <- mat[idx]

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
#' @param bin An integer specifying the number of bins. Default is 6.
#' @param sp_order An integer specifying the spline order. Must be less than `bin`. Default is 2.
#' @param bulkIdx Index to divide the task when processing many permutations. Default is 0.
#' @param permutationTimes Number of permutations for sampling. Default is 1000.
#' @param seedNum Seed for random number generation. Default is 99999999.
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
CMIBiCondimat2matPermu <- function(mat1, mat2, condi1, condi2, bin = 6, sp_order = 2, bulkIdx = 0, permutationTimes = 1000, seedNum = 99999999){
  set.seed(seedNum)

  # Validate inputs as non-empty numeric vectors
  bin <- as.integer(bin)
  sp_order <- as.integer(sp_order)

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
    permutationMIVector[p] <- getCMIBiCondi(randomExp1,randomExp2, condi1, condi2, bin = bin, sp_order = sp_order) /
      max(getMI(randomExp1,randomExp1, bin = bin, sp_order = sp_order),
          getMI(randomExp2,randomExp2, bin = bin, sp_order = sp_order))

  }
  permutationMIVector[permutationMIVector<0] = 0
  return(permutationMIVector)
}
