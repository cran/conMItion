#' @title Calculate P-Value from a Sorted Distribution
#' @description Computes the P-value for a given numeric value `x` based on its position within a sorted distribution.
#' This function utilizes a binary search approach for efficient computation.
#'
#' @param x A numeric value for which the P-value is to be calculated.
#' @param sorted_Distri A numeric vector representing a sorted distribution. This distribution must be sorted in ascending order.
#'
#' @return A numeric value indicating the P-value, representing the proportion of values in `sorted_Distri`
#' that are greater than or equal to `x`.
#'
#' @examples
#' x <- rnorm(1)
#' sorted_dist <- sort(rnorm(100))
#' getPValue(x, sorted_dist)
#'
#' @export
getPValue <- function(x,sorted_Distri) {
  sizeD = length(sorted_Distri)
  if (x>sorted_Distri[sizeD]) {return(0)}
  if (x<=sorted_Distri[1]) {return(sizeD)}
  beginIdx = 1;endIdx = sizeD;
  while ((endIdx-beginIdx)!=1) {
    midIdx = ceiling((endIdx+beginIdx)/2)
    if(sorted_Distri[midIdx]<x) {beginIdx = midIdx}
    else {endIdx = midIdx}
  }
  if (sorted_Distri[endIdx]==x) {result = (sizeD-endIdx+1)}
  else (result = (sizeD-beginIdx))

  return(result/sizeD)
}
