// info_utils.h

#ifndef INFO_UTILS_H
#define INFO_UTILS_H

#ifdef __cplusplus
extern "C" {
#endif

// Function to compute the univariate entropy.
double ComputeEntropy(const double* prob_from_b_spline, int num_samples, int num_bins);

// Function to compute the bivariate entropy.
double ComputeJointEntropy(const double* prob_from_b_spline_1, const double* prob_from_b_spline_2, int num_samples, int num_bins);

// Function to compute the trivariate entropy.
double ComputeTrivarEntropy(const double* prob_from_b_spline_1, const double* prob_from_b_spline_2, const double* prob_from_b_spline_3, int num_samples, int num_bins);

// Function to compute the quadrivariate entropy.
double ComputeQuadrivarEntropy(const double* prob_from_b_spline_1, const double* prob_from_b_spline_2, const double* prob_from_b_spline_3, const double* prob_from_b_spline_4, int num_samples, int num_bins);

#ifdef __cplusplus
}
#endif

#endif // INFO_UTILS_H
