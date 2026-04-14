#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <stdexcept>
#include <utility>

// Function to compute the base-2 logarithm of the given number
double getLogBase2(double number);

// Function to find the pair of minimum and maximum values in the vector
std::pair<double, double> findMinMaxValues(const double* double_vector, int vector_size);

// Function to normalize a data vector to the range [0, 1]
void normalizeData(const double* orig_vector, double* normalized_vector, int vector_size);

// Function to generate a knot vector
void GenerateKnotVector(double* knots, int num_bins, int spline_order);

// Function to compute the B-spline basis function value
double computeBSplineBasis(int knot_idx, int poly_degree, double number, const double* kVector, int numBins, double tolerance);

// Function to combine two probability vectors
void combin_two_prob(const double* vector_x, const double* vector_y, double* prob_vector, int vector_size, int numBins);

// Function to compute probabilities from the B-spline basis
void prob_from_b_spline(const double* orig_vector, const double* knots, double* prob_vector, int vector_size, int splineOrder, int numBins, double tolerance = 1e-9);

#endif // MATH_UTILS_H
