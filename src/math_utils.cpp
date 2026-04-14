#include "math_utils.h"

#include <cmath>  // Include the cmath library for log function
#include <stdexcept>  // Include the library for std::invalid_argument
#include <utility>    // For std::pair
#include <algorithm>
#include <vector>

/**
 * @brief Computes the base-2 logarithm of a given number.
 * @param number The number for which the base-2 logarithm is to be calculated. Must be positive.
 * @return The base-2 logarithm of the given number.
 */
double getLogBase2(double number) {
    // Validate input - check if the number is greater than zero
    if (number <= 0.0) {
        throw std::invalid_argument("Input must be a positive number.");
    }
    return std::log(number) / std::log(2.0);
}


/**
 * @brief Finds both the minimum and maximum values in an array of doubles.
 *
 * @return A pair containing the minimum and maximum values.
 * @throws std::invalid_argument If the array is empty or vector_size is less than or equal to zero.
 */
std::pair<double, double> findMinMaxValues(const double* double_vector, int vector_size) {
    if (double_vector == nullptr || vector_size <= 0) {
        throw std::invalid_argument("Array must not be null and vector_size must be greater than zero.");
    }

    double MaxValue = double_vector[0];
    double MinValue = double_vector[0];

    for (int i = 1; i < vector_size; ++i) {
        if (double_vector[i] > MaxValue) {
            MaxValue = double_vector[i];
        }
        if (double_vector[i] < MinValue) {
            MinValue = double_vector[i];
        }
    }

    return std::make_pair(MinValue, MaxValue);
}

/**
 * @brief Normalizes an array of doubles to a range of 0 to 1 based on the input array's minimum and maximum values.
 *
 * @throws std::invalid_argument If the input array is empty or null.
 */
void normalizeData(const double* orig_vector, double* normalized_vector, int vector_size) {
    if (orig_vector == nullptr || normalized_vector == nullptr || vector_size <= 0) {
        throw std::invalid_argument("Input and output arrays must not be null, and vector_size must be greater than zero.");
    }

    // Find minimum and maximum values
    auto [xMin, xMax] = findMinMaxValues(orig_vector, vector_size);

    // Prevent division by zero in case of a flat vector
    if (xMax == xMin) {
        xMax = xMin + 1.0;
    }

    // Normalize each element
    for (int i = 0; i < vector_size; ++i) {
        normalized_vector[i] = (orig_vector[i] - xMin) / (xMax - xMin);
    }
}


/**
 * @brief Generates a knot vector for a B-spline.
 *
 * This function initializes a knot vector based on the specified number of bins and spline order.
 * The knot vector is composed of three segments:
 * 1. Leading knots, set to 0.0, for the spline order count.
 * 2. Internal knots, equally spaced between 0.0 and 1.0, determined by the number of internal points (num_bins - spline_order).
 * 3. Trailing knots, set to 1.0, for the spline order count.
 *
 * @param knots Array to store the generated knot vector. It should be preallocated with sufficient size to hold 2 * spline_order + num_bins elements.
 * @param num_bins Total number of bins for the B-spline.
 * @param spline_order The spline order (degree + 1) which determines the number of leading and trailing knots set to 0.0 and 1.0, respectively.
 */
void GenerateKnotVector(double* knots, int num_bins, int spline_order) {
  int n_internal_points = num_bins - spline_order;

  // Set leading spline order knots to 0
  for (int i = 0; i < spline_order; ++i) {
    knots[i] = 0.0;
  }

  // Set internal knots with calculated divisions
  for (int i = 0; i < n_internal_points; ++i) {
    knots[spline_order + i] = static_cast<double>(i + 1) / (n_internal_points + 1);
  }

  // Set trailing spline order knots to 1
  for (int i = spline_order + n_internal_points; i < 2 * spline_order + n_internal_points; ++i) {
    knots[i] = 1.0;
  }
}

/**
 * @brief Computes the B-spline basis function value.
 *
 * This function computes the basis function value for a given index, polynomial degree,
 * and parameter value using a specified knot vector.
 *
 * @param kVector The knot vector.
 * @return The basis function value.
 */
 double computeBSplineBasis(int knot_idx, int poly_degree, double number, const double* kVector, int numBins, double tolerance) {
    // Base case for recursion
    if (poly_degree == 1) {
        if ((number >= kVector[knot_idx] && number < kVector[knot_idx + 1] && kVector[knot_idx] < kVector[knot_idx + 1]) ||
            (std::fabs(number - kVector[knot_idx + 1]) < tolerance && (knot_idx + 1 == numBins))) {
            return 1.0;
        }
        return 0.0;
    }

    // Initialize coefficients
    double e1 = 0.0, e2 = 0.0;

    // Calculate denominators and numerators for the recursive formula
    double d1 = kVector[knot_idx + poly_degree - 1] - kVector[knot_idx];
    double n1 = number - kVector[knot_idx];
    double d2 = kVector[knot_idx + poly_degree] - kVector[knot_idx + 1];
    double n2 = kVector[knot_idx + poly_degree] - number;

    // Recursive step for non-zero denominators
    if (d1 > tolerance) {
        e1 = (n1 / d1) * computeBSplineBasis(knot_idx, poly_degree - 1, number, kVector, numBins, tolerance);
    }
    if (d2 > tolerance) {
        e2 = (n2 / d2) * computeBSplineBasis(knot_idx + 1, poly_degree - 1, number, kVector, numBins, tolerance);
    }

    // Ensure result is not negative due to floating-point rounding
    double result = e1 + e2;
    return std::max(0.0, result);
}

/**
 * @brief Computes the probability vector using B-spline basis functions.
 *
 * This function computes the probability vector used for entropy,
 * and mutual information estimation.
 *
 * @param orig_vector The data vector.
 * @param knots The knot vector.
 * @param prob_vector The probability vector.
 * @return The prob_vector for later calculation.
 */
void prob_from_b_spline(const double* orig_vector, const double* knots, double* prob_vector, int vector_size, int splineOrder, int numBins, double tolerance) {
    if (orig_vector == nullptr || knots == nullptr || prob_vector == nullptr || vector_size <= 0 || splineOrder <= 0 || numBins <= 0) {
        throw std::invalid_argument("Invalid input parameters for weight computation.");
    }

    std::vector<double> normalized_vector(vector_size);

    // Normalize the sample points with given range boundaries
    normalizeData(orig_vector, normalized_vector.data(), vector_size);

    // Compute prob_vector for each data and bin
    for (int i = 0; i < vector_size; ++i) {
        for (int j = 0; j < numBins; ++j) {
            prob_vector[j * vector_size + i] = computeBSplineBasis(j, splineOrder, normalized_vector[i], knots, numBins, tolerance);
        }
    }
}

