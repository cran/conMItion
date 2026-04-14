#include "math_utils.h"
#include "info_utils.h"

/**
 * @brief Computes univariate entropy.
 * @param prob_from_b_spline Probability vector from the data.
 * @param num_samples Size of the data.
 * @param num_bins The number of bins in entropy estimation.
 * @return double entropy.
 */
double ComputeEntropy(const double* prob_from_b_spline, int num_samples, int num_bins) {
  double entropy = 0.0;

  for (int bin_idx = 0; bin_idx < num_bins; ++bin_idx) {
    double bin_probability_sum = 0.0;

    for (int sample_idx = 0; sample_idx < num_samples; ++sample_idx) {
      bin_probability_sum += prob_from_b_spline[bin_idx * num_samples + sample_idx];
    }

    double average_probability = bin_probability_sum / num_samples;

    if (average_probability > 0) {
      entropy -= average_probability * getLogBase2(average_probability);
    }
  }

  return entropy;
}

/**
 * @brief Computes bivariate joint entropy.
 * @param prob_from_b_spline_1/2 Probability vector from the data.
 * @param num_samples Size of the data.
 * @param num_bins The number of bins in entropy estimation.
 * @return double joint_entropy.
 */
double ComputeJointEntropy(const double* prob_from_b_spline_1, const double* prob_from_b_spline_2, int num_samples, int num_bins) {
  double joint_entropy = 0.0;

  for (int bin_x = 0; bin_x < num_bins; ++bin_x) {
    for (int bin_y = 0; bin_y < num_bins; ++bin_y) {
      double joint_probability_sum = 0.0;

      for (int sample = 0; sample < num_samples; ++sample) {
        joint_probability_sum += prob_from_b_spline_1[bin_x * num_samples + sample] * prob_from_b_spline_2[bin_y * num_samples + sample];
      }

      double average_joint_probability = joint_probability_sum / num_samples;

      if (average_joint_probability > 0) {
        joint_entropy -= average_joint_probability * getLogBase2(average_joint_probability);
      }
    }
  }

  return joint_entropy;
}

/**
 * @brief Computes trivariate joint entropy.
 * @param prob_from_b_spline_1/2/3 Probability vector from the data.
 * @param num_samples Size of the data.
 * @param num_bins The number of bins in entropy estimation.
 * @return double trivar_entropy.
 */
double ComputeTrivarEntropy(const double* prob_from_b_spline_1, const double* prob_from_b_spline_2, const double* prob_from_b_spline_3, int num_samples, int num_bins) {
  double trivar_entropy = 0.0;

  for (int bin_x = 0; bin_x < num_bins; ++bin_x) {
    for (int bin_y = 0; bin_y < num_bins; ++bin_y) {
      for (int bin_z = 0; bin_z < num_bins; ++bin_z) {
        double trivar_probability_sum = 0.0;

        for (int sample = 0; sample < num_samples; ++sample) {
          trivar_probability_sum += prob_from_b_spline_1[bin_x * num_samples + sample] * prob_from_b_spline_2[bin_y * num_samples + sample] * prob_from_b_spline_3[bin_z * num_samples + sample];
        }

        double average_joint_probability = trivar_probability_sum / num_samples;

        if (average_joint_probability > 0) {
          trivar_entropy -= average_joint_probability * getLogBase2(average_joint_probability);
        }
      }
    }
  }

  return trivar_entropy;
}

/**
 * @brief Computes quadrivariate joint entropy.
 * @param prob_from_b_spline_1/2/3/4 Probability vector from the data.
 * @param num_samples Size of the data.
 * @param num_bins The number of bins in entropy estimation.
 * @return double quadrivar_entropy.
 */
double ComputeQuadrivarEntropy(const double* prob_from_b_spline_1, const double* prob_from_b_spline_2, const double* prob_from_b_spline_3, const double* prob_from_b_spline_4, int num_samples, int num_bins) {
  double quadrivar_entropy = 0.0;

  for (int bin_1 = 0; bin_1 < num_bins; ++bin_1) {
    for (int bin_2 = 0; bin_2 < num_bins; ++bin_2) {
      for (int bin_3 = 0; bin_3 < num_bins; ++bin_3) {
        for (int bin_4 = 0; bin_4 < num_bins; ++bin_4) {
          double quadrivar_probability_sum = 0.0;

          for (int sample = 0; sample < num_samples; ++sample) {
            quadrivar_probability_sum += prob_from_b_spline_1[bin_1 * num_samples + sample] * prob_from_b_spline_2[bin_2 * num_samples + sample] * prob_from_b_spline_3[bin_3 * num_samples + sample] * prob_from_b_spline_4[bin_4 * num_samples + sample];
          }

          double average_joint_probability = quadrivar_probability_sum / num_samples;

          if (average_joint_probability > 0) {
            quadrivar_entropy -= average_joint_probability * getLogBase2(average_joint_probability);
          }
        }
      }
    }
  }

  return quadrivar_entropy;
}
