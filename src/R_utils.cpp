#include "info_utils.h"
#include "math_utils.h"
#include "R_utils.h"

#include <vector>
#include <memory> // For std::unique_ptr

extern "C" {
  /**
   * @brief Computes univariate entropy.
   * @param data_vector Data vector used for entropy computation.
   * @param n Size of the data.
   * @param bin The number of bins in entropy estimation.
   * @param sp_order The spline order. This value must be smaller than bin.
   *
   * @return double entropy_out.
   */
  void compute_entropy_univar(const double* data_vector, const int* n, const int* bin, const int* sp_order, double* entropy_out) {

    auto knot_vector = std::make_unique<double[]>(*bin + *sp_order);
    auto prob_vector = std::make_unique<double[]>(*bin * *n);

    GenerateKnotVector(knot_vector.get(), *bin, *sp_order);
    prob_from_b_spline(data_vector, knot_vector.get(), prob_vector.get(), *n, *sp_order, *bin);

    double entropy = ComputeEntropy(prob_vector.get(), *n, *bin);

    *entropy_out = entropy;
  }

  /**
   * @brief Computes bivariate entropy.
   * @param data_vector_1/2 Data vectors used for entropy computation.
   * @param n Size of the data.
   * @param bin The number of bins in entropy estimation.
   * @param sp_order The spline order. This value must be smaller than bin.
   *
   * @return double entropy_out.
   */
  void compute_entropy_bivar(const double* data_vector_1, const double* data_vector_2, const int* n, const int* bin, const int* sp_order, double* entropy_out) {
    // Use smart pointers to manage memory automatically
    auto knot_vector = std::make_unique<double[]>(*bin + *sp_order);
    auto prob_vector_1 = std::make_unique<double[]>(*bin * *n);
    auto prob_vector_2 = std::make_unique<double[]>(*bin * *n);

    GenerateKnotVector(knot_vector.get(), *bin, *sp_order);
    prob_from_b_spline(data_vector_1, knot_vector.get(), prob_vector_1.get(), *n, *sp_order, *bin);
    prob_from_b_spline(data_vector_2, knot_vector.get(), prob_vector_2.get(), *n, *sp_order, *bin);

    double entropy = ComputeJointEntropy(prob_vector_1.get(), prob_vector_2.get(), *n, *bin);

    *entropy_out = entropy;
  }

  /**
   * @brief Computes trivariate entropy.
   * @param data_vector_1/2/3 Data vectors used for entropy computation.
   * @param n Size of the data.
   * @param bin The number of bins in entropy estimation.
   * @param sp_order The spline order. This value must be smaller than bin.
   *
   * @return double entropy_out.
   */
  void compute_entropy_trivar(const double* data_vector_1, const double* data_vector_2, const double* data_vector_3, const int* n, const int* bin, const int* sp_order, double* entropy_out) {
    // Use smart pointers to manage memory automatically
    auto knot_vector = std::make_unique<double[]>(*bin + *sp_order);
    auto prob_vector_1 = std::make_unique<double[]>(*bin * *n);
    auto prob_vector_2 = std::make_unique<double[]>(*bin * *n);
    auto prob_vector_3 = std::make_unique<double[]>(*bin * *n);

    GenerateKnotVector(knot_vector.get(), *bin, *sp_order);
    prob_from_b_spline(data_vector_1, knot_vector.get(), prob_vector_1.get(), *n, *sp_order, *bin);
    prob_from_b_spline(data_vector_2, knot_vector.get(), prob_vector_2.get(), *n, *sp_order, *bin);
    prob_from_b_spline(data_vector_3, knot_vector.get(), prob_vector_3.get(), *n, *sp_order, *bin);

    double entropy = ComputeTrivarEntropy(prob_vector_1.get(), prob_vector_2.get(), prob_vector_3.get(), *n, *bin);

    *entropy_out = entropy;
  }

  /**
   * @brief Computes quadrivariate entropy.
   * @param data_vector_1/2/3/4 Data vectors used for entropy computation.
   * @param n Size of the data.
   * @param bin The number of bins in entropy estimation.
   * @param sp_order The spline order. This value must be smaller than bin.
   *
   * @return double entropy_out.
   */
  void compute_entropy_quadrivar(const double* data_vector_1, const double* data_vector_2, const double* data_vector_3, const double* data_vector_4, const int* n, const int* bin, const int* sp_order, double* entropy_out) {
    // Use smart pointers to manage memory automatically
    auto knot_vector = std::make_unique<double[]>(*bin + *sp_order);
    auto prob_vector_1 = std::make_unique<double[]>(*bin * *n);
    auto prob_vector_2 = std::make_unique<double[]>(*bin * *n);
    auto prob_vector_3 = std::make_unique<double[]>(*bin * *n);
    auto prob_vector_4 = std::make_unique<double[]>(*bin * *n);

    GenerateKnotVector(knot_vector.get(), *bin, *sp_order);
    prob_from_b_spline(data_vector_1, knot_vector.get(), prob_vector_1.get(), *n, *sp_order, *bin);
    prob_from_b_spline(data_vector_2, knot_vector.get(), prob_vector_2.get(), *n, *sp_order, *bin);
    prob_from_b_spline(data_vector_3, knot_vector.get(), prob_vector_3.get(), *n, *sp_order, *bin);
    prob_from_b_spline(data_vector_4, knot_vector.get(), prob_vector_4.get(), *n, *sp_order, *bin);

    double entropy = ComputeQuadrivarEntropy(prob_vector_1.get(), prob_vector_2.get(), prob_vector_3.get(),prob_vector_4.get(), *n, *bin);

    *entropy_out = entropy;
  }

  /**
   * @brief Computes mutual information between one variable and another variable.
   * @param data_vector_1/2 Data vectors used for MI computation.
   * @param n Size of the data.
   * @param bin The number of bins in entropy estimation.
   * @param sp_order The spline order. This value must be smaller than bin.
   *
   * @return double mi_out.
   */
  void compute_mutual_information(const double* data_vector_1, const double* data_vector_2, const int* n, const int* bin, const int* sp_order, double* mi_out) {
    // Use smart pointers for automatic memory management
    auto knot_vector = std::make_unique<double[]>(*bin + *sp_order);
    auto prob_vector_1 = std::make_unique<double[]>(*bin * *n);
    auto prob_vector_2 = std::make_unique<double[]>(*bin * *n);

    GenerateKnotVector(knot_vector.get(), *bin, *sp_order);
    prob_from_b_spline(data_vector_1, knot_vector.get(), prob_vector_1.get(), *n, *sp_order, *bin);
    prob_from_b_spline(data_vector_2, knot_vector.get(), prob_vector_2.get(), *n, *sp_order, *bin);

    double entropy_1 = ComputeEntropy(prob_vector_1.get(), *n, *bin);
    double entropy_2 = ComputeEntropy(prob_vector_2.get(), *n, *bin);
    double joint_entropy_12 = ComputeJointEntropy(prob_vector_1.get(), prob_vector_2.get(), *n, *bin);
    double mi = entropy_1 + entropy_2 - joint_entropy_12;

    *mi_out = mi;
  }

  /**
   * @brief Computes mutual information between two variables and the third variable.
   * @param data_vector_1/2/3 Data vectors used for MI computation.
   * @param n Size of the data.
   * @param bin The number of bins in entropy estimation.
   * @param sp_order The spline order. This value must be smaller than bin.
   *
   * @return double mi_out.
   */
  void compute_mutual_information_12_to_3(const double* data_vector_1, const double* data_vector_2, const double* data_vector_3, const int* n, const int* bin, const int* sp_order, double* mi_out) {
    auto knot_vector = std::make_unique<double[]>(*bin + *sp_order);
    auto prob_vector_1 = std::make_unique<double[]>(*bin * *n);
    auto prob_vector_2 = std::make_unique<double[]>(*bin * *n);
    auto prob_vector_3 = std::make_unique<double[]>(*bin * *n);

    GenerateKnotVector(knot_vector.get(), *bin, *sp_order);
    prob_from_b_spline(data_vector_1, knot_vector.get(), prob_vector_1.get(), *n, *sp_order, *bin);
    prob_from_b_spline(data_vector_2, knot_vector.get(), prob_vector_2.get(), *n, *sp_order, *bin);
    prob_from_b_spline(data_vector_3, knot_vector.get(), prob_vector_3.get(), *n, *sp_order, *bin);

    double entropy_12 = ComputeJointEntropy(prob_vector_1.get(),prob_vector_2.get(), *n, *bin);
    double entropy_3  = ComputeEntropy(prob_vector_3.get(), *n, *bin);
    double entropy_123 = ComputeTrivarEntropy(prob_vector_1.get(),prob_vector_2.get(),prob_vector_3.get(), *n, *bin);

    double mi = entropy_12 + entropy_3 - entropy_123;

    *mi_out = mi;
  }

  /**
   * @brief Computes conditional mutual information between two variables given the third variable.
   * @param data_vector_1/2/3 Data vectors used for CMI computation.
   * @param n Size of the data.
   * @param bin The number of bins in entropy estimation.
   * @param sp_order The spline order. This value must be smaller than bin.
   *
   * @return double mi_out.
   */
  void compute_conditional_mutual_information_12_to_3(const double* data_vector_1, const double* data_vector_2, const double* data_vector_3, const int* n, const int* bin, const int* sp_order, double* mi_out) {
    auto knot_vector = std::make_unique<double[]>(*bin + *sp_order);
    auto prob_vector_1 = std::make_unique<double[]>(*bin * *n);
    auto prob_vector_2 = std::make_unique<double[]>(*bin * *n);
    auto prob_vector_3 = std::make_unique<double[]>(*bin * *n);

    GenerateKnotVector(knot_vector.get(), *bin, *sp_order);
    prob_from_b_spline(data_vector_1, knot_vector.get(), prob_vector_1.get(), *n, *sp_order, *bin);
    prob_from_b_spline(data_vector_2, knot_vector.get(), prob_vector_2.get(), *n, *sp_order, *bin);
    prob_from_b_spline(data_vector_3, knot_vector.get(), prob_vector_3.get(), *n, *sp_order, *bin);

    double entropy_13 = ComputeJointEntropy(prob_vector_1.get(),prob_vector_3.get(), *n, *bin);
    double entropy_23 = ComputeJointEntropy(prob_vector_2.get(),prob_vector_3.get(), *n, *bin);
    double entropy_3  = ComputeEntropy(prob_vector_3.get(), *n, *bin);
    double entropy_123 = ComputeTrivarEntropy(prob_vector_1.get(),prob_vector_2.get(),prob_vector_3.get(), *n, *bin);

    double mi = entropy_13 + entropy_23 - entropy_3 - entropy_123;

    *mi_out = mi;
  }

  /**
   * @brief Computes conditional mutual information between two variables given the third and fourth variables.
   * @param data_vector_1/2/3/4 Data vectors used for CMI computation.
   * @param n Size of the data.
   * @param bin The number of bins in entropy estimation.
   * @param sp_order The spline order. This value must be smaller than bin.
   *
   * @return double mi_out.
   */
  void compute_conditional_mutual_information_12_to_34(const double* data_vector_1, const double* data_vector_2, const double* data_vector_3, const double* data_vector_4, const int* n, const int* bin, const int* sp_order, double* mi_out) {
    auto knot_vector = std::make_unique<double[]>(*bin + *sp_order);
    auto prob_vector_1 = std::make_unique<double[]>(*bin * *n);
    auto prob_vector_2 = std::make_unique<double[]>(*bin * *n);
    auto prob_vector_3 = std::make_unique<double[]>(*bin * *n);
    auto prob_vector_4 = std::make_unique<double[]>(*bin * *n);

    GenerateKnotVector(knot_vector.get(), *bin, *sp_order);
    prob_from_b_spline(data_vector_1, knot_vector.get(), prob_vector_1.get(), *n, *sp_order, *bin);
    prob_from_b_spline(data_vector_2, knot_vector.get(), prob_vector_2.get(), *n, *sp_order, *bin);
    prob_from_b_spline(data_vector_3, knot_vector.get(), prob_vector_3.get(), *n, *sp_order, *bin);
    prob_from_b_spline(data_vector_4, knot_vector.get(), prob_vector_4.get(), *n, *sp_order, *bin);

    double entropy_134 = ComputeTrivarEntropy(prob_vector_1.get(),prob_vector_3.get(),prob_vector_4.get(), *n, *bin);
    double entropy_234 = ComputeTrivarEntropy(prob_vector_2.get(),prob_vector_3.get(),prob_vector_4.get(), *n, *bin);
    double entropy_34 = ComputeJointEntropy(prob_vector_3.get(),prob_vector_4.get(), *n, *bin);
    double entropy_1234 = ComputeQuadrivarEntropy(prob_vector_1.get(),prob_vector_2.get(),prob_vector_3.get(),prob_vector_4.get(), *n, *bin);

    double mi = entropy_134 + entropy_234 - entropy_34 - entropy_1234;

    *mi_out = mi;
  }
}
