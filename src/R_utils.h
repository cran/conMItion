#ifndef R_UTILS_H
#define R_UTILS_H

#ifdef __cplusplus
extern "C" {
#endif

  // Function to compute the univariate entropy. This function is directly used by R.
  void compute_entropy_univar(const double* data_vector, const int* n, const int* bin, const int* sp_order, double* entropy_out);

  // Function to compute the bivariate entropy. This function is directly used by R.
  void compute_entropy_bivar(const double* data_vector_1, const double* data_vector_2, const int* n, const int* bin, const int* sp_order, double* entropy_out);

  // Function to compute the trivariate entropy. This function is directly used by R.
  void compute_entropy_trivar(const double* data_vector_1, const double* data_vector_2, const double* data_vector_3, const int* n, const int* bin, const int* sp_order, double* entropy_out);

  // Function to compute the quadrivariate entropy. This function is directly used by R.
  void compute_entropy_quadrivar(const double* data_vector_1, const double* data_vector_2, const double* data_vector_3, const double* data_vector_4, const int* n, const int* bin, const int* sp_order, double* entropy_out);

  // Function to compute the mutual information between two variables. This function is directly used by R.
  void compute_mutual_information(const double* data_vector_1, const double* data_vector_2, const int* n, const int* bin, const int* sp_order, double* mi_out);

  // Function to compute the mutual information between two variables and the third one. This function is directly used by R.
  void compute_mutual_information_12_to_3(const double* data_vector_1, const double* data_vector_2, const double* data_vector_3, const int* n, const int* bin, const int* sp_order, double* mi_out);

  // Function to compute the conditional mutual information between two variables given the third one. This function is directly used by R.
  void compute_conditional_mutual_information_12_to_3(const double* data_vector_1, const double* data_vector_2, const double* data_vector_3, const int* n, const int* bin, const int* sp_order, double* mi_out);

  // Function to compute the conditional mutual information between two variables given the third and fourth variables. This function is directly used by R.
  void compute_conditional_mutual_information_12_to_34(const double* data_vector_1, const double* data_vector_2, const double* data_vector_3, const double* data_vector_4, const int* n, const int* bin, const int* sp_order, double* mi_out);

#ifdef __cplusplus
}
#endif

#endif // R_UTILS_H
