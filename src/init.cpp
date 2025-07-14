#include <R.h>
#include <Rinternals.h>

// Declare the function prototype
extern "C" {

  void compute_entropy_univar(const double* data_vector,
                              const int* n, const int* bin,
                              const int* sp_order, double* entropy_out);

  void compute_entropy_bivar(const double* data_vector_1,
                             const double* data_vector_2,
                             const int* n, const int* bin,
                             const int* sp_order, double* entropy_out);

  void compute_entropy_trivar(const double* data_vector_1,
                              const double* data_vector_2,
                              const double* data_vector_3,
                              const int* n, const int* bin,
                              const int* sp_order, double* entropy_out);

  void compute_entropy_quadrivar(const double* data_vector_1,
                                 const double* data_vector_2,
                                 const double* data_vector_3,
                                 const double* data_vector_4,
                                 const int* n, const int* bin, const int* sp_order, double* entropy_out);

  void compute_mutual_information(const double* data_vector_1,
                                  const double* data_vector_2,
                                  const int* n, const int* bin,
                                  const int* sp_order, double* mi_out);

  void compute_mutual_information_12_to_3(const double* data_vector_1,
                                          const double* data_vector_2,
                                          const double* data_vector_3,
                                          const int* n, const int* bin,
                                          const int* sp_order, double* mi_out);

  void compute_conditional_mutual_information_12_to_3(const double* data_vector_1,
                                                      const double* data_vector_2,
                                                      const double* data_vector_3,
                                                      const int* n, const int* bin,
                                                      const int* sp_order, double* mi_out);

  void compute_conditional_mutual_information_12_to_34(const double* data_vector_1,
                                                       const double* data_vector_2,
                                                       const double* data_vector_3,
                                                       const double* data_vector_4,
                                                       const int* n, const int* bin,
                                                       const int* sp_order, double* mi_out);

}

static const R_CMethodDef CEntries[] = {
  {"compute_entropy_univar", (DL_FUNC) &compute_entropy_univar, 5},
  {"compute_entropy_bivar", (DL_FUNC) &compute_entropy_bivar, 6},
  {"compute_entropy_trivar", (DL_FUNC) &compute_entropy_trivar, 7},
  {"compute_entropy_quadrivar", (DL_FUNC) &compute_entropy_quadrivar, 8},
  {"compute_mutual_information", (DL_FUNC) &compute_mutual_information, 6},
  {"compute_mutual_information_12_to_3", (DL_FUNC) &compute_mutual_information_12_to_3, 7},
  {"compute_conditional_mutual_information_12_to_3", (DL_FUNC) &compute_conditional_mutual_information_12_to_3, 7},
  {"compute_conditional_mutual_information_12_to_34", (DL_FUNC) &compute_conditional_mutual_information_12_to_34, 8},
  {NULL, NULL, 0}
};

extern "C" void R_init_COATS(DllInfo *dll) {
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
