#ifndef HEADER_gslmat_
#define HEADER_gslmat_

// FUNCTIONS IN gsl_matrix_exts.cc
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <string>

using namespace std;
using std::string;

void gsl_matrix_sscanf(const string, gsl_matrix *,
                       const char *chr_delim = ",;");
void gsl_matrix_int_sscanf(const string, gsl_matrix_int *,
                           const char *chr_delim = ",;");
void gsl_compute_max_evalue_evector2(double *, gsl_vector *, const gsl_matrix *,
                                     const bool evec_flag = true);
void gsl_matrix_multiplication_vector_dbl(gsl_vector *, const gsl_matrix *,
                                          const gsl_vector *);
void gsl_matrix_multiplication_vector_product_by_row_dbl(gsl_matrix *,
                                                         const gsl_vector *);
void gsl_matrix_realloc(gsl_matrix *&, const int &, const int &);

#endif
