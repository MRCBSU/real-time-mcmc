// FUNCTIONS IN gsl_vector_exts.cc

#ifndef HEADER_gslvec_
#define HEADER_gslvec_

#include <string>
#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>

using namespace std;
using std::string;

void gsl_vector_int_set_1ton(gsl_vector_int*, int start = 0);
double gsl_vector_sum_elements(register const gsl_vector*);
void gsl_vector_subvector_memcpy(gsl_vector*, gsl_vector*, const size_t);
void gsl_vector_sscanf(const string, gsl_vector*, const char* chr_delim = ",;");
void gsl_vector_int_sscanf(const string, gsl_vector_int*, const char* chr_delim = ",;");
void gsl_vector_int_set_seq(gsl_vector_int*, const int i_from = 1, const int i_by = 1);
void gsl_vector_set_seq(gsl_vector*, const double dbl_from = 1.0, const double dbl_by = 1.0);
double gsl_double_product_of_vector_elements(const gsl_vector*);
bool in_gsl_vector(gsl_vector*, const double);
bool in_gsl_vector_int(gsl_vector_int*, const int);
void gsl_populate_vector(vector<double>&, const gsl_vector*);
gsl_complex gsl_absmax_vector_complex(const gsl_vector_complex*);

#endif
