#include "RTM_MvNorm.h"
#include "gsl_mat_ext.h"

using namespace std;
using std::string;


////// //////
mvnorm::mvnorm()
{
  // no mean and covariance matrix specified - allocate no further memory
  mean = 0;
  log_cov_det = 0;
  covariance_inverse = 0;
}
mvnorm::mvnorm(const gsl_vector* in_mean, const gsl_matrix* in_covariance)
{
  if(in_mean->size != in_covariance->size1 || in_mean->size != in_covariance->size2)
    {
      printf("Mean vector and covariance matrices not of matching dimension\n");
      exit(2);
    }

  mean = gsl_vector_alloc(in_mean->size);
  gsl_vector_memcpy(mean, in_mean);

  // Is the specified covariance matrix symmetric?
  for(int int_i = 1; int_i < in_covariance->size1; int_i++)
    for(int int_j = 0; int_j < int_i; int_j++)
      {
	if(gsl_matrix_get(in_covariance, int_i, int_j) != gsl_matrix_get(in_covariance, int_j, int_i))
	  {
	    printf("Specified covariance matrix not symmetric");
	    exit(2);
	  }
      }

  // Invert the covariance matrix - and calculate the determinant en route
  covariance_inverse = gsl_matrix_alloc(in_covariance->size1, in_covariance->size2);
  gsl_matrix_memcpy(covariance_inverse, in_covariance);
  gsl_linalg_cholesky_decomp(covariance_inverse);
  // The product of diagonal entries of the cholesky decomposition is the square root of the determinant
  log_cov_det = 0.0;
  for(int int_i = 0; int_i < in_mean->size; int_i++)
    log_cov_det += log(gsl_matrix_get(covariance_inverse, int_i, int_i));
  log_cov_det *= 2;  // Account for the square root

  gsl_linalg_cholesky_invert(covariance_inverse);

}

#define DIM_CHECK_ERROR(func_name)       {		\
    string err_string("In vector to ");			\
    err_string += func_name;				\
    err_string += " not of correct dimension\n";	\
    perror(err_string.c_str());				\
    exit(2);						\
  }

double mvnorm::ld_mvnorm_nonnorm(const gsl_vector* in_x)
{

  if( mean->size != in_x->size)
    DIM_CHECK_ERROR("mvnorm::ld_mvnorm_nonnorm");

  gsl_vector* tempvec1 = gsl_vector_alloc(in_x->size);
  gsl_vector* tempvec2 = gsl_vector_alloc(in_x->size);

  gsl_vector_memcpy(tempvec1, in_x);
  gsl_vector_sub(tempvec1, mean);

  gsl_blas_dsymv(CblasLower, -0.5, covariance_inverse, tempvec1, 0, tempvec2);

  double outval;

  gsl_blas_ddot(tempvec1, tempvec2, &outval);

  gsl_vector_free(tempvec1);
  gsl_vector_free(tempvec2);

  return outval;

}

double mvnorm::ld_mvnorm(const gsl_vector* in_x)
{

  if( mean->size != in_x->size)
    DIM_CHECK_ERROR("mvnorm::ld_mvnorm");

  double outval = ld_mvnorm_nonnorm(in_x);

  outval -= (0.5) * log_cov_det;

  outval -= (in_x->size) * log(2 * M_PI) / 2.0;

  return outval;
}

double mvnorm::ld_mvnorm_ratio(const gsl_vector* in_x1, const gsl_vector* in_x2)
{

  if( (mean->size != in_x1->size) || (mean->size != in_x2->size) )
    DIM_CHECK_ERROR("mvnorm::ld_mvnorm_ratio");

  return (ld_mvnorm_nonnorm(in_x1) - ld_mvnorm_nonnorm(in_x2));
}

mvnorm::mvnorm(const mvnorm &in)
    : log_cov_det(in.log_cov_det),
      mean(nullptr),
      covariance_inverse(nullptr) {
    if (in.mean != nullptr) {
	gsl_vector_alloc(in.mean->size);
	gsl_vector_memcpy(mean, in.mean);
    }
    if (in.covariance_inverse != nullptr) {
	gsl_matrix_alloc(in.covariance_inverse->size1, in.covariance_inverse->size2);
	gsl_matrix_memcpy(covariance_inverse, in.covariance_inverse);
    }
}

mvnorm& mvnorm::operator=(const mvnorm &in) {
    if (this != &in) {
	log_cov_det = in.log_cov_det;
	if (in.mean != nullptr) {
	    gsl_vector_alloc(in.mean->size);
	    gsl_vector_memcpy(mean, in.mean);
	}
	if (in.covariance_inverse != nullptr) {
	    gsl_matrix_alloc(in.covariance_inverse->size1, in.covariance_inverse->size2);
	    gsl_matrix_memcpy(covariance_inverse, in.covariance_inverse);
	}
    }
    return *this;
}

mvnorm::~mvnorm()
{
  if(mean != 0)
    {
      gsl_vector_free(mean);
      gsl_matrix_free(covariance_inverse);
    }
}
