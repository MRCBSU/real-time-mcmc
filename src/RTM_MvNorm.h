#ifndef HEADER_mvnorm_
#define HEADER_mvnorm_

#include "RTM_Header.h"

using namespace std;
using std::string;

// DEFINE A CLASS FOR HANDLING MULTIVARIATE NORMAL RANDOM VARIABLES

class mvnorm
{
private:
  gsl_vector* mean;
  double log_cov_det;
  gsl_matrix* covariance_inverse;
public:
  mvnorm();
  mvnorm(const gsl_vector*, const gsl_matrix*);
  ~mvnorm();

  mvnorm(const mvnorm& in);
  mvnorm& operator=(const mvnorm& in);
  
  double ld_mvnorm(const gsl_vector*);
  double ld_mvnorm_nonnorm(const gsl_vector*);
  double ld_mvnorm_ratio(const gsl_vector*, const gsl_vector*);
};


#endif
