#include "I0_patch.h"

void I0_patch(gsl_vector* Iseed)
{
  gsl_vector_set_zero(Iseed);

  gsl_vector_set(Iseed, 3, 2);
  gsl_vector_set(Iseed, 5, 1);
  gsl_vector_set(Iseed, 16, 1);
  gsl_vector_set(Iseed, 25, 1);
  gsl_vector_set(Iseed, 35, 1);

  double mat_sum = gsl_vector_sum_elements(Iseed);
  gsl_vector_scale(Iseed, 1.0 / mat_sum);

}
