#ifndef HEADER_Rlike_
#define HEADER_Rlike_

#include "string_fns.h"
#include "distributions.h"
#include "gsl_vec_ext.h"

#define MAX_ZERO_ONE_PROPOSAL_VARIANCE 0.1
#define MAX_UNBOUNDED_PROPOSAL_VARIANCE 500

// Some error messages
#define DEFAULT_NO_LIKELIHOOD printf("Unrecognised likelihood function\n"); \
  exit(2)
#define DEFAULT_NO_DISTRIBUTION   printf("Unrecognised distribution type\n"); \
    exit(2)
#define DEFAULT_NO_LINK printf("Unrecognised link function\n"); \
  exit(2)
#define MULTIVARIATE_DISTRIBUTION printf("Multivariate prior evaluated using univariate prior function\n"); \
  exit(2)

// Function prototypes
void R_gl_fac_sublevel(gsl_vector_int*, const int, const int, const int, const int);
double R_univariate_prior_log_density(const double&, const distribution_type&, const gsl_vector*);
double R_univariate_prior_log_density_nonnorm(const double&, const distribution_type&, const gsl_vector*);
double R_univariate_prior_log_density_ratio(const double&, const double&, const distribution_type&, const gsl_vector*);
int num_parameters_by_distribution(const distribution_type&);
double R_inverse_link_function(double, const link_function);
void R_generalised_linear_regression(gsl_vector*, const gsl_matrix*, const gsl_vector*, link_function);
int num_parameters_by_distribution(const distribution_type&, const int&);
double random_walk_proposal(double&, const double&, const distribution_type&, const double&, gsl_rng*, const double&, const double&);
double variance_inflation_factor(const double&, const distribution_type&, const gsl_vector*);
void R_by_sum_gsl_vector_mono_idx(gsl_vector*, gsl_vector*, const gsl_vector*, const vector<int>&);

#endif
