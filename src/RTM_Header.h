#ifndef HEADER_header_
#define HEADER_header_

// LOADING IN LIBRARIES
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstddef>
#include <cmath>
#include <cstring>
#include <string>
#include <csignal>
#include <ctime>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <cstdio>
#include <limits>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <sys/time.h>
#include <iomanip>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_trig.h>

#include <omp.h>

// GLOBAL CONSTANTS
#define NUM_AGE_GROUPS (7)
#define DAYS_IN_YEAR 365.25

// USEFUL INLINE FUNCTIONS
#define FN_MAX(a, b) (((a) > (b)) ? (a) : (b))
#define FN_MIN(a, b) (((a) > (b)) ? (b) : (a))

// THE TOLERANCE USED FOR TRUNCATING THE GAMMA DELAY DISTRIBUTIONS
#define EPSILON_TOLERANCE 0.0001

#define cTRUE 1
#define cFALSE 0

// DEBUGGING_LEVELS
#define DEBUG_ALL 20
#define DEBUG_DETAIL 17
#define DEBUG_PROGRESS  15
#define DEBUG_WARNING 10
#define DEBUG_ERROR 5
#define DEBUG_OFF 0

#ifndef DEBUG_LEVEL
#define DEBUG_LEVEL DEBUG_DETAIL
#endif

#define DEBUG(level, x) if (DEBUG_LEVEL >= level) {std::cout << __FILE__ << ":" << __LINE__ << " " << x << std::endl;}

#endif
