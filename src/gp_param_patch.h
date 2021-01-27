#ifndef HEADER_gp_patch_
#define HEADER_gp_patch_

#include "RTM_Header.h"
#include "RTM_StructDefs.h"
#include "RTM_updParams.h"


void regional_matrix_parameter_GP_PATCH(gsl_matrix*, const gsl_vector*, const regression_def&, const int, const int);

void regional_matrix_parameter_GP_PATCH(gsl_matrix* out_mat, const updParamSet &paramSet, const upd::paramIndex index, const int region_index, const int time_steps_per_day);

#endif
