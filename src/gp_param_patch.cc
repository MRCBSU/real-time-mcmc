#include "gp_param_patch.h"
#include "gsl_mat_ext.h"
#include "RTM_FunctDefs.h"

void gp_param_fn(gsl_vector* dest, const gsl_vector* src)
{
  gsl_vector_set(dest, 0, gsl_vector_get(src, 0));
  gsl_vector_set(dest, 1, gsl_vector_get(src, 1));
  gsl_vector_set(dest, 2, gsl_vector_get(src, 0) * gsl_vector_get(src, 2));
  gsl_vector_set(dest, 3, gsl_vector_get(src, 1) * gsl_vector_get(src, 3));
}


/// LOCAL VERSION OF FUNCTION IN R_like_fns.cc
void R_generalised_linear_regression_GP_PATCH(gsl_vector* out_y, const gsl_matrix *design_X,
				     const gsl_vector* param_beta, link_function g)
{
 
  gsl_vector* temporary_vector = gsl_vector_alloc(out_y->size);

  // multiply the parameter vector by the design matrix
  gsl_matrix_multiplication_vector_dbl(temporary_vector, design_X, param_beta);

  // apply the inverse of the selected link function
  for(int int_i = 0; int_i < out_y->size; int_i++)
    gsl_vector_set(temporary_vector,
		   int_i,
		   R_inverse_link_function(gsl_vector_get(temporary_vector,
							  int_i),
					   g));

  gp_param_fn(out_y, temporary_vector);

  gsl_vector_free(temporary_vector);

}



/// LOCAL VERSION OF FUNCTION IN RTM_WithinRegion.cc
void regional_matrix_parameter_GP_PATCH(gsl_matrix* out_mat, const gsl_vector* param_value, const regression_def& map_to_regional, const int region_index, const int time_steps_per_day)
{
  /// IN THIS FUNCTION, THE CODE LINES INVOLVING THE subdesign MATRICES HAVE BEEN DONE IN A VERY
  /// SAFE MANNER, ONCE WORKING, THIS MAY BE SPEEDED UP - INVESTIGATE

  // A MORE GENERAL VERSION OF THE FUNCTION regional_vector_parameter WHICH ALLOWS FOR TEMPORAL BREAKPOINTS AS WELL AS
  // BREAKPOINTS OVER AGES.


  // if there's any temporal, regional or age variation then the design matrix should be non-zero
  if(map_to_regional.design_matrix == 0)
    {
      // no variation, only first component of parameter value will be used.
      gsl_matrix_set_all(out_mat,
			 gsl_vector_get(param_value, 0));
    }
  else
    {
      // the design matrix should be set.
      // get the number of intervals over time (and age) for each region
      int dim_r = (map_to_regional.region_breakpoints == 0) ? 1 : map_to_regional.region_breakpoints->size + 1;
      int dim_a = (map_to_regional.age_breakpoints == 0) ? 1 : map_to_regional.age_breakpoints->size + 1;
      int dim_t = (map_to_regional.time_breakpoints == 0) ? 1 : map_to_regional.time_breakpoints->size + 1;

      gsl_matrix* subdesign = gsl_matrix_alloc(dim_t * dim_a, map_to_regional.design_matrix->size2);
      gsl_vector* intermediate_vec = gsl_vector_alloc(subdesign->size1);

      select_design_matrix(subdesign, map_to_regional.design_matrix, dim_r == 1,
			   region_index * dim_t * dim_a, dim_t * dim_a);

      R_generalised_linear_regression_GP_PATCH(intermediate_vec, subdesign,
				      param_value, map_to_regional.regression_link);

      // now fill out the out_mat according to the breakpoints - first need to multiply up the breakpoints to the correct time scale
      gsl_vector_int* rescaled_temporal_breakpoints = 0;      
      if(dim_t > 1){
	rescaled_temporal_breakpoints = gsl_vector_int_alloc(map_to_regional.time_breakpoints->size);
	gsl_vector_int_memcpy(rescaled_temporal_breakpoints, map_to_regional.time_breakpoints);
	gsl_vector_int_scale(rescaled_temporal_breakpoints, time_steps_per_day);
      }
      mat_breakpoint_cut(out_mat, rescaled_temporal_breakpoints, map_to_regional.age_breakpoints, intermediate_vec);
  

      // free any memory allocated within this scope
      if(rescaled_temporal_breakpoints != 0)
	gsl_vector_int_free(rescaled_temporal_breakpoints);
      gsl_vector_free(intermediate_vec);
      gsl_matrix_free(subdesign);

    }

}
