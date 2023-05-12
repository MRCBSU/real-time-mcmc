#include "RTM_StructDefs.h"
#include "RTM_FunctDefs.h"
#include "RTM_StructAssign.h"
#include "RTM_flagclass.h"
#include "gsl_mat_ext.h"
#include "gp_param_patch.h"

#define MIN_DELTA_T_TO_LENGTH_OF_STAY_RATIO 4.0

#define ERROR_INPUT_EXIT(err_string, calling_function) printf(err_string, calling_function); \
  exit(2);

#define ERROR_PARAM_INPUT(param_name) printf("Input error: redundant heterogeneity specified for parameter %s\n", param_name);


void vec_breakpoint_cut(gsl_vector* dest, const gsl_vector_int* vec_breakpoints, const gsl_vector* src)
{
  gsl_vector_view dest_sub;
  int offset = 0, sub_length = (vec_breakpoints == 0) ? dest->size : gsl_vector_int_get(vec_breakpoints, 0);
  int dim = (vec_breakpoints == 0) ? 1 : vec_breakpoints->size + 1;

  for(int inti = 0; inti < src->size; )
    {
      dest_sub = gsl_vector_subvector(dest, offset, sub_length);
      gsl_vector_set_all(&dest_sub.vector, gsl_vector_get(src, inti++));
      offset += sub_length;
      sub_length = ((inti < (dim - 1)) ? gsl_vector_int_get(vec_breakpoints, inti) : dest->size) - offset;
    }

}

void mat_breakpoint_cut(gsl_matrix* dest, const gsl_vector_int* vec_breakpoints_t, const gsl_vector_int* vec_breakpoints_a, const gsl_vector* src)
{ 
  gsl_matrix_view dest_sub;
  int offseti = 0, offsetj = 0;
  int sub_lengthi = (vec_breakpoints_t == 0) ? dest->size1 : gsl_vector_int_get(vec_breakpoints_t, 0);
  int sub_lengthj = (vec_breakpoints_a == 0) ? dest->size2 : gsl_vector_int_get(vec_breakpoints_a, 0);
  int sub_lengthj_init = sub_lengthj;

  int dim_1 = (vec_breakpoints_t == 0) ? 1 : vec_breakpoints_t->size + 1;
  int dim_2 = (vec_breakpoints_a == 0) ? 1 : vec_breakpoints_a->size + 1;

  int param_index = 0;

  for(int inti = 0; inti < dim_1; inti++)
    {

      for(int intj = 0; intj < dim_2; intj++)
	{

	  dest_sub = gsl_matrix_submatrix(dest, offseti, offsetj, sub_lengthi, sub_lengthj);
	  gsl_matrix_set_all(&dest_sub.matrix, gsl_vector_get(src, param_index++));

	  if(intj < (dim_2 - 1))
	    {
	      offsetj += sub_lengthj;
	      sub_lengthj = ((intj < (dim_2 - 2)) ? gsl_vector_int_get(vec_breakpoints_a, intj + 1) : dest->size2) - offsetj;
	    }
	  else if((intj == (dim_2 - 1)) && !(inti == (dim_1 - 1)))
	    {
	      offsetj = 0;
	      sub_lengthj = sub_lengthj_init;
	      offseti += sub_lengthi;
	      sub_lengthi = ((inti < (dim_1 - 2)) ? gsl_vector_int_get(vec_breakpoints_t, inti + 1) : dest->size1) - offseti;
	    }
	}

    }

}

void select_design_matrix(gslMatrix& dest, const gslMatrix& src, bool whole_matrix_flag,
			  int row_offset, int submatrix_nrows)
{

  // CCS TODO: Can't this be done by passing the view to the multiplication
  // fn rather than doing the memcpy?
    
  // get the sub-matrix relevant to this region
  if(whole_matrix_flag){ // use the entire design matrix
    //gsl_matrix_memcpy(dest, src);
    dest = src;
  } else { // no arbitrary splitting of the regions allowed, either full or zero heterogeneity.
    gsl_matrix_view subdesign_view = gsl_matrix_submatrix(*src,
							  row_offset,
							  0,
							  submatrix_nrows,
							  src.ncols());
    dest = &subdesign_view.matrix; // check after execution of this line if the code can be made more efficient - is subdesign necessary?
  }

}

void fixed_quantity(double& out_val, const double param_value, const string& param_name, const gslMatrix design)
{
  // TRIVIAL ERROR CHECK FOR PARAMETERS THAT ARE A FIXED QUANTITY OVER THE COURSE OF THE PANDEMIC
  if(design.nrows() > 0)
    {
      ERROR_PARAM_INPUT(param_name.c_str());
    }
  out_val = param_value;
}

void regional_scalar_parameter(double& out_val, const gsl_vector* param_value, const regression_def& map_to_regional, const int region_index)
{
  /// IN THIS FUNCTION, THE CODE LINES INVOLVING THE subdesign MATRICES HAVE BEEN DONE IN A VERY
  /// SAFE MANNER, THIS MAY BE SPEEDED UP - INVESTIGATE

  /// IF THE REGIONAL PARAMETER IS A SCALAR THEN IT REPRESENTS A QUANTITY FIXED OVER TIME AND AGE.
  /// THEREFORE ANY TEMPORAL OR AGE VARIATION SPECIFIED IN THE regression_def COMPONENT OF in_ump WILL
  /// THEREFORE BE IGNRORED

  // if there's any regional variation, then the design matrix should be non-zero
  if(map_to_regional.design_matrix.nrows() == 0)
    {
      // no regional variation, only first component of parameter value will be used
      out_val = gsl_vector_get(param_value, 0);
    }
  else
    {
      // the design matrix should be set.
      // get the number of breakpoints for each region
      int dim_r = (map_to_regional.region_breakpoints.size() == 0) ? 1 : map_to_regional.region_breakpoints.size() + 1;
      gslMatrix subdesign(1, map_to_regional.design_matrix.ncols());
      gsl_vector* intermediate_vec = gsl_vector_alloc(1);

      select_design_matrix(subdesign, map_to_regional.design_matrix, dim_r == 1, region_index, 1);

      R_generalised_linear_regression(intermediate_vec, *subdesign, param_value, map_to_regional.regression_link);

      out_val = gsl_vector_get(intermediate_vec, 0);

      // free any memory allocated within this scope
      gsl_vector_free(intermediate_vec);
      //gsl_matrix_free(subdesign);

    }
}

// Block version
void regional_scalar_parameter(double& out_val, const updParamSet &paramSet, const upd::paramIndex index, const int region_index) {
  // Design matrix is 1x1
  const updRegrDef &map = paramSet[index].map_to_regional;
  double x = map.design_matrix[region_index][0][0];
  double beta = paramSet.lookup0(index, region_index);

  out_val = R_inverse_link_function(x * beta, map.regression_link);
}

void regional_time_vector_parameter(gsl_vector* out_vec, const gsl_vector* param_value, const regression_def& map_to_regional, const int region_index, const int time_steps_per_day)
{
  // IN THIS FUNCTION, THE CODE LINES INVOLVING THE subdesign MATRICES HAVE BEEN DONE IN A VERY
  // SAFE MANNER, ONCE WORKING, THIS MAY BE SPEEDED UP - INVESTIGATE

  // PARAMETERS USING THIS FUNCTION SHOULD HAVE NO NEED OF VARIATION BY AGE VARIATION SPECIFIED
  // BY AGE IN THE regression_def COMPONENT OF in_ump WILL THEREFORE BE IGNORED.


  // if there's any temporal, regional or age variation then the design matrix should be non-zero
  if(map_to_regional.design_matrix.nrows() == 0)
    {
      // no variation, only first component of parameter value will be used.
      gsl_vector_set_all(out_vec,
			 gsl_vector_get(param_value, 0));
    }
  else
    {
      // the design matrix should be set.
      // get the number of breakpoints for each region
      int dim_r = (map_to_regional.region_breakpoints.size() == 0) ? 1 : map_to_regional.region_breakpoints.size() + 1;
      int dim_t = (map_to_regional.time_breakpoints.size() == 0) ? 1 : map_to_regional.time_breakpoints.size() + 1;

      gslMatrix subdesign(dim_t, map_to_regional.design_matrix.ncols());
      gsl_vector* intermediate_vec = gsl_vector_alloc(subdesign.nrows());

      select_design_matrix(subdesign, map_to_regional.design_matrix, dim_r == 1,
			   region_index * dim_t, dim_t);

      R_generalised_linear_regression(intermediate_vec, *subdesign,
				      param_value, map_to_regional.regression_link);

      // now fill out the out_vec according to the breakpoints
      //gsl_vector_int* rescaled_temporal_breakpoints = 0;
      gslVectorInt rescaled_temporal_breakpoints;
      if(dim_t > 1){
	//rescaled_temporal_breakpoints.alloc(map_to_regional.time_breakpoints.size());
	rescaled_temporal_breakpoints = map_to_regional.time_breakpoints;
	//gsl_vector_int_scale(rescaled_temporal_breakpoints, time_steps_per_day);
	rescaled_temporal_breakpoints *= time_steps_per_day;
      }
      vec_breakpoint_cut(out_vec, *rescaled_temporal_breakpoints, intermediate_vec);

      // free any memory allocated within this scope
      //if(rescaled_temporal_breakpoints->size != 0)
      //gsl_vector_int_free(rescaled_temporal_breakpoints);
      gsl_vector_free(intermediate_vec);
      //gsl_matrix_free(subdesign);

    }

}

//Design matrix already divided over regions, no need to divide further (?)
void regional_time_vector_parameter(gsl_vector* out_vec, const updParamSet &paramSet, upd::paramIndex index, const int region_index, const int time_steps_per_day) {
  
  gsl_vector_const_view view = paramSet.lookup(index, region_index);

  const updRegrDef &map = paramSet[index].map_to_regional;

  if(map.design_matrix[region_index].nrows() == 0)
  {
    // no variation, only first component of parameter value will be used.
    gsl_vector_set_all(out_vec, paramSet.lookup0(index, region_index));
  } else {
    int dim_t = (map.time_breakpoints.size() == 0) ? 1 : map.time_breakpoints.size() + 1;
    
    gsl_vector* intermediate_vec = gsl_vector_alloc(map.design_matrix[region_index].nrows());
    R_generalised_linear_regression(intermediate_vec, *map.design_matrix[region_index],
				    &view.vector, map.regression_link);
    
    gslVectorInt rescaled_temporal_breakpoints;
    if (dim_t > 1){
      rescaled_temporal_breakpoints = map.time_breakpoints;
      //gsl_vector_int_scale(rescaled_temporal_breakpoints, time_steps_per_day);
	rescaled_temporal_breakpoints *= time_steps_per_day;
      }
    vec_breakpoint_cut(out_vec, *rescaled_temporal_breakpoints, intermediate_vec);
    gsl_vector_free(intermediate_vec);
  }
}


// Block version
void regional_vector_parameter(gsl_vector* out_vec, const gsl_vector_const_view param_value, const updRegrDef& map_to_regional, const int region_index)
{
  /// IN THIS FUNCTION, THE CODE LINES INVOLVING THE subdesign MATRICES HAVE BEEN DONE IN A VERY
  /// SAFE MANNER, ONCE WORKING, THIS MAY BE SPEEDED UP - INVESTIGATE

  // IF THE REGIONAL PARAMETER IS A VECTOR THEN IT REPRESENTS AN INITIAL CONDITION.
  // THEREFORE WE DON'T NEED TO CONSIDER TEMPORAL VARIATION. ANY TEMPORAL VARIATION SPECIFIED
  // IN THE regression_def COMPONENT OF in_ump WILL THEREFORE BE IGNORED.


  // if there's any temporal, regional or age variation then the design matrix should be non-zero
  if(map_to_regional.design_matrix[region_index].nrows() == 0)
    {
      // no variation, only first component of parameter value will be used.
      gsl_vector_set_all(out_vec,
			 gsl_vector_get(&param_value.vector, 0));
    }
  else
    {
      // the design matrix should be set.
      // get the number of breakpoints for each region
      int dim_r = (map_to_regional.region_breakpoints.size() == 0) ? 1 : map_to_regional.region_breakpoints.size() + 1;
      int dim_a = (map_to_regional.age_breakpoints.size() == 0) ? 1 : map_to_regional.age_breakpoints.size() + 1;

      gslMatrix subdesign(dim_a, map_to_regional.design_matrix[region_index].ncols());
      gsl_vector* intermediate_vec = gsl_vector_alloc(subdesign.nrows());

      select_design_matrix(subdesign, map_to_regional.design_matrix[region_index], dim_r == 1,
			   region_index * dim_a, dim_a);

      R_generalised_linear_regression(intermediate_vec, *subdesign,
				       &param_value.vector, map_to_regional.regression_link);

      // now fill out the out_vec according to the breakpoints
      vec_breakpoint_cut(out_vec, *map_to_regional.age_breakpoints, intermediate_vec);

      // free any memory allocated within this scope
      gsl_vector_free(intermediate_vec);
      //gsl_matrix_free(subdesign);

    }

}


void regional_vector_parameter(gsl_vector* out_vec, const gsl_vector* param_value, const regression_def& map_to_regional, const int region_index)
{
  /// IN THIS FUNCTION, THE CODE LINES INVOLVING THE subdesign MATRICES HAVE BEEN DONE IN A VERY
  /// SAFE MANNER, ONCE WORKING, THIS MAY BE SPEEDED UP - INVESTIGATE

  // IF THE REGIONAL PARAMETER IS A VECTOR THEN IT REPRESENTS AN INITIAL CONDITION.
  // THEREFORE WE DON'T NEED TO CONSIDER TEMPORAL VARIATION. ANY TEMPORAL VARIATION SPECIFIED
  // IN THE regression_def COMPONENT OF in_ump WILL THEREFORE BE IGNORED.


  // if there's any temporal, regional or age variation then the design matrix should be non-zero
  if(map_to_regional.design_matrix.nrows() == 0)
    {
      // no variation, only first component of parameter value will be used.
      gsl_vector_set_all(out_vec,
			 gsl_vector_get(param_value, 0));
    }
  else
    {
      // the design matrix should be set.
      // get the number of breakpoints for each region
      int dim_r = (map_to_regional.region_breakpoints.size() == 0) ? 1 : map_to_regional.region_breakpoints.size() + 1;
      int dim_a = (map_to_regional.age_breakpoints.size() == 0) ? 1 : map_to_regional.age_breakpoints.size() + 1;

      gslMatrix subdesign(dim_a, map_to_regional.design_matrix.ncols());
      gsl_vector* intermediate_vec = gsl_vector_alloc(subdesign.nrows());

      select_design_matrix(subdesign, map_to_regional.design_matrix, dim_r == 1,
			   region_index * dim_a, dim_a);

      R_generalised_linear_regression(intermediate_vec, *subdesign,
				       param_value, map_to_regional.regression_link);

      // now fill out the out_vec according to the breakpoints
      vec_breakpoint_cut(out_vec, *map_to_regional.age_breakpoints, intermediate_vec);

      // free any memory allocated within this scope
      gsl_vector_free(intermediate_vec);
      //gsl_matrix_free(subdesign);

    }

}



// Block version - copy of old code
void regional_matrix_parameter(gsl_matrix* out_mat, gsl_vector_const_view param_value, const updRegrDef& map_to_regional, const int region_index, const int time_steps_per_day)
{

  /// IN THIS FUNCTION, THE CODE LINES INVOLVING THE subdesign MATRICES HAVE BEEN DONE IN A VERY
  /// SAFE MANNER, ONCE WORKING, THIS MAY BE SPEEDED UP - INVESTIGATE

  // A MORE GENERAL VERSION OF THE FUNCTION regional_vector_parameter WHICH ALLOWS FOR TEMPORAL BREAKPOINTS AS WELL AS
  // BREAKPOINTS OVER AGES.


  // if there's any temporal, regional or age variation then the design matrix should be non-zero
  if(map_to_regional.design_matrix[region_index].nrows() == 0)
    {
      // no variation, only first component of parameter value will be used.
      gsl_matrix_set_all(out_mat,
			 gsl_vector_get(&param_value.vector, 0));
    }
  else
    {
      // the design matrix should be set.
      // get the number of intervals over time (and age) for each region
      int dim_r = (map_to_regional.region_breakpoints.size() == 0) ? 1 : map_to_regional.region_breakpoints.size() + 1;
      int dim_a = (map_to_regional.age_breakpoints.size() == 0) ? 1 : map_to_regional.age_breakpoints.size() + 1;
      int dim_t = (map_to_regional.time_breakpoints.size() == 0) ? 1 : map_to_regional.time_breakpoints.size() + 1;

      gslMatrix subdesign(dim_t * dim_a, map_to_regional.design_matrix[region_index].ncols());
      gsl_vector* intermediate_vec = gsl_vector_alloc(subdesign.nrows());

      select_design_matrix(subdesign, map_to_regional.design_matrix[region_index], dim_r == 1,
			   region_index * dim_t * dim_a, dim_t * dim_a);

      R_generalised_linear_regression(intermediate_vec, *subdesign,
				      &param_value.vector, map_to_regional.regression_link);

      // now fill out the out_mat according to the breakpoints - first need to multiply up the breakpoints to the correct time scale
      gslVectorInt rescaled_temporal_breakpoints;
      if(dim_t > 1){
	//rescaled_temporal_breakpoints = gsl_vector_int_alloc(map_to_regional.time_breakpoints.size());
	//gsl_vector_int_memcpy(rescaled_temporal_breakpoints, map_to_regional.time_breakpoints);
	//gsl_vector_int_scale(rescaled_temporal_breakpoints, time_steps_per_day);
	rescaled_temporal_breakpoints = map_to_regional.time_breakpoints * time_steps_per_day;
      }
      mat_breakpoint_cut(out_mat, *rescaled_temporal_breakpoints, *map_to_regional.age_breakpoints, intermediate_vec);
  

      // free any memory allocated within this scope
      //if(rescaled_temporal_breakpoints != 0)
      //gsl_vector_int_free(rescaled_temporal_breakpoints);
      gsl_vector_free(intermediate_vec);
      //gsl_matrix_free(subdesign);

    }

}

void regional_matrix_parameter(gsl_matrix* out_mat, const gsl_vector* param_value, const regression_def& map_to_regional, const int region_index, const int time_steps_per_day)
{
  /// IN THIS FUNCTION, THE CODE LINES INVOLVING THE subdesign MATRICES HAVE BEEN DONE IN A VERY
  /// SAFE MANNER, ONCE WORKING, THIS MAY BE SPEEDED UP - INVESTIGATE

  // A MORE GENERAL VERSION OF THE FUNCTION regional_vector_parameter WHICH ALLOWS FOR TEMPORAL BREAKPOINTS AS WELL AS
  // BREAKPOINTS OVER AGES.


  // if there's any temporal, regional or age variation then the design matrix should be non-zero
  if(map_to_regional.design_matrix.nrows() == 0)
    {
      // no variation, only first component of parameter value will be used.
      gsl_matrix_set_all(out_mat,
			 gsl_vector_get(param_value, 0));
    }
  else
    {
      // the design matrix should be set.
      // get the number of intervals over time (and age) for each region
      int dim_r = (map_to_regional.region_breakpoints.size() == 0) ? 1 : map_to_regional.region_breakpoints.size() + 1;
      int dim_a = (map_to_regional.age_breakpoints.size() == 0) ? 1 : map_to_regional.age_breakpoints.size() + 1;
      int dim_t = (map_to_regional.time_breakpoints.size() == 0) ? 1 : map_to_regional.time_breakpoints.size() + 1;

      gslMatrix subdesign(dim_t * dim_a, map_to_regional.design_matrix.ncols());
      gsl_vector* intermediate_vec = gsl_vector_alloc(subdesign.nrows());

      select_design_matrix(subdesign, map_to_regional.design_matrix, dim_r == 1,
			   region_index * dim_t * dim_a, dim_t * dim_a);

      R_generalised_linear_regression(intermediate_vec, *subdesign,
				      param_value, map_to_regional.regression_link);

      // now fill out the out_mat according to the breakpoints - first need to multiply up the breakpoints to the correct time scale
      gslVectorInt rescaled_temporal_breakpoints;      
      if(dim_t > 1){
	//rescaled_temporal_breakpoints = gsl_vector_int_alloc(map_to_regional.time_breakpoints.size());
	//gsl_vector_int_memcpy(rescaled_temporal_breakpoints, map_to_regional.time_breakpoints);
	//gsl_vector_int_scale(rescaled_temporal_breakpoints, time_steps_per_day);
	rescaled_temporal_breakpoints = map_to_regional.time_breakpoints * time_steps_per_day;
      }
      mat_breakpoint_cut(out_mat, *rescaled_temporal_breakpoints, *map_to_regional.age_breakpoints, intermediate_vec);
  

      // free any memory allocated within this scope
      // if(rescaled_temporal_breakpoints != 0)
      //gsl_vector_int_free(rescaled_temporal_breakpoints);
      gsl_vector_free(intermediate_vec);
      //gsl_matrix_free(subdesign);

    }

}

double fn_initial_phase(const double phase_time, const double days_in_year)
{
  return cos(2 * M_PI * phase_time / days_in_year);
}

// SUB-FUNCTION FOR CALCULATING THE SEASONALITY OF R0 - THE STORING OF THE PHASE DIFFERENCES OVER THE YEAR
void Seasonal_R0_phase_difference(
				  gsl_vector* out_phase_differences,
				  const int start_day,
				  const double peak_day,
				  const double days_in_year,
				  const int time_steps_per_day)
{

  double initial_phase = fn_initial_phase((double) start_day - peak_day, days_in_year);
  for(register int int_time_steps = 0; int_time_steps < out_phase_differences->size; int_time_steps++)
    gsl_vector_set(out_phase_differences,
		   int_time_steps,
		   fn_initial_phase(((double) int_time_steps / ((double) time_steps_per_day)) + start_day - peak_day, days_in_year) - initial_phase
		   ); 

}

// SPECIAL CASE FOR THE MIXING MODEL'S PARAMETERISATION
void regional_mixmod_parameter(mixing_model& out_mix, const gsl_vector* param_value, const regression_def& map_to_regional, const mixing_model& base_mix, const int region_index, const gsl_vector* region_population, const bool evec_flag)
{
  int dim_r, dim_t, dim_a;

  // COPY THE ELEMENTS OF THE BASIC MIXING MODEL INTO THE OUTPUT MIXING MODEL
  mixing_model_memcpy(out_mix, base_mix);

  // HOPEFULLY ONLY REGIONAL VARIATION IS SPECIFIED WITHIN THE INPUT PARAMETER FILE
  if(map_to_regional.design_matrix.nrows() == 0)
    {
      // NO REGIONAL VARIATION - ALL PARAMETERS APPLY TO ALL REGIONS
      gsl_vector_memcpy(out_mix.scalants, param_value);
      dim_r = 1;
    }
  else
    {
      // the design matrix should be set
      // get the number of intervals over time (and age) for each region
      dim_r = (map_to_regional.region_breakpoints.size() == 0) ? 1 : map_to_regional.region_breakpoints.size() + 1;
      dim_a = (map_to_regional.age_breakpoints.size() == 0) ? 1 : map_to_regional.age_breakpoints.size() + 1;
      dim_t = (map_to_regional.time_breakpoints.size() == 0) ? 1 : map_to_regional.time_breakpoints.size() + 1;

      if((dim_a > 1) || (dim_t > 1))
	{
	  ERROR_INPUT_EXIT("Time and/or age breakpoints specified for contact parameters in function %s", "regional_mixmod_parameter\n");
	}

      // THE MIXING MODEL SHOULD HAVE ALL OF ITS MEMORY PRE-ALLOCATED. CHECK THAT THE SPECIFIED PARAMETER VECTOR AND THE
      // PARAMETERISATION MATRICES MATCH UP IN TERMS OF THE SPECIFIED DIMENSION - THIS SHOULD BE dim_r TIMES THE UPDATEABLE PARAMETER VALUE

      if((out_mix.scalants->size * dim_r) != param_value->size)
	{
	  ERROR_INPUT_EXIT("Mismatch between mixing matrix parameterisation matrices and the size of the contact parameter, caught in function %s", "regional_mixmod_parameter\n");
	}

      // USE THIS VALUE TO PICK OUT THE APPROPRIATE ROWS OF THE DESIGN MATRIX
      gslMatrix subdesign(out_mix.scalants->size, map_to_regional.design_matrix.ncols());
      gsl_vector* temp_subvec = gsl_vector_alloc(subdesign.nrows());
      
      // gsl_vector_view temp_subvec = gsl_vector_subvector(param_value, region_index * out_mix.scalants->size, out_mix.scalants->size);
      select_design_matrix(subdesign,
			   map_to_regional.design_matrix,
			   dim_r == 1,
			   region_index * out_mix.scalants->size,
			   out_mix.scalants->size);
      R_generalised_linear_regression(temp_subvec, *subdesign, param_value, map_to_regional.regression_link);
      
      gsl_vector_memcpy(out_mix.scalants, temp_subvec);

      gsl_vector_free(temp_subvec);
      //gsl_matrix_free(subdesign);

    }

  mixing_matrix_parameterise(out_mix);
  mixing_matrix_scale_and_normalisation(out_mix, region_population, evec_flag);

}


// Block update - no need for submatrix
void regional_mixmod_parameter(mixing_model& out_mix, const updParamSet &paramSet, const upd::paramIndex index, const mixing_model& base_mix, const int region_index, const gsl_vector* region_population, const bool evec_flag)
{
  int dim_r, dim_t, dim_a;

  const updRegrDef &map = paramSet[index].map_to_regional;
  gsl_vector_const_view view = paramSet.lookup(index, region_index);
  
  // COPY THE ELEMENTS OF THE BASIC MIXING MODEL INTO THE OUTPUT MIXING MODEL
  mixing_model_memcpy(out_mix, base_mix);

  // HOPEFULLY ONLY REGIONAL VARIATION IS SPECIFIED WITHIN THE INPUT PARAMETER FILE
  if(map.design_matrix[region_index].nrows() == 0)
    {
      // NO REGIONAL VARIATION - ALL PARAMETERS APPLY TO ALL REGIONS
      gsl_vector_memcpy(out_mix.scalants, &view.vector);
      dim_r = 1;
    }
  else
    {
      // the design matrix should be set
      // get the number of intervals over time (and age) for each region
      //dim_r = (map_to_regional.region_breakpoints.size() == 0) ? 1 : map_to_regional.region_breakpoints.size() + 1;
      dim_a = (map.age_breakpoints.size() == 0) ? 1 : map.age_breakpoints.size() + 1;
      dim_t = (map.time_breakpoints.size() == 0) ? 1 : map.time_breakpoints.size() + 1;

      if((dim_a > 1) || (dim_t > 1))
	{
	  ERROR_INPUT_EXIT("Time and/or age breakpoints specified for contact parameters in function %s", "regional_mixmod_parameter\n");
	}

      // THE MIXING MODEL SHOULD HAVE ALL OF ITS MEMORY PRE-ALLOCATED. CHECK THAT THE SPECIFIED PARAMETER VECTOR AND THE
      // PARAMETERISATION MATRICES MATCH UP IN TERMS OF THE SPECIFIED DIMENSION - THIS SHOULD BE dim_r TIMES THE UPDATEABLE PARAMETER VALUE

      /*
      if((out_mix.scalants->size * dim_r) != param_value->size)
	{
	  ERROR_INPUT_EXIT("Mismatch between mixing matrix parameterisation matrices and the size of the contact parameter, caught in function %s", "regional_mixmod_parameter\n");
	}

      
      // USE THIS VALUE TO PICK OUT THE APPROPRIATE ROWS OF THE DESIGN MATRIX
      gslMatrix subdesign(out_mix.scalants->size, map_to_regional.design_matrix.ncols());
      
      // gsl_vector_view temp_subvec = gsl_vector_subvector(param_value, region_index * out_mix.scalants->size, out_mix.scalants->size);
      select_design_matrix(subdesign,
			   map_to_regional.design_matrix,
			   dim_r == 1,
			   region_index * out_mix.scalants->size,
			   out_mix.scalants->size);
      */
      gsl_vector* temp_subvec = gsl_vector_alloc(map.design_matrix[region_index].nrows());

      R_generalised_linear_regression(temp_subvec, *map.design_matrix[region_index], &view.vector, map.regression_link);
      
      gsl_vector_memcpy(out_mix.scalants, temp_subvec);

      gsl_vector_free(temp_subvec);
      //gsl_matrix_free(subdesign);

    }

  mixing_matrix_parameterise(out_mix);
  mixing_matrix_scale_and_normalisation(out_mix, region_population, evec_flag);

}

// Block update - no need for submatrix

// Wrapper for actional regional_matrix_parameter that sorts out conversion of view to vector
/*
void regional_matrix_parameter(gsl_matrix* out_mat, const updParamSet &paramSet, const upd::paramIndex index, const int region_index, const int time_steps_per_day) {

  gsl_vector_const_view view = paramSet.lookup(index, region_index);

  regional_matrix_parameter(out_mat,
			    &(view.vector),
			    paramSet[index].map_to_regional,
			    region_index,
			    time_steps_per_day);
}
*/
/*
// Note: Actual GP_PATCH method defined in gp_param_patch.h/cc
void regional_matrix_param_GP_PATCH(gsl_matrix* out_mat, const updParamSet &paramSet, const upd::paramIndex index, const int region_index, const int time_steps_per_day) {

  gsl_vector_const_view view = paramSet.lookup(index, region_index);

  regional_matrix_parameter_GP_PATCH(out_mat,
			    &(view.vector),
			    paramSet[index].map_to_regional,
			    region_index,
			    time_steps_per_day);
}
*/

void block_regional_parameters(regional_model_params& out_rmp,
			       const updParamSet &updPars,
			       const global_model_instance_parameters& in_gmip,
			       const int& region_index,
			       const gsl_vector* population_by_age,
			       const double& total_population_size,
			       const mixing_model& in_mix,
			       flagclass& update_flags)
{

  // CCS:
  // Modify to use new param object
  // 
  // Regional params use hardcoded magic numbers to access the specific model
  // params. Unlikely to be any fast alternative.
  //
  // BUT: getFlag() also very inefficient text lookup.
  // Enum updateable_parameter_index defined in RTM_StructAssign.h
  // 
  // TODO: Rewrite regional_model_parameters to either have named flags, or
  // at the very least a vector using the same hardcoded enum.
  
  // Idea: A templated parameter class with flag, data and function ptr or
  // similar to define how each param is updated. Then just iterate over
  // vector of region_param<T>. UMP lookup can be managed by doing string
  // search at load time and storing a vector index inside region_param.
  // Will also need to store flag indicating local or global UMP.
  
  // The only things used are the param_value and map_to_regional.
  // (Why is map_to_regional a UMP member, when UMP may exist that doesn't
  // update a regional param? Shouldn't this be regional param member?)
  

  // This function should also be parallelisable over regions

  if(update_flags.getFlag("l_gp_negbin_overdispersion")){
    regional_matrix_parameter(out_rmp.l_gp_negbin_overdispersion, updPars.lookup(upd::GP_OVERDISP, region_index), updPars[upd::GP_OVERDISP].map_to_regional, region_index, 1); // Only want this per day rather per delta t
    //regional_matrix_parameter(out_rmp.l_gp_negbin_overdispersion, updPars, upd::GP_OVERDISP, region_index, 1);  // Only want this per day rather per delta t
  }
  if(update_flags.getFlag("l_hosp_negbin_overdispersion")){
    regional_matrix_parameter(out_rmp.l_hosp_negbin_overdispersion, updPars.lookup(upd::HOSP_OVERDISP, region_index), updPars[upd::HOSP_OVERDISP].map_to_regional, region_index, 1);
    //regional_matrix_parameter(out_rmp.l_hosp_negbin_overdispersion, updPars, upd::HOSP_OVERDISP, region_index, 1);
    // Again, only valid per observation, rather than per timestep.
  }
  if(update_flags.getFlag("l_day_of_week_effect"))
    regional_matrix_parameter(out_rmp.l_day_of_week_effect, updPars.lookup(upd::DOW_EFFECTS, region_index), updPars[upd::DOW_EFFECTS].map_to_regional, region_index, 1); // Only want this per day rather per delta t
  if(update_flags.getFlag("l_init_prop_sus"))
    regional_vector_parameter(out_rmp.l_init_prop_sus, updPars.lookup(upd::PROP_SUS, region_index), updPars[upd::PROP_SUS].map_to_regional, region_index);
  if(update_flags.getFlag("l_init_prop_sus_HI_geq_32"))
    regional_vector_parameter(out_rmp.l_init_prop_sus_HI_geq_32, updPars.lookup(upd::PROP_HI_GEQ_32, region_index), updPars[upd::PROP_HI_GEQ_32].map_to_regional, region_index);
  if(update_flags.getFlag("l_average_infectious_period"))
    {
      regional_matrix_parameter(out_rmp.l_average_infectious_period, updPars.lookup(upd::AIP, region_index), updPars[upd::AIP].map_to_regional, region_index, in_gmip.l_transmission_time_steps_per_day);
      gsl_matrix_add_constant(out_rmp.l_average_infectious_period, MIN_DELTA_T_TO_LENGTH_OF_STAY_RATIO / ((double) in_gmip.l_transmission_time_steps_per_day));
    }
  if(update_flags.getFlag("l_r1_period"))
    {
      regional_matrix_parameter(out_rmp.l_r1_period, updPars.lookup(upd::AR1, region_index), updPars[upd::AR1].map_to_regional, region_index, in_gmip.l_transmission_time_steps_per_day);
      gsl_matrix_add_constant(out_rmp.l_r1_period, MIN_DELTA_T_TO_LENGTH_OF_STAY_RATIO / (2 * (double) in_gmip.l_transmission_time_steps_per_day));
    }
  if(update_flags.getFlag("l_latent_period"))
    {
      regional_matrix_parameter(out_rmp.l_latent_period, updPars.lookup(upd::ALP, region_index), updPars[upd::ALP].map_to_regional, region_index, in_gmip.l_transmission_time_steps_per_day);
      gsl_matrix_add_constant(out_rmp.l_latent_period, MIN_DELTA_T_TO_LENGTH_OF_STAY_RATIO / ((double) in_gmip.l_transmission_time_steps_per_day));
    }
  if(update_flags.getFlag("l_waning_period"))
    {
      regional_matrix_parameter(out_rmp.l_waning_period, updPars.lookup(upd::IWAN, region_index), updPars[upd::IWAN].map_to_regional, region_index, in_gmip.l_transmission_time_steps_per_day);
      gsl_matrix_add_constant(out_rmp.l_waning_period, MIN_DELTA_T_TO_LENGTH_OF_STAY_RATIO / ((double) in_gmip.l_transmission_time_steps_per_day));
    }
  if(update_flags.getFlag("l_vacc1_disease"))
    regional_matrix_parameter(out_rmp.l_vacc1_disease, updPars.lookup(upd::VAC1_DISEASE, region_index), updPars[upd::VAC1_DISEASE].map_to_regional, region_index, in_gmip.l_transmission_time_steps_per_day);
  if(update_flags.getFlag("l_vaccn_disease"))
    regional_matrix_parameter(out_rmp.l_vaccn_disease, updPars.lookup(upd::VACN_DISEASE, region_index), updPars[upd::VACN_DISEASE].map_to_regional, region_index, in_gmip.l_transmission_time_steps_per_day);
  if(update_flags.getFlag("l_vaccb_disease"))
    regional_matrix_parameter(out_rmp.l_vaccb_disease, updPars.lookup(upd::VACB_DISEASE, region_index), updPars[upd::VACB_DISEASE].map_to_regional, region_index, in_gmip.l_transmission_time_steps_per_day);
  if(update_flags.getFlag("l_vacc4_disease"))
    regional_matrix_parameter(out_rmp.l_vacc4_disease, updPars.lookup(upd::VAC4_DISEASE, region_index), updPars[upd::VAC4_DISEASE].map_to_regional, region_index, in_gmip.l_transmission_time_steps_per_day);
  if(update_flags.getFlag("l_vacc1_infect"))
    regional_matrix_parameter(out_rmp.l_vacc1_infect, updPars.lookup(upd::VAC1_INFECT, region_index), updPars[upd::VAC1_INFECT].map_to_regional, region_index, in_gmip.l_transmission_time_steps_per_day);
  if(update_flags.getFlag("l_vaccn_infect"))
    regional_matrix_parameter(out_rmp.l_vaccn_infect, updPars.lookup(upd::VACN_INFECT, region_index), updPars[upd::VACN_INFECT].map_to_regional, region_index, in_gmip.l_transmission_time_steps_per_day);
  if(update_flags.getFlag("l_vaccb_infect"))
    regional_matrix_parameter(out_rmp.l_vaccb_infect, updPars.lookup(upd::VACB_INFECT, region_index), updPars[upd::VACB_INFECT].map_to_regional, region_index, in_gmip.l_transmission_time_steps_per_day);
  if(update_flags.getFlag("l_vacc4_infect"))
    regional_matrix_parameter(out_rmp.l_vacc4_infect, updPars.lookup(upd::VAC4_INFECT, region_index), updPars[upd::VAC4_INFECT].map_to_regional, region_index, in_gmip.l_transmission_time_steps_per_day);
  
  if(update_flags.getFlag("l_relative_infectiousness_I2_wrt_I1"))
    regional_matrix_parameter(out_rmp.l_relative_infectiousness_I2_wrt_I1, updPars.lookup(upd::REL_INFECT, region_index), updPars[upd::REL_INFECT].map_to_regional, region_index, in_gmip.l_transmission_time_steps_per_day);
  if(update_flags.getFlag("l_lbeta_rw")){
    regional_time_vector_parameter(out_rmp.l_lbeta_rw, updPars, upd::LBETA_RW, region_index, in_gmip.l_transmission_time_steps_per_day);
    if(0 != gsl_vector_get(out_rmp.l_lbeta_rw, 0))
      gsl_vector_scale(out_rmp.l_lbeta_rw, 1 / gsl_vector_get(out_rmp.l_lbeta_rw, 0)); // Random-walk necessarily centred on 0.
  }
  if(update_flags.getFlag("l_sensitivity"))
    { // STRICTLY A SCALAR QUANTITY
      gsl_vector_const_view view = updPars.lookup(upd::SENS, region_index);
      fixed_quantity(out_rmp.l_sensitivity,
		     gsl_vector_get(&view.vector, 0),
		     updPars[upd::SENS].param_name,
		     updPars[upd::SENS].map_to_regional.design_matrix[region_index]);
    }
  if(update_flags.getFlag("l_specificity"))
    { // STRICTLY A SCALAR QUANTITY
      gsl_vector_const_view view = updPars.lookup(upd::SPEC, region_index);
      fixed_quantity(out_rmp.l_specificity,
		     gsl_vector_get(&view.vector, 0),
		     updPars[upd::SPEC].param_name,
		     updPars[upd::SPEC].map_to_regional.design_matrix[region_index]);
    }
  if(update_flags.getFlag("l_sero_sensitivity"))
    { // NO LONGER STRICTLY A SCALAR QUANTITY
      regional_matrix_parameter(out_rmp.l_sero_sensitivity, updPars.lookup(upd::SSENS, region_index), updPars[upd::SSENS].map_to_regional, region_index, 1);
      // gsl_vector_const_view view = updPars.lookup(upd::SSENS, region_index);
      //   fixed_quantity(out_rmp.l_sero_sensitivity,
      // 		     gsl_vector_get(&view.vector, 0),
      // 		     updPars[upd::SSENS].param_name,
      // 		     updPars[upd::SSENS].map_to_regional.design_matrix[region_index]);
    }
  if(update_flags.getFlag("l_sero_specificity"))
    { // NO LONGER STRICTLY A SCALAR QUANTITY
      regional_matrix_parameter(out_rmp.l_sero_specificity, updPars.lookup(upd::SSPEC, region_index), updPars[upd::SSPEC].map_to_regional, region_index, 1);
      // gsl_vector_const_view view = updPars.lookup(upd::SSPEC, region_index);
      //   fixed_quantity(out_rmp.l_sero_specificity,
      // 		     gsl_vector_get(&view.vector, 0),
      // 		     updPars[upd::SSPEC].param_name,
      // 		     updPars[upd::SSPEC].map_to_regional.design_matrix[region_index]);
    }
  if(update_flags.getFlag("l_pr_symp"))
    regional_matrix_parameter(out_rmp.l_pr_symp, updPars.lookup(upd::PROP_SYMP, region_index), updPars[upd::PROP_SYMP].map_to_regional, region_index, in_gmip.l_reporting_time_steps_per_day);
  if(update_flags.getFlag("l_pr_onset_to_GP"))
    { // PATCHED CODE!!!
      if(in_gmip.l_GP_patch_flag)
	regional_matrix_parameter_GP_PATCH(out_rmp.l_pr_onset_to_GP, updPars, upd::PROP_GP, region_index, in_gmip.l_reporting_time_steps_per_day);	
      else regional_matrix_parameter(out_rmp.l_pr_onset_to_GP, updPars.lookup(upd::PROP_GP, region_index), updPars[upd::PROP_GP].map_to_regional, region_index, in_gmip.l_reporting_time_steps_per_day);
    } // PATCHED CODE ENDS
  if(update_flags.getFlag("l_pr_onset_to_Hosp"))
    regional_matrix_parameter(out_rmp.l_pr_onset_to_Hosp, updPars.lookup(upd::PROP_HOSP, region_index), updPars[upd::PROP_HOSP].map_to_regional, region_index, in_gmip.l_reporting_time_steps_per_day);
  if(update_flags.getFlag("l_pr_onset_to_Death"))
    regional_matrix_parameter(out_rmp.l_pr_onset_to_Death, updPars.lookup(upd::PROP_DEATH, region_index), updPars[upd::PROP_DEATH].map_to_regional, region_index, in_gmip.l_reporting_time_steps_per_day);
  if(update_flags.getFlag("l_importation_rate"))
    regional_matrix_parameter(out_rmp.l_importation_rate, updPars.lookup(upd::IMPORTATION, region_index), updPars[upd::IMPORTATION].map_to_regional, region_index, in_gmip.l_transmission_time_steps_per_day);
  if(update_flags.getFlag("l_MIXMOD")) {
    //gsl_vector_const_view view = updPars.lookupValue(upd::CONTACT, region_index);
    regional_mixmod_parameter(out_rmp.l_MIXMOD, updPars, upd::CONTACT, in_mix, region_index, population_by_age, true);
  }
  if(update_flags.getFlag("l_background_gps_counts"))
    {
      regional_matrix_parameter(out_rmp.l_background_gps_counts, updPars.lookup(upd::BGR, region_index), updPars[upd::BGR].map_to_regional, region_index, 1);
      // THE PARAMETER VECTOR GIVES CONSULTATION RATES. CONVERT THESE HERE TO BACKGROUND GP CONSULTATION TOTALS FOR EACH REGION.
      // RATES ARE NUMBER OF CONSULTATIONS PER 100K PEOPLE PER YEAR
      // FIRST WRITE AS NUMBER OF CONSULTATIONS PER PERSON PER DAY
      gsl_matrix_scale(out_rmp.l_background_gps_counts, 1 / (100000 * DAYS_IN_YEAR));
      // NOW MULTIPLY BY THE POPULATION VECTOR
      gsl_matrix_multiplication_vector_product_by_row_dbl(out_rmp.l_background_gps_counts, population_by_age);
      // out_rmp.l_background_gps_counts should now contain absolute values NOT rates
    }
  if(update_flags.getFlag("l_R0_peakday"))
    {
      if(updPars[upd::R0_PEAKDAY].map_to_regional.design_matrix[region_index].nrows() > 0)
	{
	  ERROR_PARAM_INPUT(updPars[upd::R0_PEAKDAY].param_name.c_str());
	}
      out_rmp.l_R0_peakday = updPars.lookup0(upd::R0_PEAKDAY, region_index);
    }
  if(update_flags.getFlag("l_EGR"))
    {
      // gsl_vector_const_view view = updPars.lookupValue(upd::EGR, region_index);
      regional_scalar_parameter(out_rmp.l_EGR, updPars, upd::EGR, region_index);
      //      int temp_region_index = (updPars.lookup(upd::EGR).map_to_regional.region_breakpoints == 0) ? 0 : region_index;
      //      out_rmp.l_EGR = gsl_vector_get(updPars.lookup(upd::EGR, region_index), temp_region_index);

    
    }
  // VALUES DERIVED FROM OTHER PARAMETERS
  if(update_flags.getFlag("l_R0_init"))
    {
      double d_A = gsl_matrix_get(out_rmp.l_average_infectious_period, 0, 0);
      double d_L = gsl_matrix_get(out_rmp.l_latent_period, 0, 0);

      out_rmp.l_R0_init = out_rmp.l_EGR * d_A * gsl_pow_2(((out_rmp.l_EGR * d_L / 2) + 1)) / (1 - (1 / (gsl_pow_2(((out_rmp.l_EGR * d_A) / 2) + 1))));

    }
  if(update_flags.getFlag("l_I0"))
    {
      double d_A = gsl_matrix_get(out_rmp.l_average_infectious_period, 0, 0);
      //int temp_region_index = (updPars[upd::LPL0].map_to_regional.region_breakpoints.size() == 0) ? 0 : region_index;

      // Do we need to check whether LPL0 is regional? Or is it safe to take first
      // value regardless?
      double val = updPars.lookup0(upd::LPL0, region_index);
      double nu = gsl_sf_exp(val);
      out_rmp.l_I0 = total_population_size * d_A * nu / (gsl_matrix_get(out_rmp.l_pr_onset_to_GP, 0, 0) * out_rmp.l_R0_init);
    }
  if(update_flags.getFlag("d_R0_phase_differences"))
    Seasonal_R0_phase_difference(out_rmp.d_R0_phase_differences,
				 in_gmip.l_day_of_start,
				 out_rmp.l_R0_peakday,
				 DAYS_IN_YEAR,
				 in_gmip.l_transmission_time_steps_per_day);
  if(update_flags.getFlag("l_R0_Amplitude"))
    { // STRICTLY A SCALAR QUANTITY
      if(updPars[upd::R0_AMP].map_to_regional.design_matrix[region_index].nrows() > 0)
	{
	  ERROR_PARAM_INPUT(updPars[upd::R0_AMP].param_name.c_str());
	}
      out_rmp.l_R0_Amplitude = updPars.lookup0(upd::R0_AMP, region_index) * out_rmp.l_R0_init / (1 + fn_initial_phase(in_gmip.l_day_of_start - out_rmp.l_R0_peakday, DAYS_IN_YEAR));
    }
  

}



void evaluate_regional_parameters(regional_model_params& out_rmp, const updateable_model_parameter* in_umps,
				  const global_model_instance_parameters& in_gmip, const int& region_index,
				  const gsl_vector* population_by_age, const double& total_population_size, const mixing_model& in_mix,
				  flagclass& update_flags)
{
    
  if(update_flags.getFlag("l_gp_negbin_overdispersion")){
    regional_matrix_parameter(out_rmp.l_gp_negbin_overdispersion, in_umps[GP_OVERDISP_INDEX].param_value, in_umps[GP_OVERDISP_INDEX].map_to_regional, region_index, 1); // Only want this per day rather per delta t
  }
  if(update_flags.getFlag("l_hosp_negbin_overdispersion")){
    regional_matrix_parameter(out_rmp.l_hosp_negbin_overdispersion, in_umps[HOSP_OVERDISP_INDEX].param_value, in_umps[HOSP_OVERDISP_INDEX].map_to_regional, region_index, 1); // Again, only valid per observation, rather than per timestep.
  }
  if(update_flags.getFlag("l_day_of_week_effect"))
    regional_matrix_parameter(out_rmp.l_day_of_week_effect, in_umps[DOW_EFFECTS_INDEX].param_value, in_umps[DOW_EFFECTS_INDEX].map_to_regional, region_index, 1); // Only want this per day rather per delta t
  if(update_flags.getFlag("l_init_prop_sus"))
    regional_vector_parameter(out_rmp.l_init_prop_sus, in_umps[PROP_SUS_INDEX].param_value, in_umps[PROP_SUS_INDEX].map_to_regional, region_index);
  if(update_flags.getFlag("l_init_prop_sus_HI_geq_32"))
    regional_vector_parameter(out_rmp.l_init_prop_sus_HI_geq_32, in_umps[PROP_HI_GEQ_32_INDEX].param_value, in_umps[PROP_HI_GEQ_32_INDEX].map_to_regional, region_index);
  if(update_flags.getFlag("l_average_infectious_period"))
    {
      regional_matrix_parameter(out_rmp.l_average_infectious_period, in_umps[AIP_INDEX].param_value, in_umps[AIP_INDEX].map_to_regional, region_index, in_gmip.l_transmission_time_steps_per_day);
      gsl_matrix_add_constant(out_rmp.l_average_infectious_period, MIN_DELTA_T_TO_LENGTH_OF_STAY_RATIO / ((double) in_gmip.l_transmission_time_steps_per_day));
    }
  if(update_flags.getFlag("l_r1_period"))
    {
      regional_matrix_parameter(out_rmp.l_r1_period, in_umps[AR1_INDEX].param_value, in_umps[AR1_INDEX].map_to_regional, region_index, in_gmip.l_transmission_time_steps_per_day);
      gsl_matrix_add_constant(out_rmp.l_r1_period, MIN_DELTA_T_TO_LENGTH_OF_STAY_RATIO / (2 * (double) in_gmip.l_transmission_time_steps_per_day));
    }
  if(update_flags.getFlag("l_latent_period"))
    {
      regional_matrix_parameter(out_rmp.l_latent_period, in_umps[ALP_INDEX].param_value, in_umps[ALP_INDEX].map_to_regional, region_index, in_gmip.l_transmission_time_steps_per_day);
      gsl_matrix_add_constant(out_rmp.l_latent_period, MIN_DELTA_T_TO_LENGTH_OF_STAY_RATIO / ((double) in_gmip.l_transmission_time_steps_per_day));
    }
  if(update_flags.getFlag("l_waning_period"))
    {
      regional_matrix_parameter(out_rmp.l_waning_period, in_umps[IWAN_INDEX].param_value, in_umps[IWAN_INDEX].map_to_regional, region_index, in_gmip.l_transmission_time_steps_per_day);
      gsl_matrix_add_constant(out_rmp.l_waning_period, MIN_DELTA_T_TO_LENGTH_OF_STAY_RATIO / ((double) in_gmip.l_transmission_time_steps_per_day));
    }
  if(update_flags.getFlag("l_vacc1_disease"))
    regional_matrix_parameter(out_rmp.l_vacc1_disease, in_umps[VAC1_DISEASE_INDEX].param_value, in_umps[VAC1_DISEASE_INDEX].map_to_regional, region_index, in_gmip.l_transmission_time_steps_per_day);
  if(update_flags.getFlag("l_vaccn_disease"))
    regional_matrix_parameter(out_rmp.l_vaccn_disease, in_umps[VACN_DISEASE_INDEX].param_value, in_umps[VACN_DISEASE_INDEX].map_to_regional, region_index, in_gmip.l_transmission_time_steps_per_day);
  if(update_flags.getFlag("l_vaccb_disease"))
    regional_matrix_parameter(out_rmp.l_vaccb_disease, in_umps[VACB_DISEASE_INDEX].param_value, in_umps[VACB_DISEASE_INDEX].map_to_regional, region_index, in_gmip.l_transmission_time_steps_per_day);
  if(update_flags.getFlag("l_vacc4_disease"))
    regional_matrix_parameter(out_rmp.l_vacc4_disease, in_umps[VAC4_DISEASE_INDEX].param_value, in_umps[VAC4_DISEASE_INDEX].map_to_regional, region_index, in_gmip.l_transmission_time_steps_per_day);
  if(update_flags.getFlag("l_vacc1_infect"))
    regional_matrix_parameter(out_rmp.l_vacc1_infect, in_umps[VAC1_INFECT_INDEX].param_value, in_umps[VAC1_INFECT_INDEX].map_to_regional, region_index, in_gmip.l_transmission_time_steps_per_day);
  if(update_flags.getFlag("l_vaccn_infect"))
 regional_matrix_parameter(out_rmp.l_vaccn_infect, in_umps[VACN_INFECT_INDEX].param_value, in_umps[VACN_INFECT_INDEX].map_to_regional, region_index, in_gmip.l_transmission_time_steps_per_day);
  if(update_flags.getFlag("l_vaccb_infect"))
 regional_matrix_parameter(out_rmp.l_vaccb_infect, in_umps[VACB_INFECT_INDEX].param_value, in_umps[VACB_INFECT_INDEX].map_to_regional, region_index, in_gmip.l_transmission_time_steps_per_day);
  if(update_flags.getFlag("l_vacc4_infect"))
 regional_matrix_parameter(out_rmp.l_vaccb_infect, in_umps[VAC4_INFECT_INDEX].param_value, in_umps[VAC4_INFECT_INDEX].map_to_regional, region_index, in_gmip.l_transmission_time_steps_per_day);
  if(update_flags.getFlag("l_relative_infectiousness_I2_wrt_I1"))
    regional_matrix_parameter(out_rmp.l_relative_infectiousness_I2_wrt_I1, in_umps[REL_INFECT_INDEX].param_value, in_umps[REL_INFECT_INDEX].map_to_regional, region_index, in_gmip.l_transmission_time_steps_per_day);
  if(update_flags.getFlag("l_lbeta_rw")){
    regional_time_vector_parameter(out_rmp.l_lbeta_rw, in_umps[LBETA_RW_INDEX].param_value, in_umps[LBETA_RW_INDEX].map_to_regional, region_index, in_gmip.l_transmission_time_steps_per_day);
    if(0 != gsl_vector_get(out_rmp.l_lbeta_rw, 0))
      gsl_vector_scale(out_rmp.l_lbeta_rw, 1 / gsl_vector_get(out_rmp.l_lbeta_rw, 0)); // Random-walk necessarily centred on 0.
  }
  if(update_flags.getFlag("l_sensitivity"))
    { // STRICTLY A SCALAR QUANTITY
      fixed_quantity(out_rmp.l_sensitivity,
		     gsl_vector_get(in_umps[SENS_INDEX].param_value, 0),
		     in_umps[SENS_INDEX].param_name,
		     in_umps[SENS_INDEX].map_to_regional.design_matrix);
    }
  if(update_flags.getFlag("l_specificity"))
    { // STRICTLY A SCALAR QUANTITY
      fixed_quantity(out_rmp.l_specificity,
		     gsl_vector_get(in_umps[SPEC_INDEX].param_value, 0),
		     in_umps[SPEC_INDEX].param_name,
		     in_umps[SPEC_INDEX].map_to_regional.design_matrix);
    }
  if(update_flags.getFlag("l_sero_sensitivity"))
    { // NO LONGER STRICTLY A SCALAR QUANTITY
      regional_matrix_parameter(out_rmp.l_sero_sensitivity, in_umps[SSENS_INDEX].param_value, in_umps[SSENS_INDEX].map_to_regional, region_index, 1);
      // fixed_quantity(out_rmp.l_sero_sensitivity,
      // 		     gsl_vector_get(in_umps[SSENS_INDEX].param_value, 0),
      // 		     in_umps[SSENS_INDEX].param_name,
      // 		     in_umps[SSENS_INDEX].map_to_regional.design_matrix);
    }
  if(update_flags.getFlag("l_sero_specificity"))
    { // NO LONGER STRICTLY A SCALAR QUANTITY
      regional_matrix_parameter(out_rmp.l_sero_specificity, in_umps[SSPEC_INDEX].param_value, in_umps[SSPEC_INDEX].map_to_regional, region_index, 1);
      // fixed_quantity(out_rmp.l_sero_specificity,
      // 		     gsl_vector_get(in_umps[SSPEC_INDEX].param_value, 0),
      // 		     in_umps[SSPEC_INDEX].param_name,
      // 		     in_umps[SSPEC_INDEX].map_to_regional.design_matrix);
    } 
  if(update_flags.getFlag("l_pr_symp"))
    regional_matrix_parameter(out_rmp.l_pr_symp, in_umps[PROP_SYMP_INDEX].param_value, in_umps[PROP_SYMP_INDEX].map_to_regional, region_index, in_gmip.l_reporting_time_steps_per_day);
  if(update_flags.getFlag("l_pr_onset_to_GP"))
    { // PATCHED CODE!!!
      if(in_gmip.l_GP_patch_flag)
	regional_matrix_parameter_GP_PATCH(out_rmp.l_pr_onset_to_GP, in_umps[PROP_GP_INDEX].param_value, in_umps[PROP_GP_INDEX].map_to_regional, region_index, in_gmip.l_reporting_time_steps_per_day);	
      else regional_matrix_parameter(out_rmp.l_pr_onset_to_GP, in_umps[PROP_GP_INDEX].param_value, in_umps[PROP_GP_INDEX].map_to_regional, region_index, in_gmip.l_reporting_time_steps_per_day);
    } // PATCHED CODE ENDS
  if(update_flags.getFlag("l_pr_onset_to_Hosp"))
    regional_matrix_parameter(out_rmp.l_pr_onset_to_Hosp, in_umps[PROP_HOSP_INDEX].param_value, in_umps[PROP_HOSP_INDEX].map_to_regional, region_index, in_gmip.l_reporting_time_steps_per_day);
  if(update_flags.getFlag("l_pr_onset_to_Death"))
    regional_matrix_parameter(out_rmp.l_pr_onset_to_Death, in_umps[PROP_DEATH_INDEX].param_value, in_umps[PROP_DEATH_INDEX].map_to_regional, region_index, in_gmip.l_reporting_time_steps_per_day);
  if(update_flags.getFlag("l_importation_rate"))
    regional_matrix_parameter(out_rmp.l_importation_rate, in_umps[IMPORTATION_INDEX].param_value, in_umps[IMPORTATION_INDEX].map_to_regional, region_index, in_gmip.l_transmission_time_steps_per_day);
  if(update_flags.getFlag("l_MIXMOD"))
    regional_mixmod_parameter(out_rmp.l_MIXMOD, in_umps[CONTACT_INDEX].param_value, in_umps[CONTACT_INDEX].map_to_regional, in_mix, region_index, population_by_age, true);
  if(update_flags.getFlag("l_background_gps_counts"))
    {
      regional_matrix_parameter(out_rmp.l_background_gps_counts, in_umps[BGR_INDEX].param_value, in_umps[BGR_INDEX].map_to_regional, region_index, 1);
      // THE PARAMETER VECTOR GIVES CONSULTATION RATES. CONVERT THESE HERE TO BACKGROUND GP CONSULTATION TOTALS FOR EACH REGION.
      // RATES ARE NUMBER OF CONSULTATIONS PER 100K PEOPLE PER YEAR
      // FIRST WRITE AS NUMBER OF CONSULTATIONS PER PERSON PER DAY
      gsl_matrix_scale(out_rmp.l_background_gps_counts, 1 / (100000 * DAYS_IN_YEAR));
      // NOW MULTIPLY BY THE POPULATION VECTOR
      gsl_matrix_multiplication_vector_product_by_row_dbl(out_rmp.l_background_gps_counts, population_by_age);
      // out_rmp.l_background_gps_counts should now contain absolute values NOT rates
    }
  if(update_flags.getFlag("l_R0_peakday"))
    {
      if(in_umps[R0_PEAKDAY_INDEX].map_to_regional.design_matrix.nrows() > 0)
	{
	  ERROR_PARAM_INPUT(in_umps[R0_PEAKDAY_INDEX].param_name.c_str());
	}
      out_rmp.l_R0_peakday = gsl_vector_get(in_umps[R0_PEAKDAY_INDEX].param_value, 0);
    }
  if(update_flags.getFlag("l_EGR"))
    {
      regional_scalar_parameter(out_rmp.l_EGR, in_umps[EGR_INDEX].param_value, in_umps[EGR_INDEX].map_to_regional, region_index);
      //      int temp_region_index = (in_umps[EGR_INDEX].map_to_regional.region_breakpoints == 0) ? 0 : region_index;
      //      out_rmp.l_EGR = gsl_vector_get(in_umps[EGR_INDEX].param_value, temp_region_index);
    }
  // VALUES DERIVED FROM OTHER PARAMETERS
  if(update_flags.getFlag("l_R0_init"))
    {
      double d_A = gsl_matrix_get(out_rmp.l_average_infectious_period, 0, 0);
      double d_L = gsl_matrix_get(out_rmp.l_latent_period, 0, 0);

      out_rmp.l_R0_init = out_rmp.l_EGR * d_A * gsl_pow_2(((out_rmp.l_EGR * d_L / 2) + 1)) / (1 - (1 / (gsl_pow_2(((out_rmp.l_EGR * d_A) / 2) + 1))));

    }
  if(update_flags.getFlag("l_I0"))
    {
      double d_A = gsl_matrix_get(out_rmp.l_average_infectious_period, 0, 0);
      int temp_region_index = (in_umps[LPL0_INDEX].map_to_regional.region_breakpoints.size() == 0) ? 0 : region_index;
      double val = gsl_vector_get(in_umps[LPL0_INDEX].param_value, temp_region_index);
      double nu = gsl_sf_exp(val);
      out_rmp.l_I0 = total_population_size * d_A * nu / (gsl_matrix_get(out_rmp.l_pr_onset_to_GP, 0, 0) * out_rmp.l_R0_init);
    }
  if(update_flags.getFlag("d_R0_phase_differences"))
    Seasonal_R0_phase_difference(out_rmp.d_R0_phase_differences,
				 in_gmip.l_day_of_start,
				 out_rmp.l_R0_peakday,
				 DAYS_IN_YEAR,
				 in_gmip.l_transmission_time_steps_per_day);
  if(update_flags.getFlag("l_R0_Amplitude"))
    { // STRICTLY A SCALAR QUANTITY
      if(in_umps[R0_AMP_INDEX].map_to_regional.design_matrix.nrows() > 0)
	{
	  ERROR_PARAM_INPUT(in_umps[R0_AMP_INDEX].param_name.c_str());
	}
      out_rmp.l_R0_Amplitude = gsl_vector_get(in_umps[R0_AMP_INDEX].param_value, 0) * out_rmp.l_R0_init / (1 + fn_initial_phase(in_gmip.l_day_of_start - out_rmp.l_R0_peakday, DAYS_IN_YEAR));
    }
  

}
