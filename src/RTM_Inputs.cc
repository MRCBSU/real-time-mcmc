#include "RTM_StructDefs.h" // provides RTM_Header.h
#include "RTM_FunctDefs.h"
#include "gsl_vec_ext.h"
#include "gsl_mat_ext.h"

using namespace std;
using std::string;

#define IN_FILE "rtm_input_files.txt"

#define ERROR_FILE_EXIT(err_string, ptr_filename) printf(err_string, ptr_filename);	\
  exit(2);

void input_filenames(
		     char *model_inputs,
		     char *model_parameters,
		    const int max_string_length)
{

  FILE *filenames;
  filenames = fopen(IN_FILE, "r");
  if(filenames == NULL)
    { // IF rtm_input_files DOES NOT EXIST, THEN SET *ALL* INPUTS TO DEFAULT VALUES. USER IS ALERTED.
      printf("Warning: rtm_input_files.txt not found. Default filenames to be searched for throughout.");
    }
  else
    {
      fclose(filenames);
      read_char_input(model_inputs, IN_FILE, "model_inputs:", max_string_length);
      read_char_input(model_parameters, IN_FILE, "model_parameters:", max_string_length);
    }
}

// SWAPS DEFAULT VARIABLE DEFAULT VALUES (VARIABLE NAMES
// READ IN AS A DELIMITED STRING) WITH VALUES READ IN FROM FILE
void read_variable_value(const string str_varnames, string& str_var_values, char* str_filename)
{

  // IS THE REQUIRED FILE PRESENT? IF NOT, AUTOMATICALLY USE DEFAULT VALUES
  FILE* tempfile = fopen(str_filename, "r");
  if(tempfile == NULL)
    return;
  else fclose(tempfile);

  // FIND NUMBER OF INSTANCE OF VARIABLE NAME DELIMITER IN str_varnames
  int num_instances = count_instances_in_string(str_varnames, ":");

  // CHECK THERE ARE THE SAME NUMBER OF INSTANCES OF THE DELIMITER IN str_var_values
  if(num_instances != count_instances_in_string(str_var_values, ":"))
    {
      ERROR_FILE_EXIT("Variable names and default values don't match for file %s\n", str_filename);
    }
  else if(num_instances == 0)
    {
      ERROR_FILE_EXIT("No delimiter detected when searching for variable to find in %s\n", str_filename);
    }

  // ITERATE THROUGH THESE INSTANCES
  string str_replace;
  register int i_var_start = 0, i_var_end = 0, i_var_mark;
  register int i_val_start = 0, i_val_end, i_val_mark;
  register double dbl_var;
  register string str_var;

  for(;
      i_var_end != string::npos;
      )
    {

      i_var_end = str_varnames.find(":", i_var_start);
      i_val_end = str_var_values.find(":", i_val_start);

      i_var_mark = (i_var_end == string::npos) ? str_varnames.length() : i_var_end;
      i_val_mark = (i_val_end == string::npos) ? str_var_values.length() : i_val_end;

      // GET A POINTER TO A NULL-TERMINATED STRING FROM str_varnames
      // str_var NEEDS TO BE PRE-ALLOCATED
      str_var = TrimStr(str_varnames.substr(i_var_start, i_var_mark - i_var_start));

      // SEARCH THE FILE FOR A REPLACEMENT VALUE -
      // DO NOTHING IF ONE NOT FOUND, ONLY NEED TO INCREMEMNT i_val_start AND i_val_end
      if(count_instances_in_file(str_filename, str_var.c_str()) > 0){
	read_double_input(dbl_var, str_filename, str_var.c_str());
	str_replace = dbl_to_string(dbl_var);

	// PLACE THE NEW VALUE INTO THE APPROPRIATE PLACE IN str_var_values
	str_var_values.replace(i_val_start, i_val_mark - i_val_start, str_replace);

	i_val_start = i_val_start + 1 + str_replace.length();

      } else {

	i_val_start = i_val_mark + 1;
	str_replace.clear();

      }

      i_var_start = i_var_end + 1;

    }

}

// POPULATE THE STRUCTURE FOR THE FIXED GLOBAL PARAMETERS
// (STRUCTURE CURRENTLY NEEDS NO PRIOR MEMORY ALLOCATION)
#define READ_NEXT_VARIABLE_VALUE read_from_delim_string<int>(str_vardefaults, \
							    ":",	\
							    int_delim_position);
#define READ_NEXT_DOUBLE_VARIABLE_VALUE read_from_delim_string<double>(str_vardefaults, \
								   ":", \
								   int_delim_position);



void read_mcmc_parameters(register mcmcPars &mcmc_pars,
			  char* source_file,
			  const string str_varnames,
			  string& str_vardefaults)
{

  int int_delim_position = 0;

  // READ IN FROM THE FILE ANY VALUES WHICH MAY CHANGE DUE TO USER INPUT
  read_variable_value(str_varnames, str_vardefaults, source_file);

  // VALUES WHICH CAN BE DETERMINED BY USER INPUT
  mcmc_pars.num_iterations = READ_NEXT_VARIABLE_VALUE;
  mcmc_pars.oType = (output_type) READ_NEXT_VARIABLE_VALUE;
  mcmc_pars.thin_output_every = READ_NEXT_VARIABLE_VALUE;
  mcmc_pars.thin_stats_every = READ_NEXT_VARIABLE_VALUE;
  mcmc_pars.adaptive_phase = READ_NEXT_VARIABLE_VALUE;
  mcmc_pars.adapt_every = READ_NEXT_VARIABLE_VALUE;
  mcmc_pars.burn_in = READ_NEXT_VARIABLE_VALUE;
  mcmc_pars.num_progress_reports = READ_NEXT_VARIABLE_VALUE;
  mcmc_pars.mixing_threshold_ub = READ_NEXT_DOUBLE_VARIABLE_VALUE;
  mcmc_pars.mixing_threshold_lb = READ_NEXT_DOUBLE_VARIABLE_VALUE;
  mcmc_pars.prop_var_inflation = READ_NEXT_DOUBLE_VARIABLE_VALUE;
  mcmc_pars.prop_var_shrinkage = READ_NEXT_DOUBLE_VARIABLE_VALUE;
  mcmc_pars.maximum_block_size = READ_NEXT_VARIABLE_VALUE;
  mcmc_pars.random_seed = READ_NEXT_VARIABLE_VALUE;
  mcmc_pars.max_num_threads = READ_NEXT_VARIABLE_VALUE;

  mcmcPars_alloc(mcmc_pars);

}

void read_global_fixed_parameters(register global_model_instance_parameters& fixed_pars,
				  char* source_file,
				  const string str_varnames,
				  string& str_vardefaults)
{

  int int_delim_position = 0;

  // READ IN FROM THE FILE ANY VALUES WHICH MAY CHANGE DUE TO USER INPUT
  read_variable_value(str_varnames, str_vardefaults, source_file);

  // VALUES WHICH CAN BE DETERMINED BY USER INPUT
  fixed_pars.l_gp_count_likelihood = (data_type) READ_NEXT_VARIABLE_VALUE;
  fixed_pars.l_tk = (transmission_kernel) READ_NEXT_VARIABLE_VALUE; // 0 FOR REED-FROST, 1 FOR MASS ACTION
  fixed_pars.l_transmission_time_steps_per_day = READ_NEXT_VARIABLE_VALUE;
  fixed_pars.l_reporting_time_steps_per_day = READ_NEXT_VARIABLE_VALUE;
  fixed_pars.l_duration_of_runs_in_days = READ_NEXT_VARIABLE_VALUE;
  fixed_pars.l_projection_days = READ_NEXT_VARIABLE_VALUE;
  fixed_pars.l_day_of_start = READ_NEXT_VARIABLE_VALUE;
  fixed_pars.l_start_day_of_week = READ_NEXT_VARIABLE_VALUE;
  fixed_pars.l_start_week = READ_NEXT_VARIABLE_VALUE;
  fixed_pars.l_num_regions = READ_NEXT_VARIABLE_VALUE;
  fixed_pars.l_GP_consultation_flag = READ_NEXT_VARIABLE_VALUE;
  fixed_pars.l_Hospitalisation_flag = READ_NEXT_VARIABLE_VALUE;
  fixed_pars.l_Deaths_flag = READ_NEXT_VARIABLE_VALUE;
  fixed_pars.l_Sero_data_flag = READ_NEXT_VARIABLE_VALUE;
  fixed_pars.l_Viro_data_flag = READ_NEXT_VARIABLE_VALUE;
  fixed_pars.l_GP_patch_flag = READ_NEXT_VARIABLE_VALUE;
  fixed_pars.l_Sero_delay = READ_NEXT_VARIABLE_VALUE;
  fixed_pars.l_GP_likelihood.lower = READ_NEXT_VARIABLE_VALUE;
  fixed_pars.l_GP_likelihood.upper = READ_NEXT_VARIABLE_VALUE;
  if(fixed_pars.l_GP_likelihood.upper > fixed_pars.l_duration_of_runs_in_days)
    printf("GP likelihood interval greater than number of days\n");
  fixed_pars.l_Hosp_likelihood.lower = READ_NEXT_VARIABLE_VALUE;
  fixed_pars.l_Hosp_likelihood.upper = READ_NEXT_VARIABLE_VALUE;
  if(fixed_pars.l_Hosp_likelihood.upper > fixed_pars.l_duration_of_runs_in_days)
    printf("Hospitalisations likelihood interval greater than number of days\n");
  fixed_pars.l_Deaths_likelihood.lower = READ_NEXT_VARIABLE_VALUE;
  fixed_pars.l_Deaths_likelihood.upper = READ_NEXT_VARIABLE_VALUE;
  if(fixed_pars.l_Deaths_likelihood.upper > fixed_pars.l_duration_of_runs_in_days)
    printf("Deaths likelihood interval greater than number of days\n");
  fixed_pars.l_Sero_likelihood.lower = READ_NEXT_VARIABLE_VALUE;
  fixed_pars.l_Sero_likelihood.upper = READ_NEXT_VARIABLE_VALUE;
  if(fixed_pars.l_Sero_likelihood.upper > fixed_pars.l_duration_of_runs_in_days)
    printf("Serological data likelihood interval greater than number of days\n");
  fixed_pars.l_Viro_likelihood.lower = READ_NEXT_VARIABLE_VALUE;
  fixed_pars.l_Viro_likelihood.upper = READ_NEXT_VARIABLE_VALUE;

  // ALLOC MEMORY FOR ARRAYS WITHIN THE fixed.pars STRUCTURE
  alloc_global_model_instance(fixed_pars);

  /// ASSIGN CALCULATED VALUES WITHIN THE STRUCTURE
  // d_week_numbers_by_day
  R_gl_fac_sublevel(fixed_pars.d_week_numbers_by_day,
		    fixed_pars.l_start_week,
		    fixed_pars.l_start_day_of_week,
		    7,
		    fixed_pars.l_duration_of_runs_in_days);

  //  if(fixed_pars.l_Viro_likelihood.upper > gsl_vector_int_max(fixed_pars.d_week_numbers_by_day) && fixed_pars.l_Viro_likelihood.lower < gsl_vector_int_min(fixed_pars.d_week_numbers_by_day)) // weekly data
  if(fixed_pars.l_Viro_likelihood.upper > fixed_pars.l_duration_of_runs_in_days && fixed_pars.l_Viro_likelihood.upper < fixed_pars.l_Viro_likelihood.lower)  // daily data
    printf("Invalid bounds for weekly data to be included in the virological data likelihood\n");

}

/*
 * Read a form of `<key> = <value>;` from `var_string`. Put `<value>` into `out_string`.
 * Return true on success or false on failure.
 */
bool read_string_from_instruct(string& out_string, const string& var_name, const string& var_string)
{  
  DEBUG(DEBUG_ALL, "Reading " << var_name << " from " << var_string);
  int var_dec_start = var_string.find(var_name);
  if(var_dec_start == string::npos) {
    return false;
  }
  int var_dec_end = var_string.find(";", ++var_dec_start);
  var_dec_start = var_string.find("=", var_dec_start);
  if(++var_dec_start < var_dec_end) {
     out_string = TrimStr(var_string.substr(var_dec_start, var_dec_end - var_dec_start));
     return true;
  }
  return false;
}

double read_double(const string var_name, const string var_string, const double default_value = 0.0)
{
  string var_substring;
  if (read_string_from_instruct(var_substring, var_name, var_string)) {
    DEBUG(DEBUG_ALL, "No double found for " << var_name << " returning default value of" << default_value)
    return atof(var_substring.c_str());
  } else {
    return default_value;
  }
}

int read_int(const string var_name, const string var_string, const int default_value = 0)
{
  string var_substring;
  if (read_string_from_instruct(var_substring, var_name, var_string)) {
    return atoi(var_substring.c_str());
  } else {
    DEBUG(DEBUG_ALL, "No int found for " << var_name << " returning default value of" << default_value)
    return default_value;
  }
}

void read_gsl_vector(gsl_vector* out_gslvec, const string var_name, const string var_string)
{
  string var_var_string;
  if (!read_string_from_instruct(var_var_string, var_name, var_string)) {
    DEBUG(DEBUG_ERROR, "No value found for " << var_name << " in " << var_string);
    exit(2);
  }
  gsl_vector_sscanf(var_var_string, out_gslvec);
}

void read_gsl_vector_int(gsl_vector_int* out_gslvec_int, const string var_name, const string var_string)
{
  string var_var_string;
  if (!read_string_from_instruct(var_var_string, var_name, var_string)) {
    DEBUG(DEBUG_ERROR, "No value found for " << var_name << " in " << var_string);
    exit(2);
  }
  gsl_vector_int_sscanf(var_var_string, out_gslvec_int);
}

void read_gsl_matrix(gsl_matrix* out_gslmat, const string var_name, const string var_string)
{
  string var_var_string;
  if (!read_string_from_instruct(var_var_string, var_name, var_string)) {
    DEBUG(DEBUG_ERROR, "Could not read " << var_name << " from string " << var_string);
    exit(2);
  }
  gsl_matrix_sscanf(var_var_string, out_gslmat);
}

void read_gsl_matrix_int(gsl_matrix_int* out_gslmat_int, const string var_name, const string var_string)
{
  string var_var_string;
  if (!read_string_from_instruct(var_var_string, var_name, var_string)) {
    DEBUG(DEBUG_ERROR, "Could not read " << var_name << " from string " << var_string);
    exit(2);
  }
  gsl_matrix_int_sscanf(var_var_string, out_gslmat_int);
}

void read_string_array(string* out_strarray, const int n, const string var_name, const string var_string)
{
  string var_var_string;
  if (!read_string_from_instruct(var_var_string, var_name, var_string)) {
    DEBUG(DEBUG_ERROR, "Could not read " << var_name << " from string " << var_string);
    exit(2);
  }
  int indx = 0;
  for(int i = 0; i < n; i++)
    out_strarray[i].assign(read_from_delim_string<string>(var_var_string, ",;", indx));
}

int flag_int_or_string(const string var_name, const string var_string, int default_return = 0)
{ // INDICATES WHETHER A NAME OR A VECTOR OF VALUES IS TO BE READ IN FROM THE INPUT var_string
  string var_var_string;
  if (!read_string_from_instruct(var_var_string, var_name, var_string)) {
    // VARIABLE NAME DOESN'T APPEAR IN INPUT STRING, RETURN DEFAULT VALUE
    return default_return;
  } else if(var_var_string.find_first_of("1234567890-", 0) > 1) {
    // IT'S A STRING
    return 1;
  } else {
    // IT'S A NUMERICAL VECTOR
    return 0;
  }
}

/////////// GET A PARAMETER SUBSTRING FROM A FILE STRING
int cut_out_parameter(string& out_substring, const string file_content, const char* param_name)
{

  // TAKES A SUBSTRING OF FILE STRING CONTENTS OF ALL INFORMATION RELEVANT TO PARAMETER param_name
  // RETURNS A POSITIVE INTEGER IF STRING param_name IS FOUND IN file_content, ZERO OTHERWISE
  int num_strings = count_instances_in_string(file_content, param_name);
  if(num_strings > 0)
    {
    int substring_position = nth_instance_in_string(file_content, param_name, num_strings);

    /// FIND THE START OF THE VARIABLE DECLARATION
    int var_dec_start = file_content.find("{", substring_position);
    int var_dec_end = file_content.find("}", var_dec_start);
    if(var_dec_start == string::npos && var_dec_end == string::npos)
      { // ERROR
	printf("Error defining variable %s\n", param_name);
	exit(2);
      }
    out_substring = file_content.substr(var_dec_start + 1, var_dec_end - var_dec_start);

    }

  return num_strings;
}

void alloc_and_set_breakpoint_vector(gsl_vector_int*& vector_of_breakpoints,
				     int& parameter_dimension,
				     const string str_source,
				     const string str_property_name,
				     const int pos_in_source,
				     const int default_size,
				     const bool flag_alternate_variation)
{

  // set to default value
  string temp_string("false");

  // read in the variable set value if one exists
  if(pos_in_source != string::npos)
      read_string_from_instruct(temp_string, str_property_name, str_source);

  if(temp_string.compare("false") != 0 && ((temp_string.compare("true") == 0) || (flag_alternate_variation)))
    {
      // if there are breakpoints, then we should have a true value or a numeric vector
      if(temp_string.compare("true") == 0){
	// FULL VARIATION IN THE PARAMETER
	vector_of_breakpoints = gsl_vector_int_alloc(default_size - 1);
	gsl_vector_int_set_seq(vector_of_breakpoints);
      } else if(flag_alternate_variation) {
	// THERE SHOULD BE A VECTOR OF VALUES TO READ IN
	// FIRST, GET THE NUMBER OF DISTINCT TIME INTERVALS
	// AND ALLOCATE MEMORY TO THE VECTOR
	vector_of_breakpoints = gsl_vector_int_alloc(count_delims_in_string(temp_string, ",") + 1);
	// ASSIGN THE READ VALUES INTO THE VECTOR
	gsl_vector_int_sscanf(temp_string, vector_of_breakpoints);
      }
      parameter_dimension = vector_of_breakpoints->size + 1;
    }
  else 
    {
      vector_of_breakpoints = 0;
      parameter_dimension = 1;
    }
}


void read_param_regression(regression_def& reg_def,
			   const int param_dimension,
			   const string str_source,
			   const int num_regions,
			   const int num_days,
			   const int num_ages)
{

  int dim_r = 1, dim_t = 1, dim_a = 1;
  // IF THERE IS NO SPATIAL, REGIONAL OR TEMPORAL VARIATION THEN THIS CAN BE "SET TO NULL"

  // A MODEL EXISTS IF ANY OF regional_breakpoints, time_breakpoints or age_breakpoints ARE SPECIFIED
  int reg_pos = str_source.find("region_breakpoints");
  int time_pos = str_source.find("time_breakpoints");
  int age_pos = str_source.find("age_breakpoints");

  reg_def.region_breakpoints = reg_def.time_breakpoints = reg_def.age_breakpoints = 0;
  reg_def.design_matrix = 0;

  if(((age_pos != string::npos) || (time_pos != string::npos)) || (reg_pos != string::npos))
    {
      // WE HAVE SOME FORM OF VARIATION
      int link_pos = str_source.find("regression_link");
      if(link_pos != string::npos)
	reg_def.regression_link = (link_function) read_int("regression_link", str_source);
      else reg_def.regression_link = cIDENTITY;

      alloc_and_set_breakpoint_vector(reg_def.region_breakpoints, dim_r, str_source, "region_breakpoints", reg_pos, num_regions, false);

      alloc_and_set_breakpoint_vector(reg_def.time_breakpoints, dim_t, str_source, "time_breakpoints", time_pos, num_days, true);

      alloc_and_set_breakpoint_vector(reg_def.age_breakpoints, dim_a, str_source, "age_breakpoints", age_pos, num_ages, true);


      // THE DESIGN MATRIX...
      // SHOULD HAVE NUMBER OF COLUMNS EQUAL TO THE PARAMETER DIMENSION
      // SHOULD HAVE NUMBER OR ROWS EQUAL TO THE TOTAL NUMBER OF (region, time, age) GROUP COMBINATIONS
      reg_def.design_matrix = gsl_matrix_alloc(dim_r * dim_t * dim_a, param_dimension);

      // TO ALLOCATE THE MATRIX, EITHER A FILENAME SHOULD BE GIVEN OR A VECTOR OF NUMBERS
      bool flag_infile = (bool) flag_int_or_string("regression_design", str_source);
      if(flag_infile){
        string str_filename;
        if (!read_string_from_instruct(str_filename, "regression_design", str_source)) {
          DEBUG(DEBUG_ERROR, "Could not read regression_design from string " << str_source);
          exit(2);
        }
        FILE* design_file = fopen(str_filename.c_str(), "r");

        gsl_matrix_fscanf(design_file, reg_def.design_matrix);

        fclose(design_file);

      } else if(str_source.find("regression_design") != string::npos) // READ IN THE MATRIX AS A VECTOR
	read_gsl_matrix(reg_def.design_matrix, "regression_design", str_source);

      else // DEFAULT, SET THE MATRIX TO THE IDENTITY (EVEN IF A RECTANGULAR MATRIX)
	gsl_matrix_set_identity(reg_def.design_matrix);

    }


}



/////////// READ IN A PARAMETER STRUCTURE
void read_modpar(updateable_model_parameter& modpar,
		 const char* param_name,
		 const string in_string,
		 const double default_value,
		 const string likelihood_flags,
		 const int num_regions,
		 const int num_days,
		 const int num_ages)
{
  DEBUG(DEBUG_DETAIL, "Reading parameters for " << param_name);
  string var_string, var_var_string;
  int int_marker = 0;
  int var_dimension;

  /// ASSIGN PARAMETER NAME
  modpar.param_name.assign(param_name);

  /// DOES THE PARAMETER STRING OCCUR IN THE INPUT FILE? IF SO, PUT INTO VAR_STRING THE RELEVANT PART OF THE FILE
  int num_strings = cut_out_parameter(var_string, in_string, param_name);

  if(num_strings > 0){ // THE STRING IS FOUND
    
    /// GET PARAMETER DIMENSION
    if (!read_string_from_instruct(var_var_string, "param_value", var_string)) {
      DEBUG(DEBUG_ERROR, "No value found for param_value in " << var_string);
      exit(2);
    }
    /// NUMBER OF VALUE DELIMITERS (COMMAS) + 1 SHOULD GIVE THE DIMENSION
    var_dimension = count_delims_in_string(var_var_string, ",") + 1;
    /// 

    /// CAN NOW ASSIGN SOME OF THE VECTORS IN THE UPDATEABLE_PARAMETER STRUCTURE
    modpar.param_value = gsl_vector_alloc(var_dimension);
    modpar.prior_distribution = gsl_vector_int_alloc(var_dimension);
    /// ///

    /// ASSIGN SOME OF THESE VALUES ///
    // parameter values
    gsl_vector_sscanf(var_var_string, modpar.param_value);
    read_gsl_vector_int(modpar.prior_distribution, "prior_distribution", var_string);

    // IF ONE COMPONENT OF prior_distribution IS A MULTIVARIATE DISTRIBUTION THEN ALL MUST BE
    if(in_gsl_vector_int(modpar.prior_distribution, (int) cMVNORMAL))
      gsl_vector_int_set_all(modpar.prior_distribution, (int) cMVNORMAL);

    /// IS THIS NODE TO BE UPDATED IN THIS MCMC SIMULATION?
    modpar.flag_update = (gsl_vector_int_max(modpar.prior_distribution) > (int) cCONSTANT);

    /// ALLOCATE AND ASSIGN ALL PROPERTIES WHOSE USE IS CONDITIONAL ON BEING A STOCHASTIC NODE IN THE MODEL
    if(modpar.flag_update)
      {
	int total_number_prior_params = 0;

	// ALLOCATE
	modpar.proposal_value = gsl_vector_alloc(var_dimension);
	modpar.proposal_variances = gsl_vector_alloc(var_dimension);
	modpar.posterior_mean = gsl_vector_alloc(var_dimension);
	gsl_vector_set_zero(modpar.posterior_mean);
	modpar.posterior_sumsq = gsl_vector_alloc(var_dimension);
	gsl_vector_set_zero(modpar.posterior_sumsq);
	modpar.prior_params = (gsl_vector **) calloc(var_dimension, sizeof(gsl_vector *));

	// DETERMINE WHETHER OR NOT THE PARAMETER HAS A HYPERPRIOR (it will have the property prior_hyper set to a value > 0 if so
  string tempstring;
	modpar.flag_hyperprior = read_string_from_instruct(tempstring, "prior_hyper", var_string);
  modpar.flag_hyperprior = modpar.flag_hyperprior && atoi(tempstring.c_str()) > 0;

	// ERROR IF HYPERPRIOR IS SPECIFIED FOR A MULTIVARIATE PARAMETER
	if(modpar.flag_hyperprior && (gsl_vector_int_get(modpar.prior_distribution, 0) == (int) cMVNORMAL)){
	  DEBUG(DEBUG_ERROR, "Multivariate normal selected for the hyperprior " << param_name);
	  exit(2);
	}

	// COMPONENTWISE ALLOCATATION OF MEMORY FOR THE PRIORS... HYPERPRIORS SHOULD HAVE THEIR PARAMETERS SET ELSEWHERE
	if(!modpar.flag_hyperprior){
    DEBUG(DEBUG_ALL, "Has a hyperprior")

	  for(int_marker = 0; int_marker < var_dimension;){
	    int nparam = num_parameters_by_distribution((distribution_type) gsl_vector_int_get(modpar.prior_distribution, int_marker), var_dimension);
	    modpar.prior_params[int_marker++] = (nparam > 0) ? gsl_vector_alloc(nparam) : 0;
	    total_number_prior_params += nparam; // modpar.prior_params[int_marker++]->size;
	  }

	}

	// ASSIGN PROPOSAL VARIANCES
	read_gsl_vector(modpar.proposal_variances,
			"proposal_variance",
			var_string);

	// CHECK THAT WE'RE NOT DEALING WITH A PARAMETER WHICH HAS A HYPERPRIOR - PRIOR LINKS WILL BE ESTABLISHED ELSEWHERE
	if(!modpar.flag_hyperprior){
	  gsl_vector* temp_vec = gsl_vector_alloc(total_number_prior_params);


	  read_gsl_vector(temp_vec,
			  "prior_parameters",
			  var_string);

	  int offset = 0;
	  for(int_marker = 0; int_marker < var_dimension; int_marker++){
	    if(modpar.prior_params[int_marker] != 0){
	      gsl_vector_subvector_memcpy(modpar.prior_params[int_marker],
					  temp_vec,
					  offset);
	      offset += modpar.prior_params[int_marker]->size;
	    }
	  }

	  gsl_vector_free(temp_vec);

	}

	modpar.number_accepted_moves = gsl_vector_int_alloc(var_dimension);
	modpar.number_proposed_moves = gsl_vector_int_alloc(var_dimension);
	gsl_vector_int_set_zero(modpar.number_accepted_moves);
	gsl_vector_int_set_zero(modpar.number_proposed_moves);

	// LIKELIHOOD FLAGS SHOULD BE DETERMINED IN RTM_StructAssign.h
	modpar.flag_transmission_model = atoi((likelihood_flags.substr(0, 1)).c_str());
	modpar.flag_reporting_model = atoi((likelihood_flags.substr(1, 1)).c_str());
	modpar.flag_GP_likelihood = atoi((likelihood_flags.substr(2, 1)).c_str());
	modpar.flag_Hosp_likelihood = atoi((likelihood_flags.substr(3, 1)).c_str());
	modpar.flag_Sero_likelihood = atoi((likelihood_flags.substr(4, 1)).c_str());
	modpar.flag_Viro_likelihood = atoi((likelihood_flags.substr(5, 1)).c_str());
      }


  
  } else { // THE STRING IS NOT FOUND
    DEBUG(DEBUG_WARNING, "No parameter structure found.")

    /// SET EQUAL TO THE DEFAULT (TIME AND AGE INVARIANT VALUE)
    modpar.param_value = gsl_vector_alloc(1);
    modpar.prior_distribution = gsl_vector_int_alloc(1);

    gsl_vector_set(modpar.param_value, 0, default_value);
    gsl_vector_int_set(modpar.prior_distribution, 0, 1);

  }

  // ESTABLISH THE MAPS FROM THE PARAMETER VALUES TO THE VALUES PASSED TO THE REGIONAL STRUCTURES
  read_param_regression(modpar.map_to_regional, modpar.param_value->size, var_string, num_regions, num_days, num_ages);

}


/////// INITIALISE PRIOR DENSITIES
void initialise_prior_density(updateable_model_parameter& modpar,
			      double (*prior_density_fn)(const double&, const distribution_type&, const gsl_vector*))
{

  gsl_vector_int* prior_distributions = gsl_vector_int_alloc(modpar.param_value->size);

  if(modpar.flag_hyperprior)
    gsl_vector_int_set_all(prior_distributions, gsl_vector_int_get(modpar.prior_distribution, 0));
  else
    gsl_vector_int_memcpy(prior_distributions, modpar.prior_distribution);

  modpar.log_prior_dens = 0;

  // CHECK THE PRIOR IS NOT MULTIVARIATE
  if(gsl_vector_int_get(prior_distributions, 0) != (int) cMVNORMAL)
    {
      for(int int_marker = 0; int_marker < modpar.param_value->size; int_marker++)
	if(modpar.prior_params[int_marker] != 0)
	  modpar.log_prior_dens += prior_density_fn(gsl_vector_get(modpar.param_value, int_marker),
						    (distribution_type) gsl_vector_int_get(prior_distributions, int_marker),
						    modpar.prior_params[int_marker]);
    }
  else
    {
      // SET UP THE MULTIVARIATE NORMAL MEAN AND COVARIANCE STRUCTURES
      gsl_vector* mv_mean = gsl_vector_alloc(modpar.param_value->size);
      gsl_matrix* mv_covar = gsl_matrix_alloc(modpar.param_value->size, modpar.param_value->size);

      for(int int_marker = 0; int_marker < modpar.param_value->size; int_marker++)
	{
	  gsl_vector_set(mv_mean, int_marker, gsl_vector_get(modpar.prior_params[int_marker], 0));
	  gsl_vector_view mv_covar_subrow = gsl_vector_subvector(modpar.prior_params[int_marker], 1, modpar.prior_params[int_marker]->size - 1);
	  gsl_matrix_set_row(mv_covar, int_marker, &mv_covar_subrow.vector);
	}

      // TRANSFER THIS DATA INTO THE MVNORM CLASS HELD IN THE UPDATEABLE_MODEL_PARAMETERS STRUCTURE
      modpar.prior_multivariate_norm = new mvnorm(mv_mean, mv_covar);

      // EVALUATE THE LOG-DENSITY
      if(modpar.flag_hyperprior)
	modpar.log_prior_dens = (*modpar.prior_multivariate_norm).ld_mvnorm(modpar.param_value);
      else
	modpar.log_prior_dens = (*modpar.prior_multivariate_norm).ld_mvnorm_nonnorm(modpar.param_value);

      // FREE ANY MEMORY ALLOCATED
      gsl_vector_free(mv_mean);
      gsl_matrix_free(mv_covar);
    }
  gsl_vector_int_free(prior_distributions);

}


void node_links(globalModelParams& in_pars,
		const string str_param_input)
{

  register int inti, intj;
  int int_delim_position = 0;
  string temp_param, temp_prior_param;

  // FOR EACH PARAMETER, i
  for(inti = 0; inti < in_pars.size_param_list; inti++)
    {

      in_pars.param_list[inti].flag_any_child_nodes = false;

      for(intj = 0; intj < in_pars.size_param_list; intj++)
	{
	  //FOR ALL OTHER PARAMETERS, j
	  if(inti != intj){

	    // NO PARAMETERS ARE HYPERPARAMETERS BY DEFAULT, SO NEEDS TO BE SPECIFIED IN FILE.
	    // CHECK THE VARIABLE APPEARS IN THE FILE
	    // STORE THE VARIABLE STRING
	    int num_strings = cut_out_parameter(temp_param, str_param_input, in_pars.param_list[intj].param_name.c_str());

	    if(num_strings > 0)
	      {

		// FIND THE PRIOR_PARAMS VALUE, they are expected to be missing sometimes
		(void)!read_string_from_instruct(temp_prior_param, "prior_parameters", temp_param);

		// IS THE TRIMMED STRING EQUAL TO THE NAME OF PARAMETER i
		// TRIM FIRST AND LAST CHARACTERS (i.e. QUOTATION MARKS)
		temp_prior_param.erase(0, 1);
		temp_prior_param.erase(FN_MAX(temp_prior_param.length(), 1) - 1, 1);

		if(temp_prior_param.compare(in_pars.param_list[inti].param_name) == 0){
		  // PARAMETER j IS A CHILD NODE OF PARAMETER i
		  // SET PRIOR_PARAMS OF j TO BE A POINTER TO PARAM_VALUES OF i
		  for(int intk = 0; intk < in_pars.param_list[intj].param_value->size; intk++){
		    if(gsl_vector_int_get(in_pars.param_list[intj].prior_distribution, 1) != (int) cCONSTANT)
		      in_pars.param_list[intj].prior_params[intk] = in_pars.param_list[inti].param_value;
		    

		  }
		  in_pars.param_list[inti].flag_child_nodes[intj] = true;
		  in_pars.param_list[inti].flag_any_child_nodes = true;
		}
	      }
	  }
	}
    }


}

void aggregate_delays(infection_to_data_delay& itdd)
{
  itdd.overall_delay.gamma_mean = itdd.overall_delay.gamma_sd = 0.0;
  for(int inti = 0; inti < itdd.num_components; inti++){
    itdd.overall_delay.gamma_mean += itdd.component_delays[inti].gamma_mean;
    itdd.overall_delay.gamma_sd += gsl_pow_2(itdd.component_delays[inti].gamma_sd); // FOR THE MOMENT, THE GAMMA_SD IS STORING THE AGGREGATE VARIANCE
  }
  itdd.overall_delay.gamma_shape = gsl_pow_2(itdd.overall_delay.gamma_mean) / itdd.overall_delay.gamma_sd;
  itdd.overall_delay.gamma_rate = itdd.overall_delay.gamma_mean / itdd.overall_delay.gamma_sd;
  itdd.overall_delay.gamma_scale = 1.0 / itdd.overall_delay.gamma_rate;
  itdd.overall_delay.gamma_sd = sqrt(itdd.overall_delay.gamma_sd); // NOW WE'VE USED THE VARIANCE, TRANSFORM BACK TO THE STANDARD DEVIATION

}

void calculate_fixed_distribution_function(infection_to_data_delay& itdd, double delta_t = 1)
{

  // FIND THE MAXIMAL INTEGER NECESSARY TO ADEQUATELY TRUNCATE THE DELAY DISTRIBUTIONS
  double dbl_k = gsl_cdf_gamma_Pinv(1 - EPSILON_TOLERANCE, itdd.overall_delay.gamma_shape, itdd.overall_delay.gamma_scale); // MAXIMAL NUMBER OF DAYS TO CONSIDER
  int int_k0 = (int) ceil(dbl_k / delta_t);
  double cumulative_probability = 0.0;

  itdd.distribution_function = gsl_vector_alloc(++int_k0); // INCREMENT int_k0 TO ENSURE THAT THE FINAL PROBABILITY IS LESS THAN EPSILON_TOLERANCE

  for(int int_i = 0; int_i < (int_k0 - 1); int_i++)
    {

      gsl_vector_set(itdd.distribution_function,
		     int_i,
		     gsl_cdf_gamma_P(((double) int_i + 1) * delta_t, itdd.overall_delay.gamma_shape, itdd.overall_delay.gamma_scale) - cumulative_probability);

      cumulative_probability += gsl_vector_get(itdd.distribution_function, int_i);

    }

  gsl_vector_set(itdd.distribution_function,
		 int_k0 - 1,
		 1.0 - cumulative_probability); // TRUNCATE THE DISTRIBUTION AT GREATER THAN (int_k0 - 1) DAYS

}


void load_delays(st_delay& out_delay, const string str_name, const double mean, const double sd, const string input_string)
{

  string var_string;

  out_delay.delay_name.assign(str_name);

  /// ARE VALUES FOR THIS PARAMETER PRESENT IN THE INPUT FILE
  int num_strings = cut_out_parameter(var_string, input_string, str_name.c_str());

  if(num_strings > 0){ // THE STRING IS FOUND

    out_delay.gamma_mean = read_double("gamma_mean", var_string, mean);

    out_delay.gamma_sd = read_double("gamma_sd", var_string, sd);

  }

}



  /////////// ALLOCATE MEMORY FOR ALL THE PARAMETER STRUCTURES
  /////////// READ IN ALL THE PARAMETER STRUCTURES.
  /////////// ASSIGNS INITIAL VALUES

void read_global_model_parameters(globalModelParams& in_pars,
				  const char* source_file,
				  const string str_var_names,
				  const string& str_var_initvals,
				  const string& str_var_likflags,
				  const string& str_delay_names,
				  const string& str_delay_means,
				  const string& str_delay_sds,
				  const string& str_delay_flags,
				  const int num_regions,
				  const int num_days,
				  const int num_ages,
				  const double reporting_time_steps_per_day)
{

  string str_source;
  int int_delim_name_position = 0, int_delim_val_position = 0, int_delim_flag_position = 0;
  int inti;

  // ALLOCATE THE STRUCTURE... NEED TO KNOW HOW MANY VARIABLES IT CONTAINS
  int num_instances = count_instances_in_string(str_var_names, ":");

  if((num_instances != count_instances_in_string(str_var_initvals, ":")) || 
     (num_instances != count_instances_in_string(str_var_likflags, ":")))
    {
      printf("Unequal length of input strings to fn: read_global_model_parameters\n");
      exit(2);
    }

  // ALLOCATE THE NECESSARY MEMORY
  globalModelParams_alloc(in_pars, ++num_instances);

  // LOAD THE INPUT FILE INTO A STRING VARIABLE
  fn_load_file(&str_source, source_file);

  // READ IN PARAMETER STRUCTURES OR SET TO DEFAULT.
  for(inti = 0; inti < num_instances; inti++){

    // GET VARIABLE NAME FROM DEFAULT STRING
    string str_param_name = read_from_delim_string<string>(str_var_names, ":", int_delim_name_position);

    // GET VARIABLE DEFAULT FROM DEFAULT STRING
    double dbl_param_val = read_from_delim_string<double>(str_var_initvals, ":", int_delim_val_position);

    // GET VARIABLE LIKELIHOOD INDICATORS FROM DEFAULT STRING
    string str_param_flags = read_from_delim_string<string>(str_var_likflags, ":", int_delim_flag_position);

    // POPULATE THE STRUCTURE
    read_modpar(in_pars.param_list[inti], str_param_name.c_str(), str_source, dbl_param_val, str_param_flags, num_regions, num_days, num_ages);

  }

  // LINK PARENT AND CHILD NODES
  node_links(in_pars, str_source);

  for(inti = 0; inti < num_instances; inti++)
    // CALCULATE INITIAL PRIOR DENSITY - BUT FIRST CHECK THAT THE NODE IS STOCHASTIC
    if(in_pars.param_list[inti].flag_update)
      {
	if(in_pars.param_list[inti].flag_hyperprior)
	  initialise_prior_density(in_pars.param_list[inti], R_univariate_prior_log_density);
	else
	  initialise_prior_density(in_pars.param_list[inti], R_univariate_prior_log_density_nonnorm);
      }


  // NOW TO READ IN THE PARAMETERS OF THE DELAY DISTRIBUTIONS
  // ALLOCATE THE STRUCTURE... NEED TO KNOW HOW MANY VARIABLES IT CONTAINS
  num_instances = count_instances_in_string(str_delay_names, ":");

  if((num_instances != count_instances_in_string(str_delay_means, ":")) || 
     (num_instances != count_instances_in_string(str_delay_sds, ":")) ||
     (num_instances != count_instances_in_string(str_delay_flags, ":")))
    {
      printf("Unequal length of default delay strings to fn: read_global_model_parameters\n");
      exit(2);
    }

  // ALLOCATE THE NECESSARY MEMORY

  // LOOPING
  int int_delim_mean_position = 0, int_delim_sd_position = 0;
  int_delim_name_position = int_delim_flag_position = 0;
  st_delay* temp_delay = new st_delay[++num_instances];
  string* str_param_flags = new string[num_instances];

  int gp_delay_counter = 0, hosp_delay_counter = 0, death_delay_counter = 0;

  for(inti = 0; inti < num_instances; inti++)
    {

      // get the name
      string str_param_name = read_from_delim_string<string>(str_delay_names, ":", int_delim_name_position);

      // get the gamma mean
      double dbl_param_mean = read_from_delim_string<double>(str_delay_means, ":", int_delim_mean_position);

      // get the gamma sd
      double dbl_param_sd = read_from_delim_string<double>(str_delay_sds, ":", int_delim_sd_position);

      // get the delay flags
      str_param_flags[inti] = read_from_delim_string<string>(str_delay_flags, ":", int_delim_flag_position);

      // load the individual delays into an element of temp_delay
      load_delays(temp_delay[inti], str_param_name, dbl_param_mean, dbl_param_sd, str_source);

      gp_delay_counter += atoi((str_param_flags[inti].substr(0, 1)).c_str());
      hosp_delay_counter += atoi((str_param_flags[inti].substr(1, 1)).c_str());
      death_delay_counter += atoi((str_param_flags[inti].substr(2, 1)).c_str());


    }

  // ALLOCATE THE NECESSARY MEMORY TO THE AGGREGATED DELAY STRUCTURE
  infection_to_data_delay_alloc(in_pars.gp_delay, gp_delay_counter, 0);
  infection_to_data_delay_alloc(in_pars.hosp_delay, hosp_delay_counter, 0);
  infection_to_data_delay_alloc(in_pars.death_delay, death_delay_counter, 0);

  // COPY THE ELEMENTS OF temp_delay INTO THE AGGREGATED DELAY STRUCTUES
  for(inti = gp_delay_counter = hosp_delay_counter = death_delay_counter = 0;
      inti < num_instances;
      inti++)
    {
      if((bool) atoi((str_param_flags[inti].substr(0, 1)).c_str())){
	st_delay_memcpy(in_pars.gp_delay.component_delays[gp_delay_counter], temp_delay[inti]);
	gp_delay_counter++;
      }
      if((bool) atoi((str_param_flags[inti].substr(1, 1)).c_str())){
	st_delay_memcpy(in_pars.hosp_delay.component_delays[hosp_delay_counter], temp_delay[inti]);
	hosp_delay_counter++;
      }
      if((bool) atoi((str_param_flags[inti].substr(2, 1)).c_str())){
	st_delay_memcpy(in_pars.death_delay.component_delays[death_delay_counter], temp_delay[inti]);
	death_delay_counter++;
      }

    }

  // DELETE MEMORY ASSIGNED TO TEMPORARY ARRAYS
  delete [] temp_delay;
  delete [] str_param_flags;

  // FOR EACH OF THE DIFFERENT DELAYS, CALCULATE THE AGGREGATED DELAY
  aggregate_delays(in_pars.gp_delay);
  aggregate_delays(in_pars.hosp_delay);
  aggregate_delays(in_pars.death_delay);

  // CALCULATE THE DISTRIBUTION FUNCTION FOR THE AGGREGATED DELAYS
  calculate_fixed_distribution_function(in_pars.gp_delay, 1.0 / reporting_time_steps_per_day);
  calculate_fixed_distribution_function(in_pars.hosp_delay, 1.0 / reporting_time_steps_per_day);
  calculate_fixed_distribution_function(in_pars.death_delay, 1.0 / reporting_time_steps_per_day);

}

int param_list_index(const globalModelParams src_gmp, string target_name)
{

  int int_i = 0;

  for(; int_i < src_gmp.size_param_list; int_i++)
    {
      if(src_gmp.param_list[int_i].param_name.compare(target_name) == 0)
	break;
    }

  return (int_i < src_gmp.size_param_list) ? int_i : -1;
}

/// READ IN THE DATA MATRICES AS FORMATTED BY THE HPA
void data_matrices_fscanf(const string& infilestring,
			  gsl_matrix* reported_data,
			  const int col_skip,
			  const int row_skip)
{ //NEEDS SOME ERROR HANDLING IN HERE.
  register int inti, intj;
  char str_row_date[25];
  double dbl_dummy;

  FILE *infile = fopen(infilestring.c_str(), "r");

  for(inti = 0; inti < (reported_data->size1 + row_skip); inti++){

    for(intj = 0; intj < (reported_data->size2 + col_skip); intj++){

      if((inti < row_skip) || (intj < col_skip))

	fscanf(infile, "%s", str_row_date);
    
      else
	
	fscanf(infile, "%lf", gsl_matrix_ptr(reported_data, inti - row_skip, intj - col_skip));
 
    }

  }

  fclose(infile);
}

/// READ IN THE DATA MATRICES AS FORMATTED BY THE HPA
void data_int_matrices_fscanf(string infilestring,
			  gsl_matrix_int* reported_data,
			  const int col_skip = 0,
			  const int row_skip = 0)
{ //NEEDS SOME ERROR HANDLING IN HERE.
  register int inti, intj;
  char str_row_date[25];
  double dbl_dummy;

FILE *infile = fopen(infilestring.c_str(), "r");

  for(inti = 0; inti < (reported_data->size1 + row_skip); inti++){

    for(intj = 0; intj < (reported_data->size2 + col_skip); intj++){

      if((inti < row_skip) || (intj < col_skip))

	fscanf(infile, "%s", str_row_date);
    
      else
	
	fscanf(infile, "%d", gsl_matrix_int_ptr(reported_data, inti - row_skip, intj - col_skip));
 
    }

  }

  fclose(infile);
}


// ERROR CHECKING WRAPPER FOR read_string_array
void fetch_filenames(string* out_string,
		     const int& num_strings,
		     const string var_name,
		     const string str_file)
{

  string var_filenames("");

  // find the array of filenames
  if(!var_name.empty())
    if (!read_string_from_instruct(var_filenames, var_name, str_file)) {
    // ERROR_FILE_EXIT("Required data structure %s no specified\n", var_name.c_str());
    for(int int_i = 0; int_i < num_strings;) // THIS CATERS FOR THE CASE OF NO MATCHING STRING.. WHAT IF THERE ARE INSUFFICIENT FILENAMES SPECIFIED? NEED TO CONSIDER THAT CASE.
      out_string[int_i++] = "";
  } else {
  read_string_array(out_string, num_strings, var_name, str_file);
  }

}

// READ THE FILENAMES TO A TEMP VECTOR
void read_filenames_and_data_matrices(gsl_matrix** data_matrices,
				      const int num_regions,
				      const string var_name,
				      const string var_string,
				      const int col_skip = 0,
				      const int row_skip = 0)
{

  string* str_vec_filenames = new string[num_regions];

  // find the array of filenames
  fetch_filenames(str_vec_filenames, num_regions, var_name, var_string);

  for(int int_i = 0; int_i < num_regions; int_i++)

    data_matrices_fscanf(str_vec_filenames[int_i], data_matrices[int_i], col_skip, row_skip);

  // delete or close any memory allocated
  delete [] str_vec_filenames;

}

// READ THE FILENAMES TO A TEMP VECTOR
void read_filenames_and_data_int_matrices(gsl_matrix_int** data_matrices,
				      const int num_regions,
				      const string var_name,
				      const string var_string,
				      const int col_skip = 0,
				      const int row_skip = 0)
{

  string* str_vec_filenames = new string[num_regions];

  // find the array of filenames
  fetch_filenames(str_vec_filenames, num_regions, var_name, var_string);

  for(int int_i = 0; int_i < num_regions; int_i++)

    data_int_matrices_fscanf(str_vec_filenames[int_i], data_matrices[int_i], col_skip, row_skip);

  // delete or close any memory allocated
  delete [] str_vec_filenames;

}


/// FUNCTIONS UNIQUE TO READING IN THE DATA-LIKE INPUTS
void read_metaregion_datatype(rtmData** obj_data, const gsl_matrix* population, string* countfiles, string* denomfiles, const int& num_regions, const string str_var_count, const string str_var_denom, const string str_var_agg, const string& str_source, const bool& normalise_flag){

  
  fetch_filenames(countfiles, num_regions, str_var_count, str_source);
  fetch_filenames(denomfiles, num_regions, str_var_denom, str_source);
  // Get the required level of data aggregation
  bool missing_data = (str_source.find(str_var_agg) != string::npos);

  for(int int_i = 0; int_i < num_regions; int_i++)
    {
      obj_data[int_i]->setDim(missing_data, str_source, str_var_agg, NUM_AGE_GROUPS);
      obj_data[int_i]->read(denomfiles[int_i], countfiles[int_i]);
      if(normalise_flag){
	gsl_vector_const_view population_row = gsl_matrix_const_row(population, int_i);
	obj_data[int_i]->normalise(&population_row.vector);
      } else if(missing_data){
	gsl_vector_const_view population_row = gsl_matrix_const_row(population, int_i);
	obj_data[int_i]->data_population_sizes(&population_row.vector);
      }
    }

}

void read_data_inputs(Region* meta_region, const string str_input_filename,
		      const int& num_regions)
{

  string* temp_string = new string[num_regions];
  string str_source, str_var;
  int int_i;

  // all memory used by meta_region is assumed pre-assigned EXCEPT
  // for the internal margins of the data objects which are allocated here.

  // load in the file containing the inputs
  fn_load_file(&str_source, str_input_filename.c_str());

  // cut-out the data structure
  int num_strings = cut_out_parameter(str_var, str_source, "study_region");

  if(num_strings == 0)
    {
      ERROR_FILE_EXIT("Data structure not found in file %s", str_input_filename.c_str()); 
    }
  
  // REGION NAMES
  // read in the string containing the region names
  read_string_array(temp_string, num_regions, "regions_used", str_var);

  for(int_i = 0; int_i < num_regions; int_i++)
    meta_region[int_i].name.assign(temp_string[int_i]);

  delete [] temp_string;

  // POPULATION
  gsl_matrix* tempmat = gsl_matrix_alloc(num_regions, NUM_AGE_GROUPS);

  read_gsl_matrix(tempmat, "regions_population", str_var);
  // TRANSFER FROM THE TEMPORARY MATRIX TO THE META_REGION STRUCTURE
  for(int_i = 0; int_i < num_regions; int_i++){
    gsl_vector_view tempmat_row = gsl_matrix_row(tempmat, int_i);
    gsl_vector_memcpy(meta_region[int_i].population, &tempmat_row.vector);
    meta_region[int_i].total_population = gsl_vector_sum_elements(meta_region[int_i].population);
  }

  // READ IN THE DATA..
  // CYCLING THROUGH DATA STREAMS, POPULATE THE DATA OBJECTS WITHIN EACH REGION.
  string* countfiles = new string [num_regions];
  string* denomfiles = new string [num_regions];
  rtmData** meta_data_type = new rtmData* [num_regions];

  //IF GP DATA?
  if(meta_region->GP_data != 0)
    {
      for(int_i = 0; int_i < num_regions; int_i++)
	meta_data_type[int_i] = meta_region[int_i].GP_data;
      read_metaregion_datatype(meta_data_type, tempmat, countfiles, denomfiles, num_regions,
			       "regions_gp_count_data",
			       "regions_gp_coverage_data",
			       "regions_gp_aggregation", str_var, cTRUE);
      for(int_i = 0; int_i < num_regions; int_i++) // CHECK the actions of the following:
	gsl_matrix_realloc(meta_region[int_i].det_model_params.l_gp_negbin_overdispersion,
			   meta_data_type[int_i]->getDim1(),
			   meta_data_type[int_i]->getDim2());
    }
  // IF HOSP DATA?
  if(meta_region->Hospitalisation_data != 0)
    {
      for(int_i = 0; int_i < num_regions; int_i++)
	meta_data_type[int_i] = meta_region[int_i].Hospitalisation_data;
      read_metaregion_datatype(meta_data_type, tempmat, countfiles, denomfiles, num_regions,
			       "regions_hosp_data",
			       "",
			       "regions_hosp_aggregation", str_var, cFALSE);
    }
  // IF DEATHS DATA?
  if(meta_region->Death_data != 0)
    {
      for(int_i = 0; int_i < num_regions; int_i++)
	meta_data_type[int_i] = meta_region[int_i].Death_data;
      read_metaregion_datatype(meta_data_type, tempmat, countfiles, denomfiles, num_regions,
			       "regions_death_data",
			       "",
			       "regions_death_aggregation", str_var, cFALSE);
    }  // IF SERO DATA?
  if(meta_region->Serology_data != 0)
    {
      for(int_i = 0; int_i < num_regions; int_i++)
	meta_data_type[int_i] = meta_region[int_i].Serology_data;
      read_metaregion_datatype(meta_data_type, tempmat, countfiles, denomfiles, num_regions,
			       "regions_seropositives_data",
			       "regions_serosamples_data",
			       "regions_sero_aggregation", str_var, cFALSE);
    }
  // IF VIRO DATA?
  if(meta_region->Virology_data != 0)
    {
      for(int_i = 0; int_i < num_regions; int_i++)
	meta_data_type[int_i] = meta_region[int_i].Virology_data;
      read_metaregion_datatype(meta_data_type, tempmat, countfiles, denomfiles, num_regions,
			       "regions_viropositives_data",
			       "regions_virosamples_data",
			       "regions_viro_aggregation", str_var, cFALSE);
    }
      
  // Free all allocated memory
  delete [] countfiles;
  delete [] denomfiles;
  delete [] meta_data_type;
  gsl_matrix_free(tempmat);
}

void read_mixmod_structure_inputs(mixing_model& base_mm, const string str_input_filename, const global_model_instance_parameters src_gmp)
{

  // READ THE FILE INTO A FILE CONTENTS STRING
  string str_file_contents, str_var;

  fn_load_file(&str_file_contents, str_input_filename.c_str());

  // "CUT-OUT" THE RELEVANT VARIABLE
  int num_strings = cut_out_parameter(str_var, str_file_contents, "mixing_model");

  if(num_strings == 0)
    {
      ERROR_FILE_EXIT("Mixing structure not specified in file %s", str_input_filename.c_str());
    }

  // FIRSTLY NEED TO KNOW HOW MANY TEMPORAL BREAKPOINTS THERE ARE - THIS SHOULD ALLOCATE THE MEMORY TO THE BREAKPOINT COMPONENT.
  alloc_and_set_breakpoint_vector(base_mm.breakpoints,
				  base_mm.num_breakpoints,
				  str_var,
				  "time_breakpoints",
				  str_var.find("time_breakpoints"),
				  src_gmp.l_duration_of_runs_in_days,
				  true);

  base_mm.num_breakpoints--;
  // MEMORY ALLOCATION
  base_mm.MIXMAT = new gsl_matrix*[base_mm.num_breakpoints + 1];
  base_mm.MIXMAT_param = new gsl_matrix_int*[base_mm.num_breakpoints + 1];
  base_mm.MIXMAT_scaled = new gsl_matrix*[base_mm.num_breakpoints + 1];
  for(int int_i = 0; int_i <= base_mm.num_breakpoints; )
    {
      base_mm.MIXMAT[int_i] = gsl_matrix_alloc(NUM_AGE_GROUPS, NUM_AGE_GROUPS);
      base_mm.MIXMAT_scaled[int_i] = gsl_matrix_alloc(NUM_AGE_GROUPS, NUM_AGE_GROUPS);
      base_mm.MIXMAT_param[int_i++] = gsl_matrix_int_alloc(NUM_AGE_GROUPS, NUM_AGE_GROUPS);
    }
  base_mm.scalants = base_mm.eval_dominant = 0;
  base_mm.evector_MIXMAT_normalised = 0;

  // read in the sub-variables base_matrices and multiplier_indices_matrix
  read_filenames_and_data_matrices(base_mm.MIXMAT, base_mm.num_breakpoints + 1,
				   "base_matrices", str_var);
  read_filenames_and_data_int_matrices(base_mm.MIXMAT_param, base_mm.num_breakpoints + 1,
				       "multiplier_indices_matrices", str_var);

}
