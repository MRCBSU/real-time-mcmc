#include "RTM_StructDefs.h" // provides RTM_Header.h
#include "RTM_StructAssign.h"
#include "RTM_FunctDefs.h"
#include "RTM_flagclass.h"

using namespace std;
using std::string;

int main(void){

  global_model_instance_parameters global_fixedpars;
  globalModelParams global_modpars;
  mcmcPars sim_pars;
  likelihood llhood;  
  struct tm *now;
  time_t tval;

  // SET :"TIMER" RUNNING
  tval = time(NULL);
  now = localtime(&tval);

  // GET FILES CONTAINING USER INPUT, IF ANY EXIST OR USE DEFAULT FILENAMES.
  // SPECIFIED FILENAMES SHOULD BE FOUND IN THE FILE rtm_input_files.txt
  // DEFAULTS FOUND IN ARGUMENT DECLARATION FOR input_filenames IN RTM_Inputs.cc
  char str_filename_inputs[120] = "mod_inputs.txt";
  char str_filename_modpars[120] = "mod_pars.txt";

  // INSTRUCT THE PROGRAM AS TO WHERE TO LOOK FOR VARIOUS GROUPS OF MODEL INPUTS
  input_filenames(str_filename_inputs,
  		  str_filename_modpars,
  		  120);
  // //

  // // READ IN THE VARIOUS INPUTS.

  /// ALLOCATE STRUCTURES FEATURED IN read_model_inputs
  string str_global_model_instance_members(GLOBAL_MODEL_INSTANCE_MEMBERS);
  string str_global_model_instance_defaults(GLOBAL_MODEL_INSTANCE_DEFAULTS);

  // BELOW FUNCTION WILL ALLOC MEMORY TO GLOBAL_FIXEDPARS
  read_global_fixed_parameters(global_fixedpars,
  			       str_filename_inputs,
  			       str_global_model_instance_members,
  			       str_global_model_instance_defaults);

  // // //


  /// ALLOCATE AND INITIALISE THE UPDATEABLE MODEL PARAMETERS STRUCTURE
  string str_global_model_parameters_members(GLOBAL_MODEL_PARAMETERS_MEMBERS);
  string str_global_model_parameters_initvals(GLOBAL_MODEL_PARAMETERS_DEFAULT_FIXED_VALS);
  string str_global_model_parameters_likflags(GLOBAL_MODEL_PARAMETERS_LIKELIHOOD_FLAGS);
  string str_global_model_delay_members(GLOBAL_MODEL_PARAMETERS_DELAY_NAMES);
  string str_global_model_delay_means(GLOBAL_MODEL_PARAMETERS_DELAY_MEANS);
  string str_global_model_delay_sds(GLOBAL_MODEL_PARAMETERS_DELAY_SDS);
  string str_global_model_delay_flags(GLOBAL_MODEL_PARAMETERS_DELAY_FLAGS);

  // BELOW FUNCTION WILL ALLOC MEMORY TO GLOBAL_MODPARS, AND SET TO FILE SPECIFIED VALUES OR DEFAULTS
  read_global_model_parameters(global_modpars,
  			       str_filename_modpars,
  			       str_global_model_parameters_members,
  			       str_global_model_parameters_initvals,
  			       str_global_model_parameters_likflags,
  			       str_global_model_delay_members,
  			       str_global_model_delay_means,
  			       str_global_model_delay_sds,
  			       str_global_model_delay_flags,
  			       global_fixedpars.l_num_regions,
  			       global_fixedpars.l_duration_of_runs_in_days,
  			       NUM_AGE_GROUPS,
  			       global_fixedpars.l_reporting_time_steps_per_day);
  // //

  // GOING TO READ IN THE DATA FOR EACH REGION. SET UP A META-REGION
  Region* country = new Region[global_fixedpars.l_num_regions];

  // BELOW ROUTINES WILL READ IN ALL THE DATA THAT WE'RE GOING TO USE
  // POLYMOD MATRICES ALSO TREATED AS DATA. THE POLYMOD "DATA" FILE
  // SHOULD BE SPECIFIED IN THE FILE NAMED str_filename_modpars

  // FIRST, WANT TO READ IN THE MIXING MODEL
  mixing_model mixmod_struct;
  read_mixmod_structure_inputs(mixmod_struct, str_filename_inputs, global_fixedpars);

  // ALLOCATE MEMORY TO REGIONAL SUBSTRUCTURES
  for(int int_i = 0; int_i < global_fixedpars.l_num_regions; int_i++)
    Region_alloc(country[int_i], global_fixedpars, mixmod_struct);

  // ESTABLISH THE POINTER TO ALL THE REGIONS WHICH WILL STORE ALL THE DATA.
  read_data_inputs(country, str_filename_inputs, global_fixedpars.l_num_regions);

  // NOW NEED TO INITIALISE THE MODEL TO MAKE AN INITIAL CALCULATION OF THE STARTING POINT LIKELIHOOD

  // INITIALISE EACH OF THE REGIONAL MODELS
  flagclass all_true;
  for(int int_i = 0; int_i < global_fixedpars.l_num_regions; int_i++)
    evaluate_regional_parameters(country[int_i].det_model_params, global_modpars.param_list, 
  				 global_fixedpars, int_i, country[int_i].population,
  				 country[int_i].total_population, mixmod_struct, all_true);


  // READ IN THE PARAMETERS OF THE MCMC SIMULATION
  string str_mcmc_parameter_names(MCMC_PARAMETER_NAMES);
  string str_mcmc_parameter_defaults(MCMC_PARAMETER_DEFAULT_VALUES);

  read_mcmc_parameters(sim_pars,
  		       str_filename_inputs,
  		       str_mcmc_parameter_names,
  		       str_mcmc_parameter_defaults);

  // SET THE MAXIMUM NUMBER OF PARALLEL THREADS
  // omp_set_num_threads(sim_pars.max_num_threads);

  // INITIALISE THE LIKELIHOOD STRUCTURE
  likelihood_alloc(llhood, global_fixedpars);

  // MAKE AN INITIAL EVALUATION OF THE LIKELIHOOD
  fn_log_likelihood(llhood, country, 0, true, true,
  		    global_fixedpars.l_GP_consultation_flag,
  		    global_fixedpars.l_Hospitalisation_flag,
  		    global_fixedpars.l_Viro_data_flag,
  		    global_fixedpars,
  		    global_modpars);

  // RUN THE METROPOLIS-HASTINGS SAMPLER AND SEND OUTPUT TO FILES
  metrop_hast(sim_pars,
  	      global_modpars,
  	      country,
  	      llhood,
  	      global_fixedpars,
  	      mixmod_struct,
  	      sim_pars.r);

  ////////////////////////////////////////////////////////////////////////////////////////

  // FREE THE STRUCTURE CONTAINING THE LIKELIHOOD DETAILS
  likelihood_free(llhood);

  // FREE THE STRUCTURE CONTAINING THE MCMC SIMULATION PARAMETERS
  mcmcPars_free(sim_pars);

  // FREE THE MIXING MODEL STRUCTURE INTO WHICH THE MATRICES WERE INITIALLY READ
  mixing_model_free(mixmod_struct);

  // FREE THE SUBSTRUCTURES WITHIN THE COUNTRY ARRAY
  for(int int_i = 0; int_i < global_fixedpars.l_num_regions; int_i++)
    Region_free(country[int_i], global_fixedpars);

  /// FREE THE COUNTRY ARRAY
  delete [] country;

  /// FREE STRUCTURES ALLOCATED AHEAD OF read_model_inputs CALL
  free_global_model_instance(global_fixedpars);

  /// FREE GLOBAL PARAMETERS STRUCTURE
  globalModelParams_free(global_modpars);

  // ////////////////////////////////////////////////////////////////////////////////////////

  // GET ELAPSED TIME
  printf("Elapsed program time in hours: %.4f\n", difftime(time(NULL), tval) / (60 * 60));

}
