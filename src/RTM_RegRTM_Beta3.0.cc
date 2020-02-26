#include "args.h"
#include "RTM_StructDefs.h" // provides RTM_Header.h
#include "RTM_StructAssign.h"
#include "RTM_FunctDefs.h"
#include "RTM_flagclass.h"

#include <iostream>

using namespace std;
using std::string;

/// Parse command line arguments for the input and output directories.
/// Returns a tuple of {input directory, output directory}, which will be strings.
/// Prints help string and exits if that is requested.
/// Also causes exit if error in arguments.
std::tuple<string,string> parse_command_line_arguments(int argc, char **argv) {
  args::ArgumentParser parser("Implmentation of Birrel et al. 2016");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  args::ValueFlag<string> input_dir(parser, "input",
    "The directory containing the input files mod_pars.txt and mod_input.txt",
    {'i', "input"}, "inputs"
  );
  args::ValueFlag<string> output_dir(parser, "output",
    "The directory that outputs will be placed in."
    "The directory be created if it doesn't exist.",
    {'o', "output"}, "outputs"
  );

  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help) {
    std::cout << parser;
    exit(0);
  } catch (args::Error e) {
    std::cerr << e.what() << std::endl << parser;
    exit(1);
  }
  return {input_dir.Get(), output_dir.Get()};
}

int main(int argc, char **argv) {

  global_model_instance_parameters global_fixedpars;
  globalModelParams global_modpars;
  mcmcPars sim_pars;
  likelihood llhood;  
  struct tm *now;
  time_t tval;

  // SET :"TIMER" RUNNING
  tval = time(NULL);
  now = localtime(&tval);

  auto [input_dir, output_dir] = parse_command_line_arguments(argc, argv);
  DEBUG(DEBUG_DETAIL, "Putting inputs in: " << input_dir);
  DEBUG(DEBUG_DETAIL, "Putting outputs in: " << output_dir);

  /// ALLOCATE STRUCTURES FEATURED IN read_model_inputs
  string str_global_model_instance_members(GLOBAL_MODEL_INSTANCE_MEMBERS);
  string str_global_model_instance_defaults(GLOBAL_MODEL_INSTANCE_DEFAULTS);

  string str_filename_inputs = input_dir + "/mod_inputs.txt";
  try{ 
    // BELOW FUNCTION WILL ALLOC MEMORY TO GLOBAL_FIXEDPARS
    read_global_fixed_parameters(global_fixedpars,
                str_filename_inputs.c_str(),
                str_global_model_instance_members,
                str_global_model_instance_defaults);
  } catch (std::fstream::failure e) {
    DEBUG(DEBUG_ERROR, "Cannot read " << str_filename_inputs << " :" << e.what());
    DEBUG(DEBUG_ERROR, strerror(errno));
    exit(2);
  }

  // // //


  /// ALLOCATE AND INITIALISE THE UPDATEABLE MODEL PARAMETERS STRUCTURE
  string str_global_model_parameters_members(GLOBAL_MODEL_PARAMETERS_MEMBERS);
  string str_global_model_parameters_initvals(GLOBAL_MODEL_PARAMETERS_DEFAULT_FIXED_VALS);
  string str_global_model_parameters_likflags(GLOBAL_MODEL_PARAMETERS_LIKELIHOOD_FLAGS);
  string str_global_model_delay_members(GLOBAL_MODEL_PARAMETERS_DELAY_NAMES);
  string str_global_model_delay_means(GLOBAL_MODEL_PARAMETERS_DELAY_MEANS);
  string str_global_model_delay_sds(GLOBAL_MODEL_PARAMETERS_DELAY_SDS);
  string str_global_model_delay_flags(GLOBAL_MODEL_PARAMETERS_DELAY_FLAGS);

  string str_filename_modpars = input_dir + "/mod_pars.txt";
  try{
    // BELOW FUNCTION WILL ALLOC MEMORY TO GLOBAL_MODPARS, AND SET TO FILE SPECIFIED VALUES OR DEFAULTS
    read_global_model_parameters(global_modpars,
                str_filename_modpars.c_str(),
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
  } catch (std::fstream::failure e) {
    DEBUG(DEBUG_ERROR, "Cannot read " << str_filename_modpars << " :" << e.what());
    exit(2);
  }
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
  		       str_filename_inputs.c_str(),
  		       str_mcmc_parameter_names,
  		       str_mcmc_parameter_defaults);

  // SET THE MAXIMUM NUMBER OF PARALLEL THREADS
  #ifdef USE_THREADS
  omp_set_num_threads(sim_pars.max_num_threads);
  #endif

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
  	      sim_pars.r,
          output_dir);

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
