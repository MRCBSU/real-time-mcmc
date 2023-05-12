#ifndef HEADER_FunctDefs_
#define HEADER_FunctDefs_

#include "RTM_Header.h"
#include "RTM_flagclass.h"
#include "rtm_data_class.h"
#include "RTM_StructDefs.h"
#include "RTM_updParams.h"

using namespace std;
using std::string;

// FUNCTIONS IN RTM_Inputs.cc
void input_filenames(char *,
		     char *, const int);
void read_global_fixed_parameters(global_model_instance_parameters&,
				  char*, const string, string&);
void read_global_model_parameters(globalModelParams&, updParamSet &, const char*, const string,
				  const string&, const string&,	const string&, const string&,
				  const string&, const string&, const int, const int, const int,
				  const double);
int param_list_index(const globalModelParams, string);
void read_string_from_instruct(string&, const string&, const string&);
void read_data_inputs(Region*, const string, const int&);
//TODO: Replace function so that it reads in correct regional results
void read_mixmod_structure_inputs(mixing_model&, const string, const global_model_instance_parameters);
void read_mixmod_structure_inputs_list(std::vector<std::unique_ptr<mixing_model>>&, const string, const global_model_instance_parameters);
void read_mcmc_parameters(mcmcPars&, char*, const string, string&);
void data_matrices_fscanf(const string&, gsl_matrix*, const int col_skip = 0, const int row_skip = 0);

// FUNCTIONS IN RTM_WithinRegion.cc
void select_design_matrix(gslMatrix&, const gslMatrix&, bool, int, int);
void mat_breakpoint_cut(gsl_matrix*, const gsl_vector_int*, const gsl_vector_int*, const gsl_vector*);

void block_regional_parameters(regional_model_params&, const updParamSet &, const global_model_instance_parameters&, const int&, const gsl_vector*, const double&, const mixing_model&, flagclass&);
void block_regional_parameters(regional_model_params&, const updParamSet &, const global_model_instance_parameters&, const int&, const gsl_vector*, const double&, const std::unique_ptr<mixing_model>&, flagclass&);
void evaluate_regional_parameters(regional_model_params&, const updateable_model_parameter*, const global_model_instance_parameters&,
				  const int&, const gsl_vector*, const double&, const mixing_model&,
				  flagclass&);

// FUNCTIONS IN RTM_MixingMatrices.cc
int maximal_index(const mixing_model&);
int maximal_index(const std::unique_ptr<mixing_model>&);
void mixing_matrix_scale_and_normalisation(mixing_model&, const gsl_vector*, const bool evec_flag = true);
int mix_timecat(const int, const gsl_vector_int*, const double);
void mixing_matrix_parameterise(mixing_model&);

// FUNCTIONS IN RTM_Likelihoods.cc
//void block_log_likelihood(likelihood& llhood, Region* meta_region, int region, bool flag_update_transmission_model, bool flag_update_reporting_model, bool flag_update_GP_likelihood, bool flag_update_Hosp_likelihood, bool flag_update_Viro_likelihood, bool flag_update_Sero_likelihood, bool flag_update_Prev_likelihood, const global_model_instance_parameters &gmip, updParamSet &pars, bool inBlock = false);

void fn_log_likelihood_global(glikelihood&, Region*, int, bool, bool, bool, bool, bool, bool, bool,
		       const global_model_instance_parameters &, gslVector&, gslVector&);
void fn_log_likelihood_region(rlikelihood&, Region*, int, bool, bool, bool, bool, bool, bool, bool,
		       const global_model_instance_parameters &, gslVector&, gslVector&);
void output_per_selected_period(const int&, const gsl_matrix*, gsl_matrix*);
double fn_log_lik_positivity(const gsl_matrix*, const gsl_matrix*, const gsl_matrix*);
double fn_log_lik_countdata(const gsl_matrix*, const gsl_matrix*);
double fn_log_lik_negbindata(const gsl_matrix*, const gsl_matrix*, const gsl_matrix*);
double fn_log_lik_loggaussian_fixedsd(const gsl_matrix*, const gsl_matrix*, const gsl_matrix*);

// FUNCTIONS IN RTM_MetropHast.cc
void metrop_hast(const mcmcPars&, globalModelParams&, updParamSet &, Region*, Region*, glikelihood&, const global_model_instance_parameters&, const std::vector<std::unique_ptr<mixing_model>>&, gsl_rng*);

#endif
