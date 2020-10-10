#ifndef HEADER_FunctDefs_
#define HEADER_FunctDefs_

#include "RTM_Header.h"
#include "RTM_flagclass.h"
#include "rtm_data_class.h"
#include "RTM_StructDefs.h"

using namespace std;
using std::string;

// FUNCTIONS IN RTM_Inputs.cc
void input_filenames(char *,
		     char *, const int);
void read_global_fixed_parameters(register global_model_instance_parameters&,
				  char*, const string, string&);
void read_global_model_parameters(globalModelParams&, const char*, const string,
				  const string&, const string&,	const string&, const string&,
				  const string&, const string&, const int, const int, const int,
				  const double);
int param_list_index(const globalModelParams, string);
void read_string_from_instruct(string&, const string&, const string&);
void read_data_inputs(Region*, const string, const int&);
void read_mixmod_structure_inputs(mixing_model&, const string, const global_model_instance_parameters);
void read_mcmc_parameters(register mcmcPars&, char*, const string, string&);
void data_matrices_fscanf(const string&, gsl_matrix*, const int col_skip = 0, const int row_skip = 0);

// FUNCTIONS IN RTM_WithinRegion.cc
void select_design_matrix(gsl_matrix*, gsl_matrix*, bool, int, int);
void mat_breakpoint_cut(gsl_matrix*, const gsl_vector_int*, const gsl_vector_int*, const gsl_vector*);
void evaluate_regional_parameters(regional_model_params&, const updateable_model_parameter*, const global_model_instance_parameters&,
				  const int&, const gsl_vector*, const double&, const mixing_model&,
				  flagclass&);

// FUNCTIONS IN RTM_MixingMatrices.cc
int maximal_index(const mixing_model);
void mixing_matrix_scale_and_normalisation(mixing_model&, const gsl_vector*, const bool evec_flag = true);
int mix_timecat(const int, const gsl_vector_int*, const double);
void mixing_matrix_parameterise(mixing_model&);

// FUNCTIONS IN RTM_Likelihoods.cc
void fn_log_likelihood(likelihood&, Region*, int, bool, bool, bool, bool, bool, bool, bool,
		       global_model_instance_parameters, globalModelParams);
void output_per_selected_period(const int&, const gsl_matrix*, gsl_matrix*);
double fn_log_lik_positivity(const gsl_matrix*, const gsl_matrix*, const gsl_matrix*);
double fn_log_lik_countdata(const gsl_matrix*, const gsl_matrix*);
double fn_log_lik_negbindata(const gsl_matrix*, const gsl_matrix*, const gsl_matrix*);
double fn_log_lik_loggaussian_fixedsd(const gsl_matrix*, const gsl_matrix*, const gsl_matrix*);

// FUNCTIONS IN RTM_MetropHast.cc
void metrop_hast(const mcmcPars&, globalModelParams&, Region*, likelihood&, const global_model_instance_parameters&, const mixing_model&, gsl_rng*);

#endif
