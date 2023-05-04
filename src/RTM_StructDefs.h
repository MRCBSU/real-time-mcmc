#ifndef HEADER_StructDefs_
#define HEADER_StructDefs_

#include <vector>

#include "R_like_fns.h"
#include "RTM_flagclass.h"
#include "rtm_data_class.h"
#include "RTM_MvNorm.h"
#include "RTM_modelstate.h"
#include "gslWrapper.h"

using namespace std;
using std::string;

typedef enum { cREEDFROST, cMASSACTION } transmission_kernel;
typedef enum { cMCMC, cSMC } output_type;

// MODEL INSTANCE PARAMETERS
struct global_model_instance_parameters{
  data_type l_gp_count_likelihood; // DO WE HAVE POISSON OR NEGATIVE BINOMIAL DATA?
  data_type l_hosp_count_likelihood; // AS ABOVE
  transmission_kernel l_tk; // REED-FROST OR MASS-ACTION. REED-FROST BY DEFAULT.
  int l_transmission_time_steps_per_day; // DELTA T FOR THE TRANSMISSION MODEL
  int l_reporting_time_steps_per_day; // DELTA_T FOR THE REPORTING MODEL. SHOULD BE LESS THAN l_transmission_time_steps_per_day
  int l_duration_of_runs_in_days; // T_MAX (EPIDEMIC TIME)
  int l_projection_days; // NUMBER OF DAYS OVER WHICH TO PROJECT THE EPIDEMIC
  int l_day_of_start; // NUMBERED DAY OF YEAR ON WHICH T = 1
  int l_start_day_of_week;
  int l_start_week;
  gsl_vector_int* d_week_numbers_by_day; // THE WEEK NUMBER FOR EACH DAY OF THE STUDY. CALCULATED BASED ON THE ABOVE TWO MEMBERS OF THE STRUCTURE
  int l_num_regions;
  int l_GP_consultation_flag; // FLAGS INDICATING PRESENCE/ABSENCE OF DATA TYPES, THIS ROW DOWN
  int l_Hospitalisation_flag;
  int l_Deaths_flag;
  int l_Sero_data_flag;
  int l_Viro_data_flag;
  int l_Prev_data_flag;
  int l_Vacc_data_flag;
  int l_VBoost_data_flag;
  int l_VFourth_data_flag;
  int l_GP_patch_flag;
  int l_Sero_delay; // NUMBER OF DAYS OF LAG BUILT INTO SEROLOGY DATA REPRESENTING THE TIME TAKEN FOR IMMUNOLOGICAL RESPONSE
  likelihood_bounds l_GP_likelihood; // WHICH DAYS/WEEKS OF DATA TO USE IN CALCULATING THE LIKELIHOOD, THIS ROW DOWN
  likelihood_bounds l_Hosp_likelihood;
  likelihood_bounds l_Deaths_likelihood;
  likelihood_bounds l_Sero_likelihood;
  likelihood_bounds l_Viro_likelihood;
  likelihood_bounds l_Prev_likelihood;
  likelihood_bounds l_Vacc_date_range;
  likelihood_bounds l_VBoost_date_range;
  likelihood_bounds l_VFourth_date_range;
};

// ALLOC FUNCTION FOR MEMORY ASSIGNED TO THE ABOVE TYPE OF STRUCTURE
void alloc_global_model_instance(global_model_instance_parameters&);
// FREE FUNCTION FOR MEMORY ASSIGNED TO THE ABOVE TYPE OF STRUCTURE
void free_global_model_instance(global_model_instance_parameters&);

class regression_def{
public:
  gslVectorInt region_breakpoints;
  gslVectorInt age_breakpoints;
  gslVectorInt time_breakpoints;
  link_function regression_link; // enum
  gslMatrix design_matrix;
};

// GENERAL MODEL (POSSIBLY VECTOR VALUED) PARAMETER
struct updateable_model_parameter{
  string param_name;
  gsl_vector *param_value;
  gsl_vector *proposal_value;
  gsl_vector_int *prior_distribution; // THE DISTRIBUTION OF THE INDIVIDUAL COMPONENTS
  bool flag_hyperprior;
  double log_prior_dens;
  double proposal_log_prior_dens;
  bool flag_update; // TRUE: IF PRIOR_DISTRIBUTION[i] != cDETERM FOR ALL i
  gsl_vector **prior_params; // WHAT THESE PARAMETERS REPRESENT DEPENDS ON THE prior_distribution MEMEMBER OF THE STRUCTURE. EACH ROW CORRESPONDS TO EACH OF THE DISTRIBUTIONS IN THE prior_distribution VECTOR
  gsl_vector *proposal_variances;
  mvnorm *prior_multivariate_norm;
  regression_def map_to_regional;
  gsl_vector_int* number_accepted_moves;
  gsl_vector_int* number_proposed_moves;
  gsl_vector *posterior_mean;
  gsl_vector *posterior_sumsq;
  bool flag_transmission_model; // DOES THE NUMBER OF NEW INFECTEDS NEED TO RECALCULATED WHEN UPDATING THE PARAMETER THIS PARAMETER
  bool flag_reporting_model; // DO THE REPORTING DELAYS NEED TO BE RECALCULATED WHEN UPDATING WHEN UPDATING THIS PARAMETER
  bool flag_GP_likelihood; // DOES THE LIKELIHOOD FOR G.P. CONSULTATION NEED TO BE RECALCULATED WHEN UPDATING THIS PARAMETER
  bool flag_Hosp_likelihood; // DOES THE LIKELIHOOD FOR THE HOSPITALISATIONS NEED TO BE RECALCULATED WHEN UPDATING THIS PARAMETER
  bool flag_Sero_likelihood; // DOES THE LIKELIHOOD FOR THE SEROEPIDEMIOLOGY DATA NEED TO BE RECALCULATED WHEN UPDATING THIS PARAMETER
  bool flag_Viro_likelihood; // DOES THE LIKELIHOOD FOR THE VIROLOGY (POSITIVITY) DATA NEED TO BE RECALCULATED WHEN UPDATING THIS PARAMETER
  bool flag_Prev_likelihood; // DOES THE LIKELIHOOD FOR THE PREVALENCE DATA/ESTIMATES NEED TO BE RECALCULATED WHEN UPDATING THIS PARAMETER
  bool flag_any_child_nodes; // TRUE IF ANY FLAG_CHILD_NODES ARE TRUE. FALSE OTHERWISE
  bool *flag_child_nodes; // FLAG FOR EACH OF THE OTHER PARAMETERS OF THE MODEL INDICATING WHETHER THEY ARE CHILD NODES OF THE CURRENT NODE
};

// STRUCTURE DEFINED FOR MIXING MODEL
// TODO: Use this as regional definition so limited changes needed
struct mixing_model{
  int num_breakpoints;
  gsl_vector_int *breakpoints; // TIMING OF BREAKPOINTS, ASSUMED ORDERED AND INCREASING.
  gsl_matrix **MIXMAT;
  gsl_matrix **MIXMAT_scaled;
  gsl_matrix_int **MIXMAT_param;
  gsl_vector *scalants;
  gsl_vector *eval_dominant;
  gsl_vector **evector_MIXMAT_normalised;
};

void mixing_model_alloc(mixing_model&, const int, const int, const int num_strata = NUM_AGE_GROUPS);
void mixing_model_free(mixing_model&);
void mixing_model_memcpy(mixing_model&, const mixing_model&);

// PARAMETERS DEFINING THE GAMMA DISTRIBUTION GOVERNING DELAY TIMES IN THE DISEASE MODEL
class st_delay{
public:
  string delay_name;
  double gamma_mean;
  double gamma_sd;
  double gamma_shape;
  double gamma_rate;
  double gamma_scale;
};
//void st_delay_memcpy(st_delay&, const st_delay);

class infection_to_data_delay{
public:
  size_t num_components;
  std::vector<st_delay> component_delays;
  st_delay overall_delay;
  gslVector distribution_function; // CURRENTLY DISTRIBUTIONS ARE FIXED OVER AGES AND REGIONS, SO ONLY NEED THE TIME DIMENSION HERE

  // I don't think we can implement this as a constructor, as num_delays
  // needs to be read from input file
  void setSize(size_t num_delays, size_t num_times) {
    num_components = num_delays;
    component_delays.resize(num_delays);
    distribution_function.alloc(num_times);
  }
};

//void infection_to_data_delay_alloc(infection_to_data_delay&, size_t, size_t);
//void infection_to_data_delay_free(infection_to_data_delay&);


// GLOBAL PARAMETER STRUCTURE, FOR INPUT TO THE TRANSMISSION (AND DISEASE?) MODEL
class globalModelParams{
public:
  size_t size_param_list;
  updateable_model_parameter* param_list;
  infection_to_data_delay gp_delay;
  infection_to_data_delay hosp_delay;
  infection_to_data_delay death_delay;
};

void globalModelParams_alloc(globalModelParams&, size_t);
void globalModelParams_free(globalModelParams&);

// REGION SPECIFIC PARAMETER STRUCTURE, FOR INPUT TO THE TRANSMISSION (AND DISEASE?) MODEL
#define REGIONAL_MODEL_PARAMS_MEMBERS "l_init_prop_sus, l_init_prop_sus_HI_geq_32, l_average_infectious_period, l_latent_period, l_r1_period, l_vacc1_disease, l_vaccn_disease, l_vaccb_disease, l_vacc4_disease, l_vacc1_infect, l_vaccn_infect, l_vaccb_infect, l_vacc4_infect, l_relative_infectiousness_I2_wrt_I1, l_EGR, l_lbeta_rw, l_R0_Amplitude, l_R0_peakday, l_R0_init, l_I0, l_pr_symp, l_pr_onset_to_GP, l_pr_onset_to_Hosp, l_pr_onset_to_Death, l_importation_rate, d_R0_phase_differences, l_MIXMOD, l_background_gps_counts, l_sensitivity, l_specificity, l_gp_negbin_overdispersion, l_hosp_negbin_overdispersion, l_day_of_week_effect, l_sero_sensitivity, l_sero_specificity, l_waning_period;" 
struct regional_model_params{
  gsl_vector* l_init_prop_sus; // INITIAL CONDITION, MAKES NO SENSE TO HAVE ANY TEMPORAL VARIATION
  gsl_vector* l_init_prop_sus_HI_geq_32; // INITIAL CONDITION, MAKES NO SENSE TO HAVE ANY TEMPORAL VARIATION
  gsl_matrix* l_average_infectious_period; // TEMPORALLY VARYING TO ALLOW FOR CHANGES DUE TO ACCESS TO ANTIVIRALS
  gsl_matrix* l_latent_period;
  gsl_matrix* l_r1_period;
  gsl_matrix* l_vacc1_disease;
  gsl_matrix* l_vaccn_disease;
  gsl_matrix* l_vaccb_disease;
  gsl_matrix* l_vacc4_disease;
  gsl_matrix* l_vacc1_infect;
  gsl_matrix* l_vaccn_infect;
  gsl_matrix* l_vaccb_infect;
  gsl_matrix* l_vacc4_infect;
  gsl_matrix* l_relative_infectiousness_I2_wrt_I1;
  double l_EGR; // CURRENTLY NON-AGE DEPENDENT INITIAL EXPONENTIAL GROWTH RATE
  gsl_vector* l_lbeta_rw; // RANDOM-WALK SCALING TO APPLY TO THE WHOLE MATRIX, RANDOM WALKS OVER TIME.
  double l_R0_Amplitude;
  double l_R0_peakday; // DAY OF YEAR UPON WHICH VIRUS IS MOST VIRULENT
  double l_R0_init; // INITIAL CONDITION, MAKES NO SENSE TO HAVE ANY TEMPORAL VARIATION. AGE_VARIATION IS ON HOLD
  double l_I0; // INITIAL CONDITION, MAKES NO SENSE TO HAVE ANY TEMPORAL VARIATION. AGE VARIATION IS ON HOLD
  gsl_matrix* l_pr_symp;
  gsl_matrix* l_pr_onset_to_GP;
  gsl_matrix* l_pr_onset_to_Hosp;
  gsl_matrix* l_pr_onset_to_Death;
  gsl_matrix* l_importation_rate;
  gsl_vector* d_R0_phase_differences;
  mixing_model l_MIXMOD;
  gsl_matrix* l_background_gps_counts;
  double l_sensitivity;
  double l_specificity;
  gsl_matrix* l_gp_negbin_overdispersion;
  gsl_matrix* l_hosp_negbin_overdispersion;
  gsl_matrix* l_day_of_week_effect;
  gsl_matrix* l_sero_sensitivity;
  gsl_matrix* l_sero_specificity;
  gsl_matrix* l_waning_period;
};

void regional_model_params_alloc(regional_model_params&, const unsigned int, const int, const int, const int, const mixing_model);
void regional_model_params_memcpy(regional_model_params&, const regional_model_params&, flagclass&);
void regional_model_params_free(regional_model_params&);

// STRUCTURE CONTAINING OUTPUTS
struct model_statistics{
  model_state *d_end_state; // should contain a value for S(T, a), E_1(T, a), E_2(T, a), I_1(T, a), I_2(T, a) for all a = 1, ..., NUM_AGES
  gsl_matrix *d_NNI; // SHOULD INCORPORATE PROJECTION DAYS
  gsl_matrix *d_Delta_Dis;
  gsl_matrix *d_H1N1_GP_Consultations;
  gsl_matrix *d_Reported_GP_Consultations;
  gsl_matrix *d_Reported_Hospitalisations;
  gsl_matrix *d_internal_AR;
  gsl_matrix *d_seropositivity;
  gsl_matrix *d_viropositivity;
  gsl_matrix *d_prevalence;
};

void model_statistics_alloc(model_statistics&, const int, const int);
void model_statistics_aggregate(gsl_matrix*, const model_statistics&, const int);
void model_statistics_aggregate(gsl_matrix*, gsl_matrix*, const model_statistics&, const int);
void model_statistics_memcpy(model_statistics&, const model_statistics,
			     bool NNI_flag = true, bool GP_flag = true, bool Hosp_flag = true, bool Sero_flag = true, bool Viro_flag = true, bool Prev_flag = true);
void model_statistics_free(model_statistics&);

// REGION STRUCTURE
struct Region{
  string name;
  gsl_vector* population;
  double total_population;
  regional_model_params det_model_params;
  bool data_owner;
  rtmData* GP_data;
  rtmData* Hospitalisation_data;
  rtmData* Death_data;
  rtmData* Serology_data;
  rtmData* Virology_data;
  rtmData* Prevalence_data;
  rtmData* Vaccination_data;
  rtmData* VBoosting_data;
  rtmData* VFourth_data;
  model_statistics region_modstats;

  Region()
  : GP_data(0), Hospitalisation_data(0), Death_data(0), Serology_data(0), Virology_data(0), Prevalence_data(0), Vaccination_data(0), VBoosting_data(0), VFourth_data(0)
    { }
};

void Region_memcpy(Region&, const Region&, flagclass&);
void Region_alloc(Region&, const global_model_instance_parameters, const mixing_model);
void Region_alloc(Region&, const Region&);
void Region_free(Region&, const global_model_instance_parameters);

// MCMC SIMULATION PARAMETERS
struct mcmcPars{
  int num_iterations;
  output_type oType;
  int thin_output_every;
  int thin_stats_every;
  int adaptive_phase;
  int adapt_every;
  int burn_in;
  int read_covar;
  int global_updates;
  int num_progress_reports;
  double mixing_threshold_ub;
  double mixing_threshold_lb;
  double prop_var_inflation;
  double prop_var_shrinkage;
  int maximum_block_size;
  int random_seed;
  int max_num_threads;
  gsl_rng* r;
};

void mcmcPars_alloc(mcmcPars&);
void mcmcPars_free(mcmcPars&);

// LIKELIHOOD STRUCTURE



class rlikelihood {
public:
  double region_lfx;
  double GP_lfx;
  double Hosp_lfx;
  double Deaths_lfx;
  double Sero_lfx;
  double Viro_lfx;
  double Prev_lfx;

  rlikelihood()
    : region_lfx(0), GP_lfx(0), Hosp_lfx(0), Deaths_lfx(0),
      Sero_lfx(0), Viro_lfx(0), Prev_lfx(0) { }
};

class glikelihood {
public:
  double total_lfx;
  double bar_lfx;
  double sumsq_lfx;

  std::vector<rlikelihood> rlik;

  glikelihood()
    : total_lfx(0), bar_lfx(0), sumsq_lfx(0) { }

  glikelihood(const global_model_instance_parameters &gmip)
    : total_lfx(0), bar_lfx(0), sumsq_lfx(0),
      rlik(gmip.l_num_regions) { }
};

/*
// CCS: Rewritten as a class to remove need for manual memory management
class likelihood {
public:
  double total_lfx;
  double bar_lfx;
  double sumsq_lfx;
  gslVector region_lfx;
  gslVector GP_lfx;
  gslVector Hosp_lfx;
  gslVector Deaths_lfx;
  gslVector Sero_lfx;
  gslVector Viro_lfx;
  gslVector Prev_lfx;

  // Empty constructor necessary for updParamBlock
  likelihood()
    : total_lfx(0), bar_lfx(0), sumsq_lfx(0) { }
  
  likelihood(const global_model_instance_parameters &gmip)
    : total_lfx(0), bar_lfx(0), sumsq_lfx(0) {
    int num_regions = gmip.l_num_regions;
    region_lfx.allocZero(num_regions);
    if (gmip.l_GP_consultation_flag)
      GP_lfx.allocZero(num_regions);
    if (gmip.l_Hospitalisation_flag)
      Hosp_lfx.allocZero(num_regions);
    if (gmip.l_Deaths_flag)
      Deaths_lfx.allocZero(num_regions);
    if (gmip.l_Sero_data_flag)
      Sero_lfx.allocZero(num_regions);
    if (gmip.l_Viro_data_flag)
      Viro_lfx.allocZero(num_regions);
    if (gmip.l_Prev_data_flag)
      Prev_lfx.allocZero(num_regions);
  }
};
*/
//void likelihood_alloc(likelihood&, const global_model_instance_parameters);
//void likelihood_free(likelihood&);
//void likelihood_memcpy(likelihood&, const likelihood&);

#endif
