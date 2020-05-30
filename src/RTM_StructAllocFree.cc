#include "RTM_StructDefs.h"
#include "RTM_FunctDefs.h"
#include "RTM_flagclass.h"

using namespace std;
using std::string;


void likelihood_alloc(likelihood& new_llhood, const global_model_instance_parameters gmip)
{
  int num_regions = gmip.l_num_regions;
  new_llhood.total_lfx = new_llhood.bar_lfx = new_llhood.sumsq_lfx = 0.0;
  new_llhood.GP_lfx = (gmip.l_GP_consultation_flag) ? gsl_vector_alloc(num_regions) : 0;
  if(new_llhood.GP_lfx != 0) gsl_vector_set_zero(new_llhood.GP_lfx); 
  new_llhood.Hosp_lfx = (gmip.l_Hospitalisation_flag) ? gsl_vector_alloc(num_regions) : 0;
  if(new_llhood.Hosp_lfx != 0) gsl_vector_set_zero(new_llhood.Hosp_lfx);
  new_llhood.Deaths_lfx = (gmip.l_Deaths_flag) ? gsl_vector_alloc(num_regions) : 0;
  if(new_llhood.Deaths_lfx != 0) gsl_vector_set_zero(new_llhood.Deaths_lfx);
  new_llhood.Sero_lfx = (gmip.l_Sero_data_flag) ? gsl_vector_alloc(num_regions) : 0;
  if(new_llhood.Sero_lfx != 0) gsl_vector_set_zero(new_llhood.Sero_lfx);
  new_llhood.Viro_lfx = (gmip.l_Viro_data_flag) ? gsl_vector_alloc(num_regions) : 0;
  if(new_llhood.Viro_lfx != 0) gsl_vector_set_zero(new_llhood.Viro_lfx);
}

void likelihood_free(likelihood& old_llhood)
{
  if(old_llhood.GP_lfx != 0)
    gsl_vector_free(old_llhood.GP_lfx);
  if(old_llhood.Hosp_lfx != 0)
    gsl_vector_free(old_llhood.Hosp_lfx);
  if(old_llhood.Deaths_lfx != 0)
    gsl_vector_free(old_llhood.Deaths_lfx);
  if(old_llhood.Sero_lfx != 0)
    gsl_vector_free(old_llhood.Sero_lfx);
  if(old_llhood.Viro_lfx != 0)
    gsl_vector_free(old_llhood.Viro_lfx);
}

void likelihood_memcpy(likelihood& dest_lfx, const likelihood& src_lfx)
{
  dest_lfx.total_lfx = src_lfx.total_lfx;
  dest_lfx.bar_lfx = src_lfx.bar_lfx;
  dest_lfx.sumsq_lfx = src_lfx.sumsq_lfx;
  if(src_lfx.GP_lfx != 0)
    gsl_vector_memcpy(dest_lfx.GP_lfx, src_lfx.GP_lfx);
  if(src_lfx.Hosp_lfx != 0)
    gsl_vector_memcpy(dest_lfx.Hosp_lfx, src_lfx.Hosp_lfx);
  if(src_lfx.Deaths_lfx != 0)
    gsl_vector_memcpy(dest_lfx.Deaths_lfx, src_lfx.Deaths_lfx);
  if(src_lfx.Sero_lfx != 0)
    gsl_vector_memcpy(dest_lfx.Sero_lfx, src_lfx.Sero_lfx);
  if(src_lfx.Viro_lfx != 0)
    gsl_vector_memcpy(dest_lfx.Viro_lfx, src_lfx.Viro_lfx);
}

void mcmcPars_alloc(mcmcPars& in_mP)
{
  register const gsl_rng_type* T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  in_mP.r = gsl_rng_alloc(T);
  gsl_rng_set(in_mP.r, in_mP.random_seed);
}

void mcmcPars_free(mcmcPars& out_mP)
{
  gsl_rng_free(out_mP.r);
}

void alloc_global_model_instance(global_model_instance_parameters& in_globalModelInstanceParameters)
{
  in_globalModelInstanceParameters.d_week_numbers_by_day = gsl_vector_int_alloc(in_globalModelInstanceParameters.l_duration_of_runs_in_days);
}

void free_global_model_instance(global_model_instance_parameters& in_globalModelInstanceParameters)
{
  gsl_vector_int_free(in_globalModelInstanceParameters.d_week_numbers_by_day);
}

void free_regression_def(regression_def& old_regdef)
{

  if(old_regdef.region_breakpoints != 0)
    gsl_vector_int_free(old_regdef.region_breakpoints);
  if(old_regdef.age_breakpoints != 0)
    gsl_vector_int_free(old_regdef.age_breakpoints);
  if(old_regdef.time_breakpoints != 0)
    gsl_vector_int_free(old_regdef.time_breakpoints);
  if(old_regdef.design_matrix != 0)
    gsl_matrix_free(old_regdef.design_matrix);
}

void updateable_model_parameter_free(updateable_model_parameter& old_ump)
{
  if(old_ump.flag_update)
    {
      gsl_vector_free(old_ump.proposal_value);
      gsl_vector_free(old_ump.proposal_variances);
      gsl_vector_free(old_ump.posterior_mean);
      gsl_vector_free(old_ump.posterior_sumsq);
      gsl_vector_int_free(old_ump.number_accepted_moves);
      gsl_vector_int_free(old_ump.number_proposed_moves);
      if(!old_ump.flag_hyperprior){
	for(int int_i = 0; int_i < old_ump.param_value->size; int_i++)
	  {
	    if(old_ump.prior_params[int_i] != 0)	
	      gsl_vector_free(old_ump.prior_params[int_i]);
	  }
      }
      free(old_ump.prior_params);
    }
  delete [] old_ump.flag_child_nodes;
  free_regression_def(old_ump.map_to_regional);
  gsl_vector_free(old_ump.param_value);
  gsl_vector_int_free(old_ump.prior_distribution);
}

 void infection_to_data_delay_alloc(infection_to_data_delay& new_itdd, size_t num_delays, size_t num_times)
{
  new_itdd.num_components = num_delays;
  new_itdd.component_delays = new st_delay[num_delays];
  if(num_times > 0)
    new_itdd.distribution_function = gsl_vector_alloc(num_times);
  else new_itdd.distribution_function = 0;
}

void infection_to_data_delay_free(infection_to_data_delay& old_itdd)
{
  delete [] old_itdd.component_delays;
  if(old_itdd.distribution_function != 0)
    gsl_vector_free(old_itdd.distribution_function);
}

void st_delay_memcpy(st_delay& new_st_d, const st_delay old_st_d)
{

  new_st_d.delay_name.assign(old_st_d.delay_name);
  new_st_d.gamma_mean = old_st_d.gamma_mean;
  new_st_d.gamma_sd = old_st_d.gamma_sd;
  new_st_d.gamma_shape = old_st_d.gamma_shape;
  new_st_d.gamma_rate = old_st_d.gamma_rate;

}


void globalModelParams_alloc(globalModelParams& new_gmp, size_t num_params)
{
  new_gmp.size_param_list = num_params;
  new_gmp.param_list = new updateable_model_parameter[num_params];
  for(int inti = 0; inti < num_params; inti++)
    {
      new_gmp.param_list[inti].proposal_value = 0;
      new_gmp.param_list[inti].flag_hyperprior = false;
      new_gmp.param_list[inti].log_prior_dens = 0.0;
      new_gmp.param_list[inti].proposal_log_prior_dens = 0.0;
      new_gmp.param_list[inti].flag_update = false;
      new_gmp.param_list[inti].prior_params = 0;
      new_gmp.param_list[inti].proposal_variances = 0;
      new_gmp.param_list[inti].number_accepted_moves = 0;
      new_gmp.param_list[inti].number_proposed_moves = 0;
      new_gmp.param_list[inti].posterior_mean = 0;
      new_gmp.param_list[inti].posterior_sumsq = 0;
      new_gmp.param_list[inti].flag_child_nodes = new bool[num_params];
      new_gmp.param_list[inti].flag_any_child_nodes = false;
      new_gmp.param_list[inti].prior_multivariate_norm = 0;
      for(int intj = 0; intj < num_params; intj++)
	new_gmp.param_list[inti].flag_child_nodes[intj] = false;
    }
}

void globalModelParams_free(globalModelParams& old_gmp)
{
  for(int int_i = 0; int_i < old_gmp.size_param_list; int_i++)
    {
      updateable_model_parameter_free(old_gmp.param_list[int_i]);
    }
  delete [] old_gmp.param_list;
  infection_to_data_delay_free(old_gmp.gp_delay);
  infection_to_data_delay_free(old_gmp.hosp_delay);
  infection_to_data_delay_free(old_gmp.death_delay);

}
// ---- overloaded regional_model_params_alloc function; allocation of memory for regional_model_params_alloc objects.
void regional_model_params_alloc(regional_model_params& new_rmp,
				 const int num_days,
				 const int num_ages,
				 const int transmission_time_steps_per_day,
				 const int reporting_time_steps_per_day,
				 const mixing_model src_mixing_pars)
{

  new_rmp.l_init_prop_sus = gsl_vector_calloc(num_ages);
  new_rmp.l_init_prop_sus_HI_geq_32 = gsl_vector_calloc(num_ages);
  new_rmp.l_average_infectious_period = gsl_matrix_calloc(transmission_time_steps_per_day * num_days, num_ages); /// DO I WANT NUM_DAYS AND NUM_AGES? OR THE NUMBER OF TEMPORAL AND AGE BREAKPOINTS.
  new_rmp.l_latent_period = gsl_matrix_calloc(transmission_time_steps_per_day * num_days, num_ages); /// THE SAME GOES FOR MANY OF THE OTHER MATRICES
  new_rmp.l_relative_infectiousness_I2_wrt_I1 = gsl_matrix_calloc(transmission_time_steps_per_day * num_days, num_ages); // VARIATION NOT EXPECTED TO BE USED HERE
  new_rmp.l_lbeta_rw = gsl_vector_calloc(transmission_time_steps_per_day * num_days);
  new_rmp.l_pr_symp = gsl_matrix_calloc(reporting_time_steps_per_day * num_days, num_ages);
  new_rmp.l_pr_onset_to_GP = gsl_matrix_calloc(reporting_time_steps_per_day * num_days, num_ages);
  new_rmp.l_pr_onset_to_Hosp = gsl_matrix_calloc(reporting_time_steps_per_day * num_days, num_ages);
  new_rmp.l_pr_onset_to_Death = gsl_matrix_calloc(reporting_time_steps_per_day * num_days, num_ages);
  new_rmp.l_importation_rate = gsl_matrix_calloc(transmission_time_steps_per_day * num_days, num_ages);
  new_rmp.d_R0_phase_differences = gsl_vector_calloc(transmission_time_steps_per_day * num_days);

  int num_mix_breakpoints = (src_mixing_pars.breakpoints == 0) ? 0 : src_mixing_pars.breakpoints->size;

  // GET THE SIZE OF THE PARAMETER VECTOR FROM THE MAXIMAL INDEX FOUND IN THE PARAMETERISATION MATRICES
  int param_vec_size = maximal_index(src_mixing_pars) + 1;

  mixing_model_alloc(new_rmp.l_MIXMOD,
		     num_mix_breakpoints,
		     param_vec_size);

  new_rmp.l_background_gps_counts = gsl_matrix_calloc(num_days, num_ages);
  new_rmp.l_gp_negbin_overdispersion = gsl_matrix_calloc(num_days, num_ages);
  new_rmp.l_hosp_negbin_overdispersion = gsl_matrix_calloc(num_days, num_ages);
  new_rmp.l_day_of_week_effect = gsl_matrix_calloc(num_days, num_ages);
}
void regional_model_params_alloc(regional_model_params& dest_rmp,
				 const regional_model_params& src_rmp)
{

  dest_rmp.l_init_prop_sus = gsl_vector_calloc(src_rmp.l_init_prop_sus->size);
  dest_rmp.l_init_prop_sus_HI_geq_32 = gsl_vector_calloc(src_rmp.l_init_prop_sus_HI_geq_32->size);
  dest_rmp.l_average_infectious_period = gsl_matrix_calloc(src_rmp.l_average_infectious_period->size1, src_rmp.l_average_infectious_period->size2);
  dest_rmp.l_latent_period = gsl_matrix_calloc(src_rmp.l_latent_period->size1, src_rmp.l_latent_period->size2);
  dest_rmp.l_relative_infectiousness_I2_wrt_I1 = gsl_matrix_calloc(src_rmp.l_relative_infectiousness_I2_wrt_I1->size1, src_rmp.l_relative_infectiousness_I2_wrt_I1->size2);
  dest_rmp.l_lbeta_rw = gsl_vector_calloc(src_rmp.l_lbeta_rw->size);
  dest_rmp.l_pr_symp = gsl_matrix_calloc(src_rmp.l_pr_symp->size1, src_rmp.l_pr_symp->size2);
  dest_rmp.l_pr_onset_to_GP = gsl_matrix_calloc(src_rmp.l_pr_onset_to_GP->size1, src_rmp.l_pr_onset_to_GP->size2);
  dest_rmp.l_pr_onset_to_Hosp = gsl_matrix_calloc(src_rmp.l_pr_onset_to_Hosp->size1, src_rmp.l_pr_onset_to_Hosp->size2);
  dest_rmp.l_pr_onset_to_Death = gsl_matrix_calloc(src_rmp.l_pr_onset_to_Death->size1, src_rmp.l_pr_onset_to_Death->size2);
  dest_rmp.l_importation_rate = gsl_matrix_calloc(src_rmp.l_importation_rate->size1, src_rmp.l_importation_rate->size2);
  dest_rmp.d_R0_phase_differences = gsl_vector_calloc(src_rmp.d_R0_phase_differences->size);

  mixing_model_alloc(dest_rmp.l_MIXMOD,
		     src_rmp.l_MIXMOD.num_breakpoints,
		     src_rmp.l_MIXMOD.scalants->size,
		     src_rmp.l_MIXMOD.evector_MIXMAT_normalised[0]->size);
  dest_rmp.l_background_gps_counts = gsl_matrix_calloc(src_rmp.l_background_gps_counts->size1, src_rmp.l_background_gps_counts->size2);
  dest_rmp.l_gp_negbin_overdispersion = gsl_matrix_calloc(src_rmp.l_gp_negbin_overdispersion->size1, src_rmp.l_gp_negbin_overdispersion->size2);
  dest_rmp.l_hosp_negbin_overdispersion = gsl_matrix_calloc(src_rmp.l_hosp_negbin_overdispersion->size1, src_rmp.l_hosp_negbin_overdispersion->size2);
  
  dest_rmp.l_day_of_week_effect = gsl_matrix_calloc(src_rmp.l_day_of_week_effect->size1, src_rmp.l_day_of_week_effect->size2);
}
void regional_model_params_free(regional_model_params& old_rmp)
{

  gsl_vector_free(old_rmp.l_init_prop_sus);
  gsl_vector_free(old_rmp.l_init_prop_sus_HI_geq_32);
  gsl_matrix_free(old_rmp.l_average_infectious_period);
  gsl_matrix_free(old_rmp.l_latent_period);
  gsl_matrix_free(old_rmp.l_relative_infectiousness_I2_wrt_I1);
  gsl_vector_free(old_rmp.l_lbeta_rw);
  gsl_matrix_free(old_rmp.l_pr_symp);
  gsl_matrix_free(old_rmp.l_pr_onset_to_GP);
  gsl_matrix_free(old_rmp.l_pr_onset_to_Hosp);
  gsl_matrix_free(old_rmp.l_pr_onset_to_Death);
  gsl_matrix_free(old_rmp.l_importation_rate);
  gsl_vector_free(old_rmp.d_R0_phase_differences);
  mixing_model_free(old_rmp.l_MIXMOD);
  gsl_matrix_free(old_rmp.l_background_gps_counts);
  gsl_matrix_free(old_rmp.l_gp_negbin_overdispersion);
  gsl_matrix_free(old_rmp.l_hosp_negbin_overdispersion);
  gsl_matrix_free(old_rmp.l_day_of_week_effect);
}

void regional_model_params_memcpy(regional_model_params& rmp_dest, const regional_model_params& rmp_src, flagclass& update_flags)
{
  if(update_flags.getFlag("l_init_prop_sus"))
    gsl_vector_memcpy(rmp_dest.l_init_prop_sus, rmp_src.l_init_prop_sus);
  if(update_flags.getFlag("l_init_prop_sus_HI_geq_32"))
    gsl_vector_memcpy(rmp_dest.l_init_prop_sus_HI_geq_32, rmp_src.l_init_prop_sus_HI_geq_32);
  if(update_flags.getFlag("l_average_infectious_period"))
    gsl_matrix_memcpy(rmp_dest.l_average_infectious_period, rmp_src.l_average_infectious_period);
  if(update_flags.getFlag("l_latent_period"))
    gsl_matrix_memcpy(rmp_dest.l_latent_period, rmp_src.l_latent_period);
  if(update_flags.getFlag("l_relative_infectious_period"))
    gsl_matrix_memcpy(rmp_dest.l_relative_infectiousness_I2_wrt_I1, rmp_src.l_relative_infectiousness_I2_wrt_I1);
  if(update_flags.getFlag("l_lbeta_rw"))
    gsl_vector_memcpy(rmp_dest.l_lbeta_rw, rmp_src.l_lbeta_rw);
  if(update_flags.getFlag("l_EGR"))
    rmp_dest.l_EGR = rmp_src.l_EGR;
  if(update_flags.getFlag("l_R0_Amplitude"))
    rmp_dest.l_R0_Amplitude = rmp_src.l_R0_Amplitude;
  if(update_flags.getFlag("l_R0_peakday"))
    rmp_dest.l_R0_peakday = rmp_src.l_R0_peakday;
  if(update_flags.getFlag("l_R0_init"))
    rmp_dest.l_R0_init = rmp_src.l_R0_init;
  if(update_flags.getFlag("l_I0"))
    rmp_dest.l_I0 = rmp_src.l_I0;
  if(update_flags.getFlag("l_pr_symp"))
    gsl_matrix_memcpy(rmp_dest.l_pr_symp, rmp_src.l_pr_symp);
  if(update_flags.getFlag("l_pr_onset_to_GP"))
    gsl_matrix_memcpy(rmp_dest.l_pr_onset_to_GP, rmp_src.l_pr_onset_to_GP);
  if(update_flags.getFlag("l_pr_onset_to_Hosp"))
    gsl_matrix_memcpy(rmp_dest.l_pr_onset_to_Hosp, rmp_src.l_pr_onset_to_Hosp);
  if(update_flags.getFlag("l_pr_onset_to_Death"))
    gsl_matrix_memcpy(rmp_dest.l_pr_onset_to_Death, rmp_src.l_pr_onset_to_Death);
  if(update_flags.getFlag("l_importation_rate"))
    gsl_matrix_memcpy(rmp_dest.l_importation_rate, rmp_src.l_importation_rate);
  if(update_flags.getFlag("d_R0_phase_differences"))
    gsl_vector_memcpy(rmp_dest.d_R0_phase_differences, rmp_src.d_R0_phase_differences);
  if(update_flags.getFlag("l_MIXMOD"))
    mixing_model_memcpy(rmp_dest.l_MIXMOD, rmp_src.l_MIXMOD);
  if(update_flags.getFlag("l_background_gps_counts"))
    gsl_matrix_memcpy(rmp_dest.l_background_gps_counts, rmp_src.l_background_gps_counts);
  if(update_flags.getFlag("l_sensitivity"))
    rmp_dest.l_sensitivity = rmp_src.l_sensitivity;
  if(update_flags.getFlag("l_specificity"))
    rmp_dest.l_specificity = rmp_src.l_specificity;
  if(update_flags.getFlag("l_sero_sensitivity"))
    rmp_dest.l_sero_sensitivity = rmp_src.l_sero_sensitivity;
  if(update_flags.getFlag("l_sero_specificity"))
    rmp_dest.l_sero_specificity = rmp_src.l_sero_specificity;
  if(update_flags.getFlag("l_gp_negbin_overdispersion"))
    gsl_matrix_memcpy(rmp_dest.l_gp_negbin_overdispersion, rmp_src.l_gp_negbin_overdispersion);
  if(update_flags.getFlag("l_hosp_negbin_overdispersion"))
    gsl_matrix_memcpy(rmp_dest.l_hosp_negbin_overdispersion, rmp_src.l_hosp_negbin_overdispersion);
  if(update_flags.getFlag("l_day_of_week_effect"))
    gsl_matrix_memcpy(rmp_dest.l_day_of_week_effect, rmp_src.l_day_of_week_effect);
}

void model_statistics_alloc(model_statistics &ms, const int times, const int age_classes, const int reporting_timesteps)
{
  ms.d_end_state = new model_state(age_classes);
  ms.d_NNI = gsl_matrix_alloc(times * reporting_timesteps, age_classes);
  ms.d_H1N1_GP_Consultations = gsl_matrix_alloc(times, age_classes);
  ms.d_Reported_GP_Consultations = gsl_matrix_alloc(times, age_classes);
  ms.d_Reported_Hospitalisations = gsl_matrix_alloc(times, age_classes);
  ms.d_seropositivity = gsl_matrix_alloc(times, age_classes);
  ms.d_viropositivity = gsl_matrix_alloc(times, age_classes);
}

void model_statistics_aggregate(gsl_matrix* output_NNI, const model_statistics& ms, const int contraction)
{
  output_per_selected_period(contraction, ms.d_NNI, output_NNI);
}

void model_statistics_memcpy(model_statistics &ms_dest, const model_statistics ms_src,
			     bool NNI_flag, bool GP_flag, bool Hosp_flag, bool Sero_flag, bool Viro_flag)
{
  *ms_dest.d_end_state = *ms_src.d_end_state;
if(NNI_flag)
    gsl_matrix_memcpy(ms_dest.d_NNI, ms_src.d_NNI);
  if(GP_flag)
    {
      gsl_matrix_memcpy(ms_dest.d_Reported_GP_Consultations, ms_src.d_Reported_GP_Consultations);
      gsl_matrix_memcpy(ms_dest.d_H1N1_GP_Consultations, ms_src.d_H1N1_GP_Consultations);
    }
  if(Hosp_flag)
    gsl_matrix_memcpy(ms_dest.d_Reported_Hospitalisations, ms_src.d_Reported_Hospitalisations);
  if(Sero_flag)
    gsl_matrix_memcpy(ms_dest.d_seropositivity, ms_src.d_seropositivity);
  if(Viro_flag)
    gsl_matrix_memcpy(ms_dest.d_viropositivity, ms_src.d_viropositivity);
}

void model_statistics_free(struct model_statistics &ms)
{
  delete ms.d_end_state;
  gsl_matrix_free(ms.d_NNI);
  gsl_matrix_free(ms.d_Reported_GP_Consultations);
  gsl_matrix_free(ms.d_H1N1_GP_Consultations);
  gsl_matrix_free(ms.d_Reported_Hospitalisations);
  gsl_matrix_free(ms.d_seropositivity);
  gsl_matrix_free(ms.d_viropositivity);
}

// ---- Overloaded region allocation function
void Region_alloc(Region& new_reg,
		  const global_model_instance_parameters src_gmip,
		  const mixing_model src_mixmod)
{
  int num_ages = NUM_AGE_GROUPS;
  int num_days = src_gmip.l_duration_of_runs_in_days;

  new_reg.data_owner = true;
  new_reg.population = gsl_vector_calloc(num_ages);
  if(src_gmip.l_GP_consultation_flag) new_reg.GP_data = new rtmData(src_gmip.l_GP_likelihood, src_gmip.l_gp_count_likelihood);
  else new_reg.GP_data = 0;
  if(src_gmip.l_Hospitalisation_flag) new_reg.Hospitalisation_data = new rtmData(src_gmip.l_Hosp_likelihood, src_gmip.l_hosp_count_likelihood);
  else new_reg.Hospitalisation_data = 0;
  if(src_gmip.l_Deaths_flag) new_reg.Death_data = new rtmData(src_gmip.l_Deaths_likelihood, cPOISSON);
  else new_reg.Death_data = 0;
  if(src_gmip.l_Sero_data_flag) new_reg.Serology_data = new rtmData(src_gmip.l_Sero_likelihood, cBINOMIAL);
  else new_reg.Serology_data = 0;
  if(src_gmip.l_Viro_data_flag) new_reg.Virology_data = new rtmData(src_gmip.l_Viro_likelihood, cBINOMIAL);
  else new_reg.Virology_data = 0;

  // FUNCTIONS ALLOCATING THE regional_model_params AND model_statistics STRUCTURES NEEDED HERE
  regional_model_params_alloc(new_reg.det_model_params, num_days, num_ages, src_gmip.l_transmission_time_steps_per_day, src_gmip.l_reporting_time_steps_per_day, src_mixmod);
  model_statistics_alloc(new_reg.region_modstats, num_days, num_ages, src_gmip.l_reporting_time_steps_per_day);
}

// ---- alternate form of the above
void Region_alloc(Region& new_reg,
		  const Region& old_reg)
{
  new_reg.data_owner = false;
  regional_model_params_alloc(new_reg.det_model_params, old_reg.det_model_params);
  model_statistics_alloc(new_reg.region_modstats,
			 old_reg.region_modstats.d_seropositivity->size1,
			 old_reg.region_modstats.d_seropositivity->size2,
			 old_reg.region_modstats.d_NNI->size1 / old_reg.region_modstats.d_seropositivity->size1);
}
void Region_free(Region& old_reg, const global_model_instance_parameters src_gmip)
{

  if(old_reg.data_owner){
    gsl_vector_free(old_reg.population);
    if(src_gmip.l_GP_consultation_flag) delete old_reg.GP_data;
    if(src_gmip.l_Hospitalisation_flag) delete old_reg.Hospitalisation_data;
    if(src_gmip.l_Deaths_flag) delete old_reg.Death_data;
    if(src_gmip.l_Sero_data_flag) delete old_reg.Serology_data;
    if(src_gmip.l_Viro_data_flag) delete old_reg.Virology_data;
  }
  // FUNCTIONS FREEING THE regional_model_params AND model_statistics STRUCTURES NEEDED HERE
  regional_model_params_free(old_reg.det_model_params);
  model_statistics_free(old_reg.region_modstats);
}

void Region_memcpy(Region& reg_dest, const Region& reg_src, flagclass& update_flags)
{
  reg_dest.name.assign(reg_src.name);
  reg_dest.population = reg_src.population;
  reg_dest.total_population = reg_src.total_population;
  reg_dest.data_owner = false;
  if(reg_src.GP_data != 0) // Datasets should be invariant, so ok to equate pointers.
      reg_dest.GP_data = reg_src.GP_data;
  if(reg_src.Hospitalisation_data != 0)
    reg_dest.Hospitalisation_data = reg_src.Hospitalisation_data;
  if(reg_src.Death_data != 0)
      reg_dest.Death_data = reg_src.Death_data;
  if(reg_src.Serology_data != 0)
      reg_dest.Serology_data = reg_src.Serology_data;
  if(reg_src.Virology_data != 0)
      reg_dest.Virology_data = reg_src.Virology_data;


  regional_model_params_memcpy(reg_dest.det_model_params, reg_src.det_model_params, update_flags);
  model_statistics_memcpy(reg_dest.region_modstats, reg_src.region_modstats);
}
