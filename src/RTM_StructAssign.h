#ifndef HEADER_StructAssign_
#define HEADER_StructAssign_

#define GLOBAL_MODEL_INSTANCE_MEMBERS "gp_count_likelihood : hosp_count_likelihood : transmission_kernel : transmission_time_steps_per_day : reporting_time_steps_per_day : duration_of_runs_in_days : projection_days : day_of_start : start_day_of_week : start_week : num_regions : GP_consultation_flag : Hospitalisation_flag : Deaths_flag : Sero_data_flag : Viro_data_flag : Use_GP_reparam_patch : Sero_delay : GP_likelihood_start_day : GP_likelihood_end_day : Hosp_likelihood_start_day : Hosp_likelihood_end_day : Deaths_likelihood_start_day : Deaths_likelihood_end_day : Sero_likelihood_start_day : Sero_likelihood_end_day : Viro_likelihood_start_day : Viro_likelihood_end_day"



#define GLOBAL_MODEL_INSTANCE_DEFAULTS "0 : 0 : 0 : 4 : 1 : 192 : 250 : 121 : 5 : 18 : 1 : 1 : 1 : 0 : 1 : 0 : 1 : 14 : 1 : 192 : 1 : 50 : 1 : 1 : 124 : 124 : 1 : 192"



// UPDATEABLE PARAMETERS

typedef enum { EGR_HYPER_INDEX, LPL0_HYPER_INDEX, PROP_SUS_HYPER_INDEX, GP_OVERDISP_INDEX, HOSP_OVERDISP_INDEX, ALP_INDEX, AIP_INDEX, REL_INFECT_INDEX, PROP_SYMP_INDEX, CONTACT_INDEX, R0_AMP_INDEX, R0_PEAKDAY_INDEX, EGR_INDEX, LPL0_INDEX, PROP_SUS_INDEX, PROP_HI_GEQ_32_INDEX, PROP_GP_INDEX, PROP_HOSP_INDEX, PROP_DEATH_INDEX, IMPORTATION_INDEX, BGR_INDEX, SENS_INDEX, SPEC_INDEX, DOW_EFFECTS_INDEX } updateable_parameter_index;

#define GLOBAL_MODEL_PARAMETERS_MEMBERS "exponential_growth_rate_hyper : l_p_lambda_0_hyper : prop_susceptible_hyper : gp_negbin_overdispersion : hosp_negbin_overdispersion : latent_period : infectious_period : relative_infectiousness : prop_symptomatic : contact_parameters : R0_amplitude_kA : R0_seasonal_peakday : exponential_growth_rate : log_p_lambda_0 : prop_susceptible : prop_HI_32_to_HI_8 : prop_case_to_GP_consultation : prop_case_to_hosp : prop_case_to_death : importation_rates : background_GP : test_sensitivity : test_specificity : day_of_week_effects"

#define GLOBAL_MODEL_PARAMETERS_DEFAULT_FIXED_VALS "0.15 : -15.0 : 0.8 : 1.0 : 0.01 : 2.0 : 2.5 : 1.0 : 0.45 : 1.0 : 0.0 : 355 : 0.15 : -15.0 : 0.8 : 0.5 : 0.0 : 0.1 : 0.001 : 0.0 : 5.7 : 0.95 : 0.95 : 1.0"

#define GLOBAL_MODEL_PARAMETERS_LIKELIHOOD_FLAGS "000000 : 000000 : 000000 : 001000 : 000100 : 111111 : 111111 : 111111 : 011101 : 111111 : 111111 : 111111 : 111111 : 111111 : 111111 : 000010 : 111111 : 010100 : 010000 : 111111 : 001001 : 000001 : 000001 : 001000 " // FLAGS CORRESPOND TO... TRANSMISSION MODEL, REPORTING MODEL, GP LIKELIHOOD, HOSP LIKELIHOOD, SEROLOGICAL DATA LIKELIHOOD, VIROLOGICAL DATA LIKELIHOOD



// NON-UPDATEABLE PARAMETERS OF THE DELAY DISTRIBUTIONS

#define GLOBAL_MODEL_PARAMETERS_DELAY_NAMES "incubation_distribution : symp_to_gp_distribution : symp_to_hosp_distribution : symp_to_death_distribution : gp_reporting_distribution : hosp_reporting_distribution : death_reporting_distribution"

#define GLOBAL_MODEL_PARAMETERS_DELAY_MEANS " 1.620 : 1.968 : 0.0 : 10.75 : 0.5 : 6.577 : 75.89"

#define GLOBAL_MODEL_PARAMETERS_DELAY_SDS " 1.791 : 1.199 : 0.0 : 8.336 : 0.5 : 3.722 : 133.8"

#define GLOBAL_MODEL_PARAMETERS_DELAY_FLAGS "111 : 100 : 010 : 001 : 100 : 010 : 001"


// PARAMETERS OF THE MCMC SIMULATION

#define MCMC_PARAMETER_NAMES " num_iterations : output_type : thin_output_every : thin_stats_every : adaptive_phase : adapt_every : burn_in : num_progress_reports : mixing_threshold_ub : mixing_threshold_lb : prop_var_inflation : prop_var_shrinkage : maximum_block_size : random_seed : max_threads"

#define MCMC_PARAMETER_DEFAULT_VALUES " 1100000 : 0 : 100 : 10000 : 50000 : 100 : 100000 : 10 : 0.3 : 0.2 : 1.1 : 0.9 : 32 : 32 : 1"

#endif
