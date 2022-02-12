#ifndef HEADER_StructAssign_
#define HEADER_StructAssign_

#define GLOBAL_MODEL_INSTANCE_MEMBERS "gp_count_likelihood : hosp_count_likelihood : transmission_kernel : transmission_time_steps_per_day : reporting_time_steps_per_day : duration_of_runs_in_days : projection_days : day_of_start : start_day_of_week : start_week : num_regions : GP_consultation_flag : Hospitalisation_flag : Deaths_flag : Sero_data_flag : Viro_data_flag : Prev_data_flag : Vacc_data_flag : VBoost_data_flag : Use_GP_reparam_patch : Sero_delay : GP_likelihood_start_day : GP_likelihood_end_day : Hosp_likelihood_start_day : Hosp_likelihood_end_day : Deaths_likelihood_start_day : Deaths_likelihood_end_day : Sero_likelihood_start_day : Sero_likelihood_end_day : Viro_likelihood_start_day : Viro_likelihood_end_day : Prev_likelihood_start_day : Prev_likelihood_end_day : Vacc_start_day : Vacc_end_day : VBoost_start_day : VBoost_end_day"



#define GLOBAL_MODEL_INSTANCE_DEFAULTS "0 : 0 : 0 : 4 : 1 : 192 : 250 : 121 : 5 : 18 : 1 : 1 : 1 : 0 : 1 : 0 : 0 : 0 : 0 : 1 : 14 : 1 : 192 : 1 : 50 : 1 : 1 : 124 : 124 : 1 : 192 : 1 : 192 : 1 : 192 : 1 : 192"



// UPDATEABLE PARAMETERS

typedef enum { EGR_HYPER_INDEX, LPL0_HYPER_INDEX, PROP_SUS_HYPER_INDEX, LBETA_RW_SD_INDEX, GP_OVERDISP_INDEX, HOSP_OVERDISP_INDEX, ALP_INDEX, AIP_INDEX, AR1_INDEX, VAC1_DISEASE_INDEX, VACN_DISEASE_INDEX, VACB_DISEASE_INDEX, VAC1_INFECT_INDEX, VACN_INFECT_INDEX, VACB_INFECT_INDEX, REL_INFECT_INDEX, PROP_SYMP_INDEX, CONTACT_INDEX, LBETA_RW_INDEX, R0_AMP_INDEX, R0_PEAKDAY_INDEX, EGR_INDEX, LPL0_INDEX, PROP_SUS_INDEX, PROP_HI_GEQ_32_INDEX, PROP_GP_INDEX, PROP_HOSP_INDEX, PROP_DEATH_INDEX, IMPORTATION_INDEX, BGR_INDEX, SENS_INDEX, SPEC_INDEX, DOW_EFFECTS_INDEX, SSENS_INDEX, SSPEC_INDEX, IWAN_INDEX } updateable_parameter_index;

#define GLOBAL_MODEL_PARAMETERS_MEMBERS "exponential_growth_rate_hyper : l_p_lambda_0_hyper : prop_susceptible_hyper : log_beta_rw_sd : gp_negbin_overdispersion : hosp_negbin_overdispersion : latent_period : infectious_period : r1_period : vacc_1st_disease : vacc_nth_disease : vacc_boost_disease : vacc_1st_infect : vacc_nth_infect : vacc_boost_infect : relative_infectiousness : prop_symptomatic : contact_parameters : log_beta_rw : R0_amplitude_kA : R0_seasonal_peakday : exponential_growth_rate : log_p_lambda_0 : prop_susceptible : prop_HI_32_to_HI_8 : prop_case_to_GP_consultation : prop_case_to_hosp : prop_case_to_death : importation_rates : background_GP : test_sensitivity : test_specificity : day_of_week_effects : sero_test_sensitivity : sero_test_specificity : immunity_period"

#define GLOBAL_MODEL_PARAMETERS_DEFAULT_FIXED_VALS "0.15 : -15.0 : 0.8 : 0.05 : 1.0 : 0.01 : 2.0 : 2.5 : 1.5 : 0.8 : 0.95 : 0.99 : 0.8 : 0.95 : 0.99 : 1.0 : 0.45 : 1.0 : 0.0 : 0.0 : 355 : 0.15 : -15.0 : 0.8 : 0.5 : 0.0 : 0.1 : 0.001 : 0.0 : 5.7 : 0.95 : 0.95 : 1.0 : 1.0 : 1.0 : 0.002"

#define GLOBAL_MODEL_PARAMETERS_LIKELIHOOD_FLAGS "0000000 : 0000000 : 0000000 : 0000000 : 0010000 : 0001000 : 1111111 : 1111111 : 1000001 : 1111111 : 1111111 : 1111111 : 1111111 : 1111111 : 1111111 : 0111010 : 1111111 : 1111111 : 1111111 : 1111111 : 1111111 : 1111111 : 1111111 : 0000100 : 1111111 : 0101000 : 0100000 : 1111111 : 0010010 : 0000010 : 0000010 : 0010000 : 0000100 : 0000100 : 1111111" // FLAGS CORRESPOND TO... TRANSMISSION MODEL, REPORTING MODEL, GP LIKELIHOOD, HOSP LIKELIHOOD, SEROLOGICAL DATA LIKELIHOOD, VIROLOGICAL DATA LIKELIHOOD, PREV_DATA_LIKELIHOOD



// NON-UPDATEABLE PARAMETERS OF THE DELAY DISTRIBUTIONS

#define GLOBAL_MODEL_PARAMETERS_DELAY_NAMES "incubation_distribution : symp_to_gp_distribution : symp_to_hosp_distribution : symp_to_death_distribution : gp_reporting_distribution : hosp_reporting_distribution : death_reporting_distribution"

#define GLOBAL_MODEL_PARAMETERS_DELAY_MEANS " 1.620 : 1.968 : 0.0 : 10.75 : 0.5 : 6.577 : 75.89"

#define GLOBAL_MODEL_PARAMETERS_DELAY_SDS " 1.791 : 1.199 : 0.0 : 8.336 : 0.5 : 3.722 : 133.8"

#define GLOBAL_MODEL_PARAMETERS_DELAY_FLAGS "111 : 100 : 010 : 001 : 100 : 010 : 001"


// PARAMETERS OF THE MCMC SIMULATION

#define MCMC_PARAMETER_NAMES " num_iterations : output_type : thin_output_every : thin_stats_every : adaptive_phase : adapt_every : burn_in : read_covar : global_updates : num_progress_reports : mixing_threshold_ub : mixing_threshold_lb : prop_var_inflation : prop_var_shrinkage : maximum_block_size : random_seed : max_threads"

#define MCMC_PARAMETER_DEFAULT_VALUES " 1100000 : 0 : 100 : 10000 : 50000 : 100 : 100000 : 0 : 1 : 10 : 0.3 : 0.2 : 1.1 : 0.9 : 32 : 32 : 1"

#endif
