#include "RTM_StructDefs.h"
#include "RTM_FunctDefs.h"
#include "RTM_flagclass.h"
#include "gsl_vec_ext.h"
#include "string_fns.h"

using namespace std;
using std::string;

#define OUTPUT_VECTOR_DELIMITER(i, j) (((i) + 1) == simulation_parameters.num_iterations) && (((j) + 1) == theta_i->param_value->size) ? "" : \
  ((((j) + 1) == theta_i->param_value->size) ? ",\n" : ", ")
#define OUTPUT_MATRIX_DELIMITER(i, j, k)  (((i) + 1 == simulation_parameters.num_iterations) && ((j) + 1 == NUM_AGE_GROUPS) && ((k) + 1 == gmip.l_duration_of_runs_in_days)) ? "" : (((j) + 1 == NUM_AGE_GROUPS) && ((k) + 1 == gmip.l_duration_of_runs_in_days) ? ",\n\n" : (((j) + 1 == NUM_AGE_GROUPS) ? ",\n" : ", "))

#define CHAIN_LENGTH (simulation_parameters.num_iterations - simulation_parameters.burn_in)

// //
// INITIALISES model_state OUTPUT FILES
void fstream_model_statistics_open(ofstream *&infiles, const string& filenames, const int& num_files, const Region* reg_name)
{

  infiles = new ofstream[num_files];
  for(int int_reg = 0; int_reg < num_files; int_reg++)
    {
      string full_filename = filenames + "_" + reg_name[int_reg].name;
      infiles[int_reg].open(full_filename.c_str(), ios::out|ios::trunc|ios::binary);
    }
}

// TERMINATE AND CLOSE model_statistics OUTPUT_FILES, fstream VERSION
void fstream_model_statistics_close(ofstream *&infiles, const int& nregions, const int& dim_ages, const int&dim_days, const int& dim_iters)
{
  for(int int_reg = 0; int_reg < nregions; int_reg++)
    infiles[int_reg] << dim_ages << " " << dim_days << " " << dim_iters << endl;

  delete [] infiles;
}
// //

// CALCULATES A RATIO OF THE PRIOR EVALUATED AT THE CANDIDATE PROPOSAL AND AT THE CURRENT STATE OF THE CHAIN //
double univariate_prior_ratio(const updateable_model_parameter& theta,
		   const int& component)
{

  int dist_component = (theta.flag_hyperprior) ? 1 : component;

  return R_univariate_prior_log_density_ratio(gsl_vector_get(theta.proposal_value, component),
					      gsl_vector_get(theta.param_value, component),
					      (distribution_type) gsl_vector_int_get(theta.prior_distribution, dist_component),
					      theta.prior_params[component]);

}
// //


// ADAPT PROPOSAL VARIANCES //
void adapt_proposal_variance(updateable_model_parameter& ump, const mcmcPars mcmc_pars)
{

  for(int int_i = 0; int_i < ump.param_value->size; int_i++)
    {

      int distribution_index = (ump.flag_hyperprior) ? 0 : int_i;

      if(gsl_vector_int_get(ump.prior_distribution, distribution_index) != (int) cCONSTANT){
	double acceptance_ratio = ((double) gsl_vector_int_get(ump.number_accepted_moves, int_i)) / ((double) gsl_vector_int_get(ump.number_proposed_moves, int_i));

	if(acceptance_ratio < mcmc_pars.mixing_threshold_lb) // NOT ENOUGH JUMPS, DECREASE PROPOSAL VARIANCE
	  gsl_vector_set(ump.proposal_variances, int_i, gsl_vector_get(ump.proposal_variances, int_i) * mcmc_pars.prop_var_shrinkage);
	else if(acceptance_ratio > mcmc_pars.mixing_threshold_ub){ // TOO MANY SMALL JUMPS, INCREASE PROPOSAL VARIANCE IF PROPOSAL VARIANCE NOT ALREADY TOO LARGE
	  double dbl_var_adj = variance_inflation_factor(gsl_vector_get(ump.proposal_variances, int_i) * mcmc_pars.prop_var_inflation,
							 (distribution_type) gsl_vector_int_get(ump.prior_distribution, int_i),
							 ump.prior_params[int_i]);

	  gsl_vector_set(ump.proposal_variances, int_i, dbl_var_adj);

	}
      }
    }

}
// //

// RESET THE ACCEPTANCE RATIO COUNTERS //
void reset_counters(globalModelParams& theta)
{
  for(int int_param = 0; int_param < theta.size_param_list; int_param++)
    if(theta.param_list[int_param].flag_update){
      gsl_vector_int_set_zero(theta.param_list[int_param].number_accepted_moves);
      gsl_vector_int_set_zero(theta.param_list[int_param].number_proposed_moves);
    }
}

//// FUNCTIONS TO PRODUCE PROGRESS REPORTS ////

// TRANSFORM DYNAMICALLY UPDATED POSTERIOR MEAN AND POSTERIOR MEAN SUM OF SQUARES TO POSTERIOR STANDARD DEVIATION
void fn_get_sd_from_sum_squares(gsl_vector* sd, const gsl_vector* mean, const gsl_vector* mean_ss)
{

  gsl_vector_memcpy(sd, mean);
  gsl_vector_mul(sd, mean); // sd contains xbar^2
  gsl_vector_sub(sd, mean_ss);
  gsl_vector_scale(sd, -1.0);

  for(int i = 0; i < mean->size; i++)
    gsl_vector_set(sd, i, sqrt(gsl_vector_get(sd, i)));

}

void print_param_name(FILE* ifile, const string& param_name)
{
  fprintf(ifile, "\n%s\n", param_name.c_str());
  for(int i = 0; i++ < param_name.length(); )
    fprintf(ifile, "=");
  fprintf(ifile, "\n");
}


// ALL-PURPOSE REPORT FUNCTION //
void write_progress_report(const string& report_type, const int& int_report_no, const int& int_iter, const int& chain_length,
			   const globalModelParams& gmp, const likelihood& llhood,
			   bool posterior_stats_flag, bool acceptance_flag, bool propvar_flag)
{
  string filename(report_type);
  filename += int_to_string(int_report_no);
  FILE* ofile = fopen(filename.c_str(), "w");
  double scale_val = ((double) chain_length) / ((double) int_iter);

  for(int int_param = 0; int_param < gmp.size_param_list; int_param++)
    {
      if(gmp.param_list[int_param].flag_update)
	{
	  updateable_model_parameter theta = gmp.param_list[int_param];

	  print_param_name(ofile, theta.param_name);
	  
	  if(posterior_stats_flag)
	    {

	      gsl_vector* posterior_mean = gsl_vector_alloc(theta.posterior_mean->size);
	      gsl_vector* posterior_sumsq = gsl_vector_alloc(theta.posterior_mean->size);
	      gsl_vector* posterior_sd = gsl_vector_alloc(theta.posterior_mean->size);
	      gsl_vector_memcpy(posterior_mean, theta.posterior_mean);
	      gsl_vector_memcpy(posterior_sumsq, theta.posterior_sumsq);
	      gsl_vector_scale(posterior_mean, scale_val);
	      gsl_vector_scale(posterior_sumsq, scale_val);
	      fn_get_sd_from_sum_squares(posterior_sd, posterior_mean, posterior_sumsq);
	      for(int int_i = 0; int_i < posterior_sd->size; int_i++)
		fprintf(ofile, "posterior_mean%d = %f (%f)\n", int_i, gsl_vector_get(posterior_mean, int_i), gsl_vector_get(posterior_sd, int_i));
	      gsl_vector_free(posterior_mean);
	      gsl_vector_free(posterior_sumsq);
	      gsl_vector_free(posterior_sd);
	    }	 

	  if(acceptance_flag)
	    {
	      for(int int_i = 0; int_i < theta.param_value->size; int_i++)
		fprintf(ofile, "acceptance_ratio%d = %f\n", int_i + 1, ((double) gsl_vector_int_get(theta.number_accepted_moves, int_i)) / ((double) gsl_vector_int_get(theta.number_proposed_moves, int_i))); 
	      fprintf(ofile, "\n");
	    }

	  if(propvar_flag)
	    {
	      for(int int_i = 0; int_i < theta.param_value->size; int_i++)
		fprintf(ofile, "proposal_variance%d = %f\n", int_i + 1, gsl_vector_get(theta.proposal_variances, int_i));
	      fprintf(ofile, "\n");
	    }
	}

    }
  // IF PROVIDING A POSTERIOR SUMMARY, ALSO ADD IN SOMETHING FOR THE LIKELIHOOD
  if(posterior_stats_flag)
    {
      double temp_mean = llhood.bar_lfx, temp_sumsq = llhood.sumsq_lfx, temp_sd;
      temp_mean *= scale_val;
      temp_sumsq *= scale_val;
      temp_sd = temp_sumsq - gsl_pow_2(temp_mean);
      print_param_name(ofile, "likelihood");
      fprintf(ofile, "lfx_bar = %f (%f)\n", temp_mean, sqrt(temp_sd));
    }
    
  fclose(ofile);
}
//// ////







// FUNCTION TO BE CALLED FROM THE MAIN ROUTINE //
void metrop_hast(const mcmcPars& simulation_parameters,
		 globalModelParams& theta,
		 Region* state_country,
		 likelihood& lfx,
		 const global_model_instance_parameters& gmip,
		 const mixing_model& base_mix,
		 gsl_rng* r)
{

  register int int_iter = 0;
  register int int_param = 0;
  gsl_vector_int* num_component_updates = gsl_vector_int_alloc(theta.size_param_list);
  gsl_vector_int* block_size = gsl_vector_int_alloc(theta.size_param_list);
  gsl_vector_int **a = new gsl_vector_int*[theta.size_param_list];
  gsl_vector_int **b = new gsl_vector_int*[theta.size_param_list];
  gsl_vector_int* adaptive_progress_report_iterations = gsl_vector_int_alloc(simulation_parameters.num_progress_reports);
  gsl_vector_int* chain_progress_report_iterations = gsl_vector_int_alloc(simulation_parameters.num_progress_reports);
  int int_progress_report = 0;
  int nregions = gmip.l_num_regions;

  // TAKE A COPY OF THE REGION STRUCTURE AND THE LIKELIHOOD STRUCTURE
  flagclass all_pos;

  // Region
  Region* prop_country = new Region[nregions];
  for(int int_reg = 0; int_reg < nregions; ++int_reg)
    {
      Region_alloc(prop_country[int_reg], state_country[int_reg]);

      Region_memcpy(prop_country[int_reg], state_country[int_reg], all_pos);
    }

  // Likelihood
  likelihood prop_lfx;
  likelihood_alloc(prop_lfx, gmip);
  likelihood_memcpy(prop_lfx, lfx);
  // //

  int maximum_block_size = simulation_parameters.maximum_block_size;

  // Outputs
  ofstream* output_codas = new ofstream[num_component_updates->size];
  ofstream output_coda_lfx("coda_lfx", ios::out|ios::trunc|ios::binary);
  ofstream *file_NNI, *file_GP, *file_Hosp, *file_Sero, *file_Viro, *file_state;

  fstream_model_statistics_open(file_NNI, "NNI", gmip.l_num_regions, state_country);
  if(simulation_parameters.oType == cMCMC)
    {
      fstream_model_statistics_open(file_Sero, "Sero", gmip.l_num_regions, state_country);
      if(gmip.l_GP_consultation_flag){
	fstream_model_statistics_open(file_GP, "GP", gmip.l_num_regions, state_country);
	fstream_model_statistics_open(file_Viro, "Viro", gmip.l_num_regions, state_country);
      }
      if(gmip.l_Hospitalisation_flag) 
	fstream_model_statistics_open(file_Hosp, "Hosp", gmip.l_num_regions, state_country);
    }
  if(simulation_parameters.oType == cSMC)
    fstream_model_statistics_open(file_state, "state", gmip.l_num_regions, state_country);

  gsl_matrix** output_NNI = new gsl_matrix*[nregions];
  for(int int_reg = 0; int_reg < nregions; ++int_reg)
    output_NNI[int_reg] = gsl_matrix_alloc(gmip.l_duration_of_runs_in_days, NUM_AGE_GROUPS);

  gsl_vector_int_set_1ton(adaptive_progress_report_iterations, 1);
  gsl_vector_int_scale(adaptive_progress_report_iterations, simulation_parameters.adaptive_phase / simulation_parameters.num_progress_reports); // NOTE THE INTEGRAL DIVISION
  gsl_vector_int_set(adaptive_progress_report_iterations, adaptive_progress_report_iterations->size - 1, simulation_parameters.adaptive_phase); // Ensure the adaptive phase runs for the full specified number of iterations.
  gsl_vector_int_set_1ton(chain_progress_report_iterations, 1);
  gsl_vector_int_scale(chain_progress_report_iterations, CHAIN_LENGTH / simulation_parameters.num_progress_reports); // NOTE THE INTEGRAL DIVISION
  gsl_vector_int_add_constant(chain_progress_report_iterations, simulation_parameters.burn_in);
  // //

  // Determine block sizes for the various updates
  for(int_param = 0; int_param < num_component_updates->size; int_param++)
    {

      int theta_i_size = theta.param_list[int_param].param_value->size;

      gsl_vector_int_set(num_component_updates, int_param, theta_i_size / maximum_block_size); // NOTE: INTEGER DIVISION
      int int_temp = theta_i_size % maximum_block_size;
      gsl_vector_int_set(num_component_updates, int_param, gsl_vector_int_get(num_component_updates, int_param) + ((int_temp == 0) ? 0 : 1));

      // What is the size of the component blocks proposed
      gsl_vector_int_set(block_size, int_param, ((theta_i_size - 1) / gsl_vector_int_get(num_component_updates, int_param)) + 1); // NOTE: INTEGER DIVISION

      b[int_param] = gsl_vector_int_alloc(theta_i_size);
      a[int_param] = gsl_vector_int_alloc(gsl_vector_int_get(block_size, int_param));

      gsl_vector_int_set_1ton(b[int_param]);
      if(theta_i_size <= maximum_block_size)
	gsl_vector_int_set_1ton(a[int_param]);

      // OPEN PARAMETER OUTPUT FILES...
      string filename("coda_");
      filename += theta.param_list[int_param].param_name;
      if(theta.param_list[int_param].flag_update)
	{output_codas[int_param].open(filename.c_str(), ios::out|ios::trunc|ios::binary);}
    }

  // Central Loop //
  for(; int_iter < simulation_parameters.num_iterations; int_iter++)
    {

      // LOOP THROUGH EACH OF THE UPDATEABLE PARAMETERS
      for(int_param = 0; int_param < theta.size_param_list; int_param++)
	{

	  if(theta.param_list[int_param].flag_update)
	    {

	      // SET UP THE FLAGCLASS TO TELL evaluate_regional_parameters WHICH COMPONENTS NEED UPDATING
	      // SHOULD BE CONSTANT WHILST UPDATING THIS ONE PARAMETER GROUP
	      flagclass update_flags(int_param);

	      // KEEP AN EYE OUT FOR WHICH PARTS OF THIS LOOP CAN BE TAKEN OUTSIDE AND NOT REPEATED EVERY ITERATION

	      // ISOLATE THE PARAMETER WHOSE COMPONENTS WE WILL BE UPDATING
	      updateable_model_parameter *theta_i = theta.param_list + int_param;

	      // CARRY OUT BLOCK UPDATES ON THIS PARAMETER
	      for(int int_i = 0; int_i < gsl_vector_int_get(num_component_updates, int_param); int_i++)
		{

		  double log_accep = 0.0;

		  gsl_vector_memcpy(theta_i->proposal_value, theta_i->param_value);

		  // IS THE PARAMETER VECTOR TOO LARGE FOR A SINGLE BLOCK UPDATE?
		  // IF SO, SET a TO SOME RANDOM INDICES TO UPDATE
		  if(theta_i->param_value->size > maximum_block_size)
		    gsl_ran_choose(r, a[int_param]->data, gsl_vector_int_get(block_size, int_param), b[int_param]->data, theta_i->param_value->size, sizeof(int));
		  // ELSE a = (1, ..., theta_i.param_value->size) AND THIS SHOULD BE PRESET.
		  theta_i->proposal_log_prior_dens = 0;

		  for(int int_component = 0; int_component < gsl_vector_int_get(block_size, int_param); int_component++)
		    {

		      int a_component = gsl_vector_int_get(a[int_param], int_component);
		      // UPDATE THE PROPOSAL COUNTER FOR THIS INDEX
		      gsl_vector_int_set(theta_i->number_proposed_moves, a_component, gsl_vector_int_get(theta_i->number_proposed_moves, a_component) + 1);

		      if(gsl_vector_int_get(theta_i->prior_distribution, a_component) != cCONSTANT)
			{ // PROPOSE NEW VALUES - STORE RANDOM WALK PROPOSAL RATIO

			  double trunclo = (gsl_vector_int_get(theta_i->prior_distribution, a_component) == cUNIFORM) ? gsl_vector_get(theta_i->prior_params[a_component], 0) : 0.0;
			  double trunchi = (gsl_vector_int_get(theta_i->prior_distribution, a_component) == cUNIFORM) ? gsl_vector_get(theta_i->prior_params[a_component], 1) : 1.0;

			  log_accep += random_walk_proposal(*(gsl_vector_ptr(theta_i->proposal_value, a_component)), // IS THE *gsl_vector_ptr SYNTAX NECESSARY - COULD gsl_vector_get BE USED INSTEAD?
							    *(gsl_vector_ptr(theta_i->param_value, a_component)),
							    (distribution_type) gsl_vector_int_get(theta_i->prior_distribution, a_component),
							    gsl_vector_get(theta_i->proposal_variances, a_component),
							    r, trunclo, trunchi);

			  // EVALUATE THE LOG PRIOR DENSITY - CAN ONLY DO THIS HERE IF UNIVARIATE PRIOR IS SPECIFIED
			  if(log_accep > GSL_NEGINF){
			    if(gsl_vector_int_get(theta_i->prior_distribution, a_component) != (int) cMVNORMAL)
			      theta_i->proposal_log_prior_dens += univariate_prior_ratio(*theta_i, a_component);
			  }
			}

		    }

		  // EVALUATE THE LOG PRIOR DENSITY - MULTIVARIATE CASE
		  if(log_accep > GSL_NEGINF){
		    if(gsl_vector_int_get(theta_i->prior_distribution, 0) == (int) cMVNORMAL)
		      theta_i->proposal_log_prior_dens = (*theta_i->prior_multivariate_norm).ld_mvnorm_ratio(theta_i->proposal_value, theta_i->param_value);

		    log_accep += theta_i->proposal_log_prior_dens;

		    if(log_accep > GSL_NEGINF){
		      // IS THIS PARAMETER A MODEL PARAMETER OR A HYPERPARAMETER?
		      // HYPERPARAMETER IF IT HAS ANY CHILD NODES
		      if(theta_i->flag_any_child_nodes)
			{

			  // THE "LIKELIHOOD" COMES FROM THE PRIOR DENSITIES OF THE CHILD NODES
			  for(int int_child_node = 0; int_child_node < theta.size_param_list; int_child_node++)
			    {
			      if(int_child_node != int_param)
				if(theta_i->flag_child_nodes[int_child_node])
				  {
				    updateable_model_parameter* theta_child = theta.param_list + int_child_node;

				    theta_child->proposal_log_prior_dens = 0;

				    for(int int_j = 0; int_j < theta_child->param_value->size; int_j++)
				      theta_child->proposal_log_prior_dens += R_univariate_prior_log_density(gsl_vector_get(theta_child->param_value, int_j),
													     (distribution_type) gsl_vector_int_get(theta_child->prior_distribution, 0),
													     theta_i->proposal_value);

				    log_accep += theta_child->proposal_log_prior_dens - theta_child->log_prior_dens;
				  }

			    }
			}
		    
		      else
			{
		       
			  // ADJUST THE COPIED REGION STRUCTURES TO THE PROPOSED VALUE
			  // evaluate_regional_parameters work by using the param_value... need to temporarily switch the two
			  gsl_vector* tempvec = theta_i->param_value;
			  theta_i->param_value = theta_i->proposal_value; // temporarily switch the value of the pointer
			  for(int int_reg = 0; int_reg < nregions; ++int_reg)
    
			    evaluate_regional_parameters(prop_country[int_reg].det_model_params, theta.param_list, gmip, int_reg, prop_country[int_reg].population, prop_country[int_reg].total_population, base_mix, update_flags);
    
			  theta_i->param_value = tempvec; // and then re-set to the original address

			  // EVALUATE THE NEW LIKELIHOOD (ESSENTIALLY, ADJUST THE LIKELIHOOD STRUCTURE FOR THE PROPOSED VALUES)
			  fn_log_likelihood(prop_lfx, 
					    prop_country,
					    0,
					    theta_i->flag_transmission_model,
					    theta_i->flag_reporting_model,
					    theta_i->flag_GP_likelihood,
					    theta_i->flag_Hosp_likelihood,
					    theta_i->flag_Viro_likelihood,
					    theta_i->flag_Sero_likelihood,
					    gmip,
					    theta);

			  log_accep += prop_lfx.total_lfx - lfx.total_lfx;

			}
		    }
		  }
		  // dbl_A and dbl_U variables need to be defined to govern acceptance.
		  double dbl_A = (log_accep < 0.0) ? log_accep : 0.0;
		  double dbl_U = (log_accep < 0.0) ? gsl_sf_log(gsl_rng_uniform(r)) : -1.0;

		  if(dbl_U < dbl_A) // PROPOSAL IS ACCEPTED
		    {
		      
		      // UPDATE THE UPDATEABLE_MODEL_PARAMETER STRUCTURE
		      for(int int_component = 0; int_component < gsl_vector_int_get(block_size, int_param); int_component++)
			{
			  int a_component = gsl_vector_int_get(a[int_param], int_component);
			  gsl_vector_int_set(theta_i->number_accepted_moves, a_component, gsl_vector_int_get(theta_i->number_accepted_moves, a_component) + 1);
			  gsl_vector_set(theta_i->param_value, a_component, gsl_vector_get(theta_i->proposal_value, a_component));
			}
		      theta_i->log_prior_dens += theta_i->proposal_log_prior_dens;

		      // COPY ELEMENTS OF THE PROPOSAL LIKELIHOOD TO THE MODEL LIKELIHOOD STRUCTURE OR!!! UPDATE CHILD NODES PRIOR DENSITY
		      if(theta_i->flag_any_child_nodes)
			{ // NEED TO UPDATE THE PRIOR DENSITY OF THE CHILD NODES
			  for(int int_child_node = 0; int_child_node < theta.size_param_list; int_child_node++)
			    if(theta_i->flag_child_nodes[int_child_node])
			      theta.param_list[int_child_node].log_prior_dens = theta.param_list[int_child_node].proposal_log_prior_dens;
			}
		      else {
			// MEMCPY THE LIKELIHOOD STRUCTURE
			likelihood_memcpy(lfx, prop_lfx);

			for(int int_reg = 0; int_reg < nregions; int_reg++){
			  // COPY ELEMENTS OF THE PROPOSAL REGION TO THE ACCEPTED REGION
			  // SPECIFICALLY THE det_model_params MEMBER AND THE region_modstats
			  // det_model_params:
			  regional_model_params_memcpy(state_country[int_reg].det_model_params, prop_country[int_reg].det_model_params, update_flags);
			  
			  // region_modstats:
			  model_statistics_memcpy(state_country[int_reg].region_modstats, prop_country[int_reg].region_modstats,
						  theta_i->flag_transmission_model,
						  theta_i->flag_GP_likelihood && (bool) gmip.l_GP_consultation_flag,
						  theta_i->flag_Hosp_likelihood && (bool) gmip.l_Hospitalisation_flag,
						  theta_i->flag_Sero_likelihood,
						  theta_i->flag_Viro_likelihood && (bool) gmip.l_Viro_data_flag);
			  
			}
			
		      }
		      
		    }
		  else {
		    // PRETTY MUCH NEED TO TAKE THE OPPOSITE STEPS OF THE ACCEPTANCE CLAUSE
		    for(int int_component = 0; int_component < gsl_vector_int_get(block_size, int_param); int_component++)
		      {
			int a_component = gsl_vector_int_get(a[int_param], int_component);
			gsl_vector_set(theta_i->proposal_value, a_component, gsl_vector_get(theta_i->param_value, a_component));
		      }
		    for(int int_reg = 0; int_reg < nregions; int_reg++){
		      // COPY ELEMENTS OF THE PROPOSAL REGION TO THE ACCEPTED REGION
		      // SPECIFICALLY THE det_model_params MEMBER AND THE region_modstats
		      // det_model_params:
		      regional_model_params_memcpy(prop_country[int_reg].det_model_params, state_country[int_reg].det_model_params, update_flags);
			
		      // region_modstats:
		      model_statistics_memcpy(prop_country[int_reg].region_modstats, state_country[int_reg].region_modstats,
					      theta_i->flag_transmission_model,
					      theta_i->flag_GP_likelihood && (bool) gmip.l_GP_consultation_flag,
					      theta_i->flag_Hosp_likelihood && (bool) gmip.l_Hospitalisation_flag,
					      theta_i->flag_Sero_likelihood,
					      theta_i->flag_Viro_likelihood && (bool) gmip.l_Viro_data_flag);

		    }

		    // COPY ELEMENTS OF THE PROPOSAL LIKELIHOOD TO THE MODEL LIKELIHOOD STRUCTURE OR!!! UPDATE CHILD NODES PRIOR DENSITY
		    if(!theta_i->flag_any_child_nodes)
		      likelihood_memcpy(prop_lfx, lfx);

		  }

		}
	      // iterate the posterior_mean and posterior_sumsq members of the updateable_model_parameter structure regardless of whether the move is accepted or not
	      if(int_iter >= simulation_parameters.burn_in)
		{
		  for(int int_i = 0; int_i < theta_i->param_value->size; int_i++)
		    {
		      gsl_vector_set(theta_i->posterior_mean, // ADJUST THE ESTIMATE FOR THE POSTERIOR MEAN
				     int_i,
				     gsl_vector_get(theta_i->posterior_mean, int_i) +
				     (gsl_vector_get(theta_i->param_value, int_i) / (CHAIN_LENGTH)));

		      gsl_vector_set(theta_i->posterior_sumsq, // ADJUST THE POSTERIOR ESTIMATE FOR \mathbb{E} X^2
				     int_i,
				     gsl_vector_get(theta_i->posterior_sumsq, int_i) +
				     (gsl_pow_2(gsl_vector_get(theta_i->param_value, int_i)) / (CHAIN_LENGTH)));

		      // IF THE ITERATION NUMBER IS ALSO DIVIBILSE BY THE THINNING INTERVAL, THEN WE CAN SEND THE STATE OF THE CHAIN TO OUTPUT
		      if(!((int_iter + 1 - simulation_parameters.burn_in) % simulation_parameters.thin_output_every))
			output_codas[int_param].write(reinterpret_cast<char const*>(gsl_vector_ptr(theta_i->param_value, int_i)),
						      sizeof(double)); // SHOULDN'T NEED TO CHECK THAT THE PARAMETER IS UPDATED...
		    }
		}

	      // IF THE ITERATION NUMBER IS LESS THAN THE BURN_IN AND IS DIVISIBLE BY THE ADAPTIVE INTERVAL, THEN ADAPT PROPOSAL VARIANCES
	      if((int_iter < simulation_parameters.adaptive_phase) && int_iter >= 0 && !((int_iter + 1) % simulation_parameters.adapt_every))
		adapt_proposal_variance(*theta_i, simulation_parameters);

	    } // IF PARAMETER IS UPDATEABLE

	} // END FOR(int_param < size_param_list)

	  // OUTPUT TO THE LIKELIHOOD CHAIN
      if(int_iter >= simulation_parameters.burn_in && !((int_iter + 1 - simulation_parameters.burn_in) % simulation_parameters.thin_output_every))
	output_coda_lfx.write(reinterpret_cast<char const*>(&lfx.total_lfx), sizeof(double));

      // OUTPUT MODEL STATISTICS
      if(int_iter >= simulation_parameters.burn_in && !((int_iter + 1 - simulation_parameters.burn_in) % simulation_parameters.thin_stats_every))
	{
	  for(int int_reg = 0; int_reg < nregions; int_reg++)
	    {
	      model_statistics_aggregate(output_NNI[int_reg], state_country[int_reg].region_modstats, gmip.l_reporting_time_steps_per_day);
	      if(simulation_parameters.oType == cMCMC)
		{	 
		  for(int int_i = 0; int_i < output_NNI[int_reg]->size1; int_i++)
		    for(int int_j = 0; int_j < output_NNI[int_reg]->size2; int_j++)
		      {
			file_NNI[int_reg].write(reinterpret_cast<const char*>(gsl_matrix_ptr(output_NNI[int_reg], int_i, int_j)), sizeof(double));
			file_Sero[int_reg].write(reinterpret_cast<const char*>(gsl_matrix_ptr(state_country[int_reg].region_modstats.d_seropositivity, int_i, int_j)), sizeof(double));
			if(gmip.l_GP_consultation_flag)
			  {
			    file_GP[int_reg].write(reinterpret_cast<const char*>(gsl_matrix_ptr(state_country[int_reg].region_modstats.d_Reported_GP_Consultations, int_i, int_j)), sizeof(double));
			    file_Viro[int_reg].write(reinterpret_cast<const char*>(gsl_matrix_ptr(state_country[int_reg].region_modstats.d_viropositivity, int_i, int_j)), sizeof(double));
			  }
			if(gmip.l_Hospitalisation_flag)
			  file_Hosp[int_reg].write(reinterpret_cast<const char*>(gsl_matrix_ptr(state_country[int_reg].region_modstats.d_Reported_Hospitalisations, int_i, int_j)), sizeof(double));
		      }
		}
	      else if(simulation_parameters.oType == cSMC)
		{
		  for(int int_i = 0; int_i < state_country[int_reg].region_modstats.d_NNI->size1; int_i++)
		    for(int int_j = 0; int_j < state_country[int_reg].region_modstats.d_NNI->size2; int_j++)		  
		      file_NNI[int_reg].write(reinterpret_cast<const char*>(gsl_matrix_ptr(state_country[int_reg].region_modstats.d_NNI, int_i, int_j)), sizeof(double));
		  state_country[int_reg].region_modstats.d_end_state->write(file_state[int_reg]);
		}
	    }
	}

      // UPDATE LIKELIHOOD POSTERIOR STATISTICS..
      if(int_iter >= simulation_parameters.burn_in)
	{
	  lfx.bar_lfx += lfx.total_lfx / ((double) CHAIN_LENGTH);
	  lfx.sumsq_lfx += gsl_pow_2(lfx.total_lfx) / ((double) CHAIN_LENGTH);
	  prop_lfx.bar_lfx = lfx.bar_lfx;
	  prop_lfx.sumsq_lfx = lfx.sumsq_lfx;
	}

      // OUTPUT MCMC SAMPLER PROGRESS REPORTS...
	  if (int_progress_report < simulation_parameters.num_progress_reports) {
		  if(int_iter + 1 == gsl_vector_int_get(adaptive_progress_report_iterations, int_progress_report))
		write_progress_report("adaptive_report", ++int_progress_report, int_iter + 1, CHAIN_LENGTH,
					  theta, lfx, false, true, true);
		  else if(int_iter + 1 == gsl_vector_int_get(chain_progress_report_iterations, int_progress_report))
		write_progress_report("posterior_report", ++int_progress_report, int_iter + 1 - simulation_parameters.burn_in, CHAIN_LENGTH,
					  theta, lfx, true, true, false);
	  }
	  if(int_iter + 1 == simulation_parameters.adaptive_phase) int_progress_report = 0;

      // RESET COUNTERS WHERE NECESSARY - if start of a new adaptive phase
      // or the end of the burn-in
      if(((!((int_iter + 1) % simulation_parameters.adapt_every)) && (int_iter < simulation_parameters.adaptive_phase)) || ((int_iter + 1) == simulation_parameters.burn_in))
	reset_counters(theta);
      
    } // END FOR(int_iter < num_iterations)


  likelihood_free(prop_lfx);
  for(int int_i = 0; int_i < nregions; int_i++){
    gsl_matrix_free(output_NNI[int_i]);
    // TERMINATE STATISTIC OUTPUT FILES

    Region_free(prop_country[int_i], gmip);
  }
  delete [] prop_country;
  for(int_param = 0; int_param < theta.size_param_list; int_param++)
    {
      if(theta.param_list[int_param].flag_update)
	{
	  int iTemp = (int) ((CHAIN_LENGTH) / simulation_parameters.thin_output_every);
	  output_codas[int_param].write(reinterpret_cast<char const*>(&theta.param_list[int_param].param_value->size), sizeof(int));
	  output_codas[int_param].write(reinterpret_cast<char const*>(&iTemp), sizeof(int));   // TERMINATE PARAMETER OUTPUT FILES - by adding dimensions.
	  output_codas[int_param].close();
	}
      gsl_vector_int_free(a[int_param]);
      gsl_vector_int_free(b[int_param]);
    }
  delete [] output_NNI;
  fstream_model_statistics_close(file_NNI, nregions, NUM_AGE_GROUPS, gmip.l_duration_of_runs_in_days * (simulation_parameters.oType == cMCMC ? 1 : gmip.l_reporting_time_steps_per_day),
				 (CHAIN_LENGTH) / simulation_parameters.thin_stats_every);
  if(simulation_parameters.oType == cMCMC)
    {
      fstream_model_statistics_close(file_Sero, nregions, NUM_AGE_GROUPS, gmip.l_duration_of_runs_in_days, (CHAIN_LENGTH) / simulation_parameters.thin_stats_every);
      if(gmip.l_GP_consultation_flag){
	fstream_model_statistics_close(file_GP, nregions, NUM_AGE_GROUPS, gmip.l_duration_of_runs_in_days, (CHAIN_LENGTH) / simulation_parameters.thin_stats_every);
	fstream_model_statistics_close(file_Viro, nregions, NUM_AGE_GROUPS, gmip.l_duration_of_runs_in_days, (CHAIN_LENGTH) / simulation_parameters.thin_stats_every);
      }
      if(gmip.l_Hospitalisation_flag)
	fstream_model_statistics_close(file_Hosp, nregions, NUM_AGE_GROUPS, gmip.l_duration_of_runs_in_days, (CHAIN_LENGTH) / simulation_parameters.thin_stats_every);
    }
  if(simulation_parameters.oType == cSMC)
    delete [] file_state;

  delete [] a;
  delete [] b;
  delete [] output_codas;
  gsl_vector_int_free(num_component_updates);
  gsl_vector_int_free(block_size);
  gsl_vector_int_free(adaptive_progress_report_iterations);
  gsl_vector_int_free(chain_progress_report_iterations);
}



