#include "RTM_StructDefs.h"
#include "RTM_FunctDefs.h"
#include "RTM_updParams.h"
#include "gsl_vec_ext.h"
#include <cfloat>

using namespace std;
using std::string;

// FUNCTION PROTOTYPE
void fn_reporting_model(gsl_matrix*, const gsl_matrix*, const gsl_matrix*, const gsl_matrix*,
			global_model_instance_parameters, int, const gsl_vector*,
			gsl_vector*, const int int_num_threads = 1);

// OUTPUT_PER_SELECTED_PERIOD TAKES THE OUTPUTS OF THE MODEL AND CONVERTS THEM TO OUTPUTS COMPARABLE WITH THE INPUT DATA
void output_per_selected_period(
				const int &time_steps_per_day, // HOW MANY TIME STEPS DO YOU WISH TO CONTRACT OVER
				const gsl_matrix* matrix_uncontracted, // THE INPUT: PER TIME STEPS PER DAY.
				gsl_matrix* matrix_observed) // THE MATRIX OF RESULTS BUT WITH OUTPUTS PER TIME STEPS, AS THESE COMPUTED IN THE MODEL.
{// WE NEED TO KNOW ITS DIMENSIONS BEFORE WE PUT IN! THE 1-DIMENSION IS [maxiter/(time_steps_per_day*output_days)]
  gsl_matrix_set_zero(matrix_observed); // INITIALIZATION

  for (int i = 0; i < (matrix_observed->size1); ++i)
    {
      for (int j = 0; j < matrix_observed->size2; ++j)
	{
	  for (int k = 0; k <= (time_steps_per_day - 1); ++k) // MAKE SURE THAT YOU ARE PASSING A ZERO MATRIX INSIDE HERE
	    {
	      gsl_matrix_set(matrix_observed, i, j, gsl_matrix_get(matrix_observed, i, j) + gsl_matrix_get(matrix_uncontracted, (time_steps_per_day * i) + k, j)); // YOU NEED TO GO UP TO -1
	    }// FOR
	}// FOR
    }//FOR
}



// SEASONAL R0_FUNCTION GIVES THE VECTOR OF THE SEASONAL R0 FOR THE LENGTH OF THE RUN OF THE MODEL BASED ON THE R0 AT THE START OF THE SIMULATION (R0INIT) AND THE AMPLITUDE OF OSCILLATION OF R0
void Seasonal_R0_function(
			  gsl_vector* Seasonal_R0, // memory is assumed to have already been allocated for this vector
			  const double R0init,
			  const double amplitude,
			  const gsl_vector* phase_shift)
{

  gsl_vector_memcpy(Seasonal_R0, phase_shift); // R0(t) = phase_difference(t) (this line should also provide any necessary error checking)
  gsl_vector_scale(Seasonal_R0, amplitude); // R0(t) = A * phase_difference(t)
  gsl_vector_add_constant(Seasonal_R0, R0init); // R0(t) = R0(0) + (A * phase_difference(t))

}

void initialising_Deterministic_SE1E2I1I2R_AG_RF_CIN2(regional_model_params in_dmp,
						      global_model_instance_parameters global_params,
						      const gsl_vector* evector_MIXMAT_normalised,
						      gsl_vector* l_E1_0, gsl_vector* l_E2_0, gsl_vector* l_I1_0, gsl_vector* l_I2_0)

{
  // REDEFINE VARIABLES
  double timestepsperday = double(global_params.l_transmission_time_steps_per_day);

  gsl_vector* rho = gsl_vector_alloc(in_dmp.l_average_infectious_period->size2);
  gsl_vector* sigma = gsl_vector_alloc(in_dmp.l_latent_period->size2);

  for(int int_i = 0; int_i < rho->size; int_i++)
    {
      gsl_vector_set(rho, int_i, 2 / (gsl_matrix_get(in_dmp.l_average_infectious_period, 0, int_i) * timestepsperday));
      gsl_vector_set(sigma, int_i, 2 / (gsl_matrix_get(in_dmp.l_latent_period, 0, int_i) * timestepsperday));
    }
  double alpha = gsl_sf_exp(in_dmp.l_EGR/timestepsperday)-1; // ALWAYS GREATER THAN 1

  // SPECIFY I0 & IMPORTATION RATE AS IN THE OBJECTIVE FUNCTION
  gsl_vector* I0_vec = gsl_vector_alloc(NUM_AGE_GROUPS);
 
  for (int i = 0; i < NUM_AGE_GROUPS; ++i)
    {
      gsl_vector_set(I0_vec, i, gsl_vector_get(evector_MIXMAT_normalised, i) * in_dmp.l_I0);
      gsl_vector_set(l_I1_0, i, gsl_vector_get(I0_vec, i) / (1 + ((gsl_matrix_get(in_dmp.l_relative_infectiousness_I2_wrt_I1, 0, i) * gsl_vector_get(rho, i)) / (alpha + gsl_vector_get(rho, i)))));
      gsl_vector_set(l_I2_0, i, gsl_vector_get(l_I1_0, i) * gsl_vector_get(rho, i)/(alpha + gsl_vector_get(rho, i)));
      gsl_vector_set(l_E2_0, i, gsl_vector_get(l_I1_0, i) * (alpha + gsl_vector_get(rho, i))/gsl_vector_get(sigma, i));
      gsl_vector_set(l_E1_0, i, gsl_vector_get(l_E2_0, i) * (alpha + gsl_vector_get(sigma, i))/gsl_vector_get(sigma, i));
    }

  gsl_vector_free(I0_vec);
  gsl_vector_free(rho);
  gsl_vector_free(sigma);
}

// Collapsing age- and vaccination state-specific vectors into an age only-dependent quantity.
void collapse_vac(gsl_vector* vec_age, const gsl_vector* vec_vac_age){
  gsl_vector_set_zero(vec_age);
  for(int intA = 0; intA < NUM_AGE_GROUPS; ++intA)
    for(int intV = 0; intV <= MAX_VAX_DOSES; ++intV)
      gsl_vector_set(vec_age, intA, gsl_vector_get(vec_age, intA) + gsl_vector_get(vec_vac_age, STATE_IDX(intV, intA)));
}

void prob_infection_RF_MA(
			  gsl_vector* prob_infection, // This is an output
			  gsl_matrix* force_infectious_contact_matrix,
			  const gsl_vector* epsilon,
			  const gsl_vector* importations,
			  const gsl_vector* I_1,
			  const gsl_vector* I_2,
			  const transmission_kernel& tk,
			  const double delta_t)
{

  
  gsl_vector* Itemp = gsl_vector_calloc(NUM_AGE_GROUPS);
  gsl_vector* Iages = gsl_vector_calloc(NUM_AGE_GROUPS);

  collapse_vac(Itemp, I_1);
  collapse_vac(Iages, I_2);
  
  // gsl_vector_memcpy(I, I_2);
  gsl_vector_mul(Iages, epsilon);
  gsl_vector_add(Iages, Itemp);
  
  if(tk == cREEDFROST){
    for (int a = 0; a < prob_infection->size; ++a)
      {
	// gsl_vector* p_beta = gsl_vector_alloc(prob_infection->size);
	double lprod = 0;
	for (int b = 0; b < prob_infection->size; ++b)
	  {
	    lprod += gsl_vector_get(Iages, b) * gsl_sf_log(1 - gsl_matrix_get(force_infectious_contact_matrix, a, b));
	    // gsl_vector_set(p_beta,
	    // 		   b,
	    // 		   pow(1 - gsl_matrix_get(force_infectious_contact_matrix, a, b), gsl_vector_get(I, b)));
	  }
	// gsl_vector_set(prob_infection, a, 1 - gsl_double_product_of_vector_elements(p_beta));
	// gsl_vector_free(p_beta);
	gsl_vector_set(prob_infection, a, 1 - gsl_sf_exp(lprod));
      }
  } else if(tk == cMASSACTION){
    
    gsl_blas_dsymv(CblasLower, 1.0, force_infectious_contact_matrix, Iages, 0.0, prob_infection);

  }
  gsl_vector_free(Iages);
  gsl_vector_free(Itemp);
  
  // Allow for outside importations and scale according to the time-step size.
  gsl_vector_add(prob_infection, importations);
  gsl_vector_scale(prob_infection, delta_t);

}

void Deterministic_S_E1_E2_I1_I2_R_AG_RF(					 // THE MODEL MODIFIES ALL THE PARAMETERS THAT ARE PASSED BY REFERENCE AND ARE NOT CONSTANT
					 gsl_matrix* l_S,
					 gsl_matrix* l_E_1, 
					 gsl_matrix* l_E_2, 
					 gsl_matrix* l_I_1, 
					 gsl_matrix* l_I_2, 
					 gsl_matrix* l_R_pos, 
					 gsl_matrix* l_R_neg,
					 gsl_matrix* l_W,
					 gsl_matrix* l_WS,
					 gsl_matrix* l_NNI,
					 gsl_matrix* l_Delta_Dis,
					 gsl_matrix* l_Seropositivity,
					 gsl_matrix* l_internal_AR,
					 gsl_matrix* l_Prevalence,
					 model_state* l_end_state,
					 const regional_model_params& in_dmp,
					 const gsl_vector* l_R0_t,
					 const gsl_vector* regional_population_by_age,
					 const gsl_vector* E1_0,
					 const gsl_vector* E2_0,
					 const gsl_vector* I1_0,
					 const gsl_vector* I2_0,
					 const global_model_instance_parameters& gmip,
					 const mixing_model& l_MIXMOD_ADJUSTED2,
					 const rtmData* Vaccination,
					 const rtmData* VBooster,
					 const rtmData* VFourth)
{
  // FUNCTION DEFINITION 
  int timestepsperday = gmip.l_transmission_time_steps_per_day;
  int time_points = l_S->size1;
  // ARRAYS USED TO STORE THE CALCULATIONS OF SUMMARY STATISTICS
  gsl_matrix* Number_New_Infected = gsl_matrix_calloc(time_points, NUM_AGE_GROUPS);
  gsl_matrix* Delta_Disease = (Vaccination) ? gsl_matrix_calloc(time_points, NUM_AGE_GROUPS) : Number_New_Infected;
  gsl_matrix* p_lambda = gsl_matrix_calloc(time_points, NUM_AGE_GROUPS);
  // ARRAY TO HOLD SELECTED CONTACT MATRICES AT TIME t.
  gsl_matrix* P = gsl_matrix_alloc(NUM_AGE_GROUPS, NUM_AGE_GROUPS);
  
  int mix_interval;

  // INITIALISING: SEIR view variables
  // CHECK EVERYTHING HAS BEEN INITIALISED TO ZERO!
  gsl_vector_view S_view = gsl_matrix_row(l_S, 0);
  gsl_vector_view I1_view = gsl_matrix_row(l_I_1, 0);
  gsl_vector_view I2_view = gsl_matrix_row(l_I_2, 0);
  gsl_vector_view E1_view = gsl_matrix_row(l_E_1, 0);
  gsl_vector_view E2_view = gsl_matrix_row(l_E_2, 0);
  gsl_vector_view R_pos_view = gsl_matrix_row(l_R_pos, 0);
  // gsl_vector_set_all(&R_pos_view.vector, 0);
  gsl_vector_view R_neg_view = gsl_matrix_row(l_R_neg, 0);
  // gsl_vector_set_all(&R_neg_view.vector, 0);
  gsl_vector_view W_view = gsl_matrix_row(l_W, 0);
  // gsl_vector_set_all(&W_pos_view.vector, 0);
  gsl_vector_view WS_view = gsl_matrix_row(l_WS, 0);

  // INITIALISING: Other view variables needed
  gsl_vector_view P_view;
  gsl_vector_view dA_view;

  // RE-SETTING COUNTERS OF POSITIVITY (MEASURED AS THE FRACTION CURRENTLY UNSUSCEPTIBLE) AND INFECTION ATTACK RATE
  gsl_matrix_set_all(l_Seropositivity, 1.0);
  gsl_matrix_set_all(l_internal_AR, 1.0);
  gsl_matrix_set_all(l_Prevalence, 0.0);

  // ASSUME ONLY UNVACCINATED STATES ARE INITIALLY NON-ZERO  

  // SELECT THE CORRECT (time-0) SCALED MIXING MATRIX
  gsl_matrix_memcpy(P, l_MIXMOD_ADJUSTED2.MIXMAT_scaled[0]);

  // GET SOME UNVACCINATED VIEWS OF THE FIRST ROW OF THE STATE MATRICES
  // AND FILL THOSE ROWS FOR THE E1, E2, I1, I2 STATES // HERE!! 20220131
  gsl_vector_view S_V0_view = gsl_vector_subvector(&S_view.vector, 0, NUM_AGE_GROUPS);
  gsl_vector_view E1_V0_view = gsl_vector_subvector(&E1_view.vector, 0, NUM_AGE_GROUPS);
  gsl_vector_memcpy(&E1_V0_view.vector, E1_0);
  gsl_vector_view E2_V0_view = gsl_vector_subvector(&E2_view.vector, 0, NUM_AGE_GROUPS);
  gsl_vector_memcpy(&E2_V0_view.vector, E2_0);
  gsl_vector_view I1_V0_view = gsl_vector_subvector(&I1_view.vector, 0, NUM_AGE_GROUPS);
  gsl_vector_memcpy(&I1_V0_view.vector, I1_0);
  gsl_vector_view I2_V0_view = gsl_vector_subvector(&I2_view.vector, 0, NUM_AGE_GROUPS);
  gsl_vector_memcpy(&I2_V0_view.vector, I2_0);
  for (int a = 0; a < NUM_AGE_GROUPS; ++a)
    {
      gsl_matrix_set(l_S, 0, a, gsl_vector_get(regional_population_by_age, a) * gsl_vector_get(in_dmp.l_init_prop_sus, a));
      gsl_matrix_set(l_S, 0, a, gsl_matrix_get(l_S, 0, a) -
		     gsl_matrix_get(l_I_1, 0, a) -
		     gsl_matrix_get(l_E_1, 0, a) -
		     gsl_matrix_get(l_I_2, 0, a) -
		     gsl_matrix_get(l_E_2, 0, a));
      gsl_matrix_set(l_R_neg, 0, a, gsl_vector_get(regional_population_by_age, a) * (1 - gsl_vector_get(in_dmp.l_init_prop_sus, a))); // OLD CODE, l_R WAS IGNORED, SO HOPEFULLY CAN RE-PURPOSE IT FOR THE NEW SWAB POSITIVE STATE
      P_view = gsl_matrix_row(P, a);
      dA_view = gsl_matrix_row(in_dmp.l_average_infectious_period, 0);
      gsl_vector_div(&P_view.vector, &dA_view.vector);
    }
  gsl_matrix_scale(P, gsl_vector_get(l_R0_t, 0));

  // calculate the probability of infection
  gsl_vector_view p_lambda_view = gsl_matrix_row(p_lambda, 0);
  gsl_vector_view epsilon_view = gsl_matrix_row(in_dmp.l_relative_infectiousness_I2_wrt_I1, 0);
  gsl_vector_view importations_view = gsl_matrix_row(in_dmp.l_importation_rate, 0);
  gsl_vector_view Number_New_Infected_view = gsl_matrix_row(Number_New_Infected, 0);

  prob_infection_RF_MA(&p_lambda_view.vector,
		       P,
		       &epsilon_view.vector,
		       &importations_view.vector,
		       &I1_view.vector,
		       &I2_view.vector,
		       gmip.l_tk,
		       1 / ((double) timestepsperday));
  
  // Get number of new infections.. probability of infection multiplied by number of susceptibles.
  gsl_matrix_set_row(Number_New_Infected, 0, &S_V0_view.vector);
  gsl_vector_mul(&Number_New_Infected_view.vector, &p_lambda_view.vector);
  gsl_matrix_set_row(Delta_Disease, 0, &Number_New_Infected_view.vector);
  
  // NEED TO AVOID THIS LOOP IF S(t) has any negative values

  // UPDATING THE MODEL:
  int nday = 0;
  
  for (int t = 0; t < time_points - 1; ++t)// YOU PUT -1 BECAUSE YOUR VECTORS HAVE LENGTHS UP TO TIME_POINTS AND WE ARE WORKING WITH t+1
    { 

      // PICK OUT MIXING MATRIX FOR THIS TIME-STEP
      mix_interval = (l_MIXMOD_ADJUSTED2.num_breakpoints == 0) ? 0 : mix_timecat(t + 1, l_MIXMOD_ADJUSTED2.breakpoints, timestepsperday);
      gsl_matrix_memcpy(P, l_MIXMOD_ADJUSTED2.MIXMAT_scaled[mix_interval]);

      for (int a = 0; a < NUM_AGE_GROUPS; ++a)
	{
	  for(int v = 0; v <= MAX_VAX_DOSES; ++v)
	    {
	      double pi = (v == 0) ? 0 :
		((v == 1) ? gsl_matrix_get(in_dmp.l_vacc1_infect, t, a) :
		 ((v == 2) ? gsl_matrix_get(in_dmp.l_vaccn_infect, t, a) :
		  ((v == 3) ? gsl_matrix_get(in_dmp.l_vaccb_infect, t, a) : gsl_matrix_get(in_dmp.l_vacc4_infect, t, a)
		   )
		  )
		 );
	      double alpha = (v == 0) ? 0 :
		((v == 1) ? gsl_matrix_get(in_dmp.l_vacc1_disease, t, a) :
		 ((v == 2) ? gsl_matrix_get(in_dmp.l_vaccn_disease, t, a) :
		  ((v == 3) ? gsl_matrix_get(in_dmp.l_vaccb_disease, t, a) : gsl_matrix_get(in_dmp.l_vacc4_disease, t, a)
		   )
		  )
		 );
	      double vacc_in = ((v == 0) ? 0 :
				((v == 1) ? Vaccination->getCount(nday, a) :
				 ((v == 2) ? Vaccination->getDenom(nday, a) :
				  ((v == 3) ? VBooster->getCount(nday, a) : VFourth->getCount(nday, a)
				   )
				  )
				 )) / timestepsperday;
	      double vacc_out = ((v == 0) ? Vaccination->getCount(nday, a) :
				 ((v == 1) ? Vaccination->getDenom(nday, a) :
				  ((v == 2) ? VBooster->getCount(nday, a) :
				   ((v == 3) ? VFourth->getCount(nday, a) : 0
				    )
				   )
				  )) / timestepsperday;

	      // NUMBER OF NEW INFECTIONS FROM CURRENT SUSCEPTIBLE GROUP DEFINED BY (v,a)
	      double adj_S = (1 - pi) * (1 - vacc_out) *
		(gsl_matrix_get(l_S, t, STATE_IDX(v, a)) + gsl_matrix_get(l_WS, t, STATE_IDX(v, a))) * gsl_matrix_get(p_lambda, t, a);

	      gsl_matrix_set(Number_New_Infected, t + 1, a, gsl_matrix_get(Number_New_Infected, t + 1, a) + adj_S);
	      gsl_matrix_set(Delta_Disease, t + 1, a, gsl_matrix_get(Delta_Disease, t + 1, a) + (1 - alpha) * adj_S);

	      gsl_matrix_set(l_S, t + 1, STATE_IDX(v, a),
			     (gsl_matrix_get(l_S, t, STATE_IDX(FN_MAX(v - 1, 0), a)) * vacc_in) +
			     gsl_matrix_get(l_S, t, STATE_IDX(v, a)) * (1 - vacc_out) * (1 - ((1 - pi) * gsl_matrix_get(p_lambda, t, a)))
			     );

	      gsl_matrix_set(l_E_1, t + 1, STATE_IDX(v, a),
			     (gsl_matrix_get(l_E_1, t, STATE_IDX(FN_MAX(v - 1, 0), a)) * vacc_in) +
			     (gsl_matrix_get(l_E_1, t, STATE_IDX(v, a)) * (1 - vacc_out) * (1 - (2 / (gsl_matrix_get(in_dmp.l_latent_period, t, a) * timestepsperday)))) +
			     adj_S);

	      gsl_matrix_set(l_E_2, t + 1, STATE_IDX(v, a),
			     (gsl_matrix_get(l_E_2, t, STATE_IDX(FN_MAX(v - 1, 0), a)) * vacc_in) +
			     (gsl_matrix_get(l_E_2, t, STATE_IDX(v, a)) * (1 - vacc_out) * (1 - (2 / (gsl_matrix_get(in_dmp.l_latent_period, t, a) * timestepsperday)))) +
			     (gsl_matrix_get(l_E_1, t, STATE_IDX(v, a)) * (1 - vacc_out) * 2 / (gsl_matrix_get(in_dmp.l_latent_period, t, a) * timestepsperday)));

	      gsl_matrix_set(l_I_1, t + 1, STATE_IDX(v, a),
			     (gsl_matrix_get(l_I_1, t, STATE_IDX(FN_MAX(v - 1, 0), a)) * vacc_in) +
			     (gsl_matrix_get(l_I_1, t, STATE_IDX(v, a)) * (1 - vacc_out) * (1 - (2 / (gsl_matrix_get(in_dmp.l_average_infectious_period, t, a) * timestepsperday)))) +
			     (gsl_matrix_get(l_E_2, t, STATE_IDX(v, a)) * (1 - vacc_out) * 2 / (gsl_matrix_get(in_dmp.l_latent_period, t, a) * timestepsperday)));

	      gsl_matrix_set(l_I_2, t + 1, STATE_IDX(v, a),
			     (gsl_matrix_get(l_I_2, t, STATE_IDX(FN_MAX(v - 1, 0), a)) * vacc_in) +
			     (gsl_matrix_get(l_I_2, t, STATE_IDX(v, a)) * (1 - vacc_out) * (1 - (2 / (gsl_matrix_get(in_dmp.l_average_infectious_period, t, a) * timestepsperday)))) +
			     (gsl_matrix_get(l_I_1, t, STATE_IDX(v, a)) * (1 - vacc_out) * 2 / (gsl_matrix_get(in_dmp.l_average_infectious_period, t, a) * timestepsperday)));
	      
	      gsl_matrix_set(l_R_pos, t + 1, STATE_IDX(v, a),
			     (gsl_matrix_get(l_R_pos, t, STATE_IDX(FN_MAX(v - 1, 0), a)) * vacc_in) +
			     (gsl_matrix_get(l_R_pos, t, STATE_IDX(v, a)) * (1 - vacc_out) * (1 - (1 / (gsl_matrix_get(in_dmp.l_r1_period, t, a) * timestepsperday)))) +
			     (gsl_matrix_get(l_I_2, t, STATE_IDX(v, a)) * (1 - vacc_out) * 2 / (gsl_matrix_get(in_dmp.l_average_infectious_period, t, a) * timestepsperday)));

	      gsl_matrix_set(l_R_neg, t + 1, STATE_IDX(v, a),
			     (gsl_matrix_get(l_R_neg, t, STATE_IDX(FN_MAX(v - 1, 0), a)) * vacc_in) +
			     (gsl_matrix_get(l_R_neg, t, STATE_IDX(v, a)) * (1 - vacc_out) * (1 - (2 / (gsl_matrix_get(in_dmp.l_waning_period, t, a) * timestepsperday)))) +
			     (gsl_matrix_get(l_R_pos, t, STATE_IDX(v, a)) * (1 - vacc_out) / (gsl_matrix_get(in_dmp.l_r1_period, t, a) * timestepsperday)));

	      gsl_matrix_set(l_W, t + 1, STATE_IDX(v, a),
			     (gsl_matrix_get(l_W, t, STATE_IDX(FN_MAX(v - 1, 0), a)) * vacc_in) +
			     (gsl_matrix_get(l_W, t, STATE_IDX(v, a)) * (1 - vacc_out) * (1 - (2 / (gsl_matrix_get(in_dmp.l_waning_period, t, a) * timestepsperday)))) +
			     (gsl_matrix_get(l_R_neg, t, STATE_IDX(v, a)) * (1 - vacc_out) * 2 / (gsl_matrix_get(in_dmp.l_waning_period, t, a) * timestepsperday)));

	      gsl_matrix_set(l_WS, t + 1, STATE_IDX(v, a),
			     (gsl_matrix_get(l_WS, t, STATE_IDX(FN_MAX(v - 1, 0), a)) * vacc_in) +
			     (gsl_matrix_get(l_WS, t, STATE_IDX(v, a)) * (1 - vacc_out) * (1 - ((1 - pi) * gsl_matrix_get(p_lambda, t, a)))) +
			     (gsl_matrix_get(l_W, t, STATE_IDX(v, a)) * (1 - vacc_out) * 2 / (gsl_matrix_get(in_dmp.l_waning_period, t, a) * timestepsperday)));

	      if((t + 1) % timestepsperday){ // THESE ARE OUTPUT MATRICES AND ARE NOT CALCULATED EVERY (DELTA T) DAYS.. THESE MATRICES ARE DESIGNED FOR DAILY VALUES
		int tindex = t / timestepsperday;
		// Serology
		gsl_matrix_set(l_internal_AR, tindex, a, gsl_matrix_get(l_internal_AR, tindex, a) - (gsl_matrix_get(l_S, t + 1, STATE_IDX(v, a)) / gsl_vector_get(regional_population_by_age, a)));
		gsl_matrix_set(l_Seropositivity, tindex, a, gsl_matrix_get(l_Seropositivity, tindex, a) - ((1 - pi) * (gsl_matrix_get(l_S, t + 1, STATE_IDX(v, a)) + gsl_matrix_get(l_WS, t + 1, STATE_IDX(v, a))) / gsl_vector_get(regional_population_by_age, a)));
		// Virological positivity
		gsl_matrix_set(l_Prevalence, tindex, a, gsl_matrix_get(l_Prevalence, tindex, a) + gsl_matrix_get(l_I_1, t + 1, STATE_IDX(v, a)) + gsl_matrix_get(l_I_2, t + 1, STATE_IDX(v, a)) + gsl_matrix_get(l_R_pos, t + 1, STATE_IDX(v, a)));

	      }
	      
	    } // FOR v

	  P_view = gsl_matrix_row(P, a);
	  dA_view = gsl_matrix_row(in_dmp.l_average_infectious_period, 0);  // Time-0 for the infectious period to ensure standardisation
	  gsl_vector_div(&P_view.vector, &dA_view.vector);
	   
	} // FOR a
    
      // UPDATE THE PROBABILITY OF INFECTION

      // scaled mixing matrix
      gsl_matrix_scale(P, gsl_vector_get(l_R0_t, t + 1));

      // update vector views
      p_lambda_view = gsl_matrix_row(p_lambda, t + 1);
      epsilon_view = gsl_matrix_row(in_dmp.l_relative_infectiousness_I2_wrt_I1, t + 1);
      importations_view = gsl_matrix_row(in_dmp.l_importation_rate, t + 1);

      I1_view = gsl_matrix_row(l_I_1, t + 1);
      I2_view = gsl_matrix_row(l_I_2, t + 1);
      
      // update probability of infection
      prob_infection_RF_MA(&p_lambda_view.vector,
			   P,
			   &epsilon_view.vector,
			   &importations_view.vector,
			   &I1_view.vector,
			   &I2_view.vector,
			   gmip.l_tk,
			   1 / ((double) timestepsperday));

      // update the day number
      nday = nday + ((t % timestepsperday) == 0 ? 1 : 0);
    } // FOR t

  // // DEBUGGING //
  // // ofstream output_debug("debug.txt", ios::out | ios::trunc | ios::binary);
  // gsl_matrix* tempmat = gsl_matrix_calloc(l_S->size1, l_S->size2);
  // gsl_vector* tempvec = gsl_vector_calloc(l_S->size1);
  // gsl_matrix_add(tempmat, l_S);
  // gsl_matrix_add(tempmat, l_E_1);
  // gsl_matrix_add(tempmat, l_E_2);
  // gsl_matrix_add(tempmat, l_I_1);
  // gsl_matrix_add(tempmat, l_I_2);
  // gsl_matrix_add(tempmat, l_R_pos);
  // gsl_matrix_add(tempmat, l_R_neg);
  // gsl_matrix_add(tempmat, l_W);
  // gsl_matrix_add(tempmat, l_WS);
  // for(int inti = 0; inti < tempmat->size1; inti++){
  //   for(int inta = 0; inta < tempmat->size2; inta++)
  //     gsl_vector_set(tempvec, inti, gsl_vector_get(tempvec, inti) + gsl_matrix_get(tempmat, inti, inta));
  //   std::cout << "t = " << inti << "; people in model = " << gsl_vector_get(tempvec, inti) << std::endl;
  // }
  // gsl_matrix_free(tempmat);
  // gsl_vector_free(tempvec);
  // // END DEBUGGING //
  
  S_view = gsl_matrix_row(l_S, time_points - 1);
  E1_view = gsl_matrix_row(l_E_1, time_points - 1);
  E2_view = gsl_matrix_row(l_E_2, time_points - 1);      
  // I1_view = gsl_matrix_row(l_I_1, time_points - 1); // should already be set
  // I2_view = gsl_matrix_row(l_I_2, time_points - 1); // should already be set
  R_pos_view = gsl_matrix_row(l_R_pos, time_points - 1);
  R_neg_view = gsl_matrix_row(l_R_neg, time_points - 1);
  W_view = gsl_matrix_row(l_W, time_points - 1);
  WS_view = gsl_matrix_row(l_WS, time_points - 1);

  l_end_state->fill(&S_view.vector,
		    &E1_view.vector,
		    &E2_view.vector,
		    &I1_view.vector,
		    &I2_view.vector,
		    &R_pos_view.vector,
		    &R_neg_view.vector,
		    &W_view.vector,
		    &WS_view.vector,
		    &p_lambda_view.vector);

  // AGGREGATE THE NUMBER_NEW_INFECTED MATRIX TO THE TIME STEPS USED BY THE REPORTING MODEL
  output_per_selected_period(timestepsperday / gmip.l_reporting_time_steps_per_day, Number_New_Infected, l_NNI);
  output_per_selected_period(timestepsperday / gmip.l_reporting_time_steps_per_day, Delta_Disease, l_Delta_Dis);
  
  double dbl_isnan_check = gsl_matrix_max(l_NNI);
  if(std::isnan(dbl_isnan_check)){
    gsl_matrix_set_all(l_NNI, GSL_NEGINF);
    gsl_matrix_set_all(l_Seropositivity, -1.0);
    gsl_matrix_set_all(l_internal_AR, -1.0);
  }

  gsl_matrix_free(Delta_Disease);
  gsl_matrix_free(Number_New_Infected);
  gsl_matrix_free(p_lambda);
  gsl_matrix_free(P);

}

void propagate_SEEIIR(regional_model_params in_dmp, const gsl_vector* regional_population_by_age,
		      const global_model_instance_parameters& global_params,
		      gsl_vector* R0_t, const gsl_vector* E1_0, const gsl_vector* E2_0,
		      const gsl_vector* I1_0, const gsl_vector* I2_0,
		      gsl_matrix* d_NNI, gsl_matrix* d_Delta_Dis, gsl_matrix* d_Seropositivity, gsl_matrix* d_internal_AR, gsl_matrix* d_Prevalence,
		      model_state* d_end_state,
		      const rtmData* vaccination,
		      const rtmData* vbooster,
		      const rtmData* vfourth)
{
  int num_days = global_params.l_duration_of_runs_in_days, step_size = global_params.l_transmission_time_steps_per_day;

  const mixing_model l_MIXMOD_ADJUSTED2 = in_dmp.l_MIXMOD;

  int narr_rows = num_days * step_size;
  int narr_cols = (MAX_VAX_DOSES + 1) * NUM_AGE_GROUPS;
  gsl_matrix* S = gsl_matrix_calloc(narr_rows, narr_cols); 
  gsl_matrix* E_1 = gsl_matrix_calloc(narr_rows, narr_cols); 
  gsl_matrix* E_2 = gsl_matrix_calloc(narr_rows, narr_cols); 
  gsl_matrix* I_1 = gsl_matrix_calloc(narr_rows, narr_cols); 
  gsl_matrix* I_2 = gsl_matrix_calloc(narr_rows, narr_cols); 
  gsl_matrix* R_pos = gsl_matrix_calloc(narr_rows, narr_cols);
  gsl_matrix* R_neg = gsl_matrix_calloc(narr_rows, narr_cols);
  gsl_matrix* W = gsl_matrix_calloc(narr_rows, narr_cols);
  gsl_matrix* WS = gsl_matrix_calloc(narr_rows, narr_cols);
  
  Deterministic_S_E1_E2_I1_I2_R_AG_RF( // OUTPUTS THE NUMBER OF NEW INFECTEDS, THE INPUT TO THE REPORTING MODEL
				      S,
				      E_1, 
				      E_2, 
				      I_1, 
				      I_2, 
				      R_pos,
				      R_neg,
				      W,
				      WS,
				      d_NNI,
				      d_Delta_Dis,
				      d_Seropositivity,
				      d_internal_AR,
				      d_Prevalence,
				      d_end_state,
				      in_dmp,
				      R0_t,
				      regional_population_by_age,
				      E1_0, 
				      E2_0, 
				      I1_0, 
				      I2_0, 
				      global_params,
				      l_MIXMOD_ADJUSTED2,
				      vaccination,
				      vbooster,
				      vfourth);

  gsl_matrix_free(S);
  gsl_matrix_free(E_1);
  gsl_matrix_free(E_2);
  gsl_matrix_free(I_1);
  gsl_matrix_free(I_2);
  gsl_matrix_free(R_pos);
  gsl_matrix_free(R_neg);
  gsl_matrix_free(W);
  gsl_matrix_free(WS);
}

void fn_transmission_model(regional_model_params in_dmp,
			   global_model_instance_parameters global_params,
			   const gsl_vector* regional_population_by_age,
			   model_statistics& mod_stats,
			   const rtmData* vaccination,
			   const rtmData* vbooster,
			   const rtmData* vfourth)
{

  // WRITE DOWN THE FLUCTUATION IN R_0 OVER TIME
  gsl_vector* R0_t = gsl_vector_alloc(in_dmp.d_R0_phase_differences->size);
  Seasonal_R0_function(R0_t, in_dmp.l_R0_init, in_dmp.l_R0_Amplitude, in_dmp.d_R0_phase_differences);
  // ACCOUNT FOR THE FLUCTUATING BETA PARAMETER...
  gsl_vector_mul(R0_t, in_dmp.l_lbeta_rw);
  
  // INITIALISE THE SEIR MODEL
  gsl_vector* l_E1_0 = gsl_vector_alloc(NUM_AGE_GROUPS);
  gsl_vector* l_E2_0 = gsl_vector_alloc(NUM_AGE_GROUPS);
  gsl_vector* l_I1_0 = gsl_vector_alloc(NUM_AGE_GROUPS);
  gsl_vector* l_I2_0 = gsl_vector_alloc(NUM_AGE_GROUPS);

  // This will initialise the top, unvaccinated layer
  initialising_Deterministic_SE1E2I1I2R_AG_RF_CIN2(in_dmp, global_params,
						   *in_dmp.l_MIXMOD.evector_MIXMAT_normalised,
						   l_E1_0, l_E2_0, l_I1_0, l_I2_0);

  propagate_SEEIIR(in_dmp,
		   regional_population_by_age,
		   global_params,
		   R0_t,
		   l_E1_0,
		   l_E2_0,
		   l_I1_0,
		   l_I2_0,
		   mod_stats.d_NNI,
		   mod_stats.d_Delta_Dis,
		   mod_stats.d_seropositivity,
		   mod_stats.d_internal_AR,
		   mod_stats.d_prevalence,
		   mod_stats.d_end_state,
		   vaccination,
		   vbooster,
		   vfourth);

  ///////////////////////////////////////////////
  // free whatever isn't kept
  gsl_vector_free(R0_t);

  gsl_vector_free(l_E1_0);
  gsl_vector_free(l_E2_0);  
  gsl_vector_free(l_I1_0);
  gsl_vector_free(l_I2_0);

}





void fn_reporting_model(gsl_matrix* expected_counts, const gsl_matrix* NNI_transmission, const gsl_matrix* prop_symp, const gsl_matrix* severity_ratio,
			global_model_instance_parameters in_gmip, int max_days_data, const gsl_vector* population,
			gsl_vector* distribution_function,
			const int int_num_threads)
{

  // the input to the reporting model is the number of new infections.
  // these need to be aggregated up to the new delta_t, checking that this is
  // compatible with the delta_t used by the transmission model (i.e. one divides the other)
  // WANT TO STORE THE NUMBER OF NEW INFECTEDS, BUT ONLY AS A DAILY NUMBER, CONTRACT OVER DAYS...
  if((in_gmip.l_transmission_time_steps_per_day % in_gmip.l_reporting_time_steps_per_day) != 0)
    {
      printf("Incorrect specification of time steps per day for transmission model and reporting model\n");
      exit(2);
    }
  // allocate the memory we'll need
  gsl_matrix* NNI_rep_model_input = gsl_matrix_alloc(in_gmip.l_reporting_time_steps_per_day * in_gmip.l_duration_of_runs_in_days, NUM_AGE_GROUPS);
  // copy the input from the transmission model into a matrix whose elements we can edit
  gsl_matrix_memcpy(NNI_rep_model_input, NNI_transmission);

  // so far, have assumed that all infections are symptomatic. not so. multiply by the proportion of symptomatics
  gsl_matrix_mul_elements(NNI_rep_model_input, prop_symp);
  gsl_matrix_mul_elements(NNI_rep_model_input, severity_ratio);

  // similarly, multiply by proportion of individuals who will experience the severe event
  gsl_matrix_mul_elements(NNI_rep_model_input, severity_ratio);
  
  gsl_matrix* modelled_events = gsl_matrix_calloc(in_gmip.l_reporting_time_steps_per_day * in_gmip.l_duration_of_runs_in_days, NUM_AGE_GROUPS);

  // convolve the number of new infections over the (already calculated) delay distribution ##### NEED TO ADD A NUMBER OF THREADS CLAUSE!!!!!!
  // #pragma omp parallel for default(shared) schedule(static) num_threads(int_num_threads)
  for (int k = 0; k < (max_days_data * in_gmip.l_reporting_time_steps_per_day); ++k)
    {  
      for (int j = 0; j < NUM_AGE_GROUPS; ++j)
	{
	  double modelled_out = 0;
	  for (int m = 0; m <= FN_MIN(distribution_function->size - 1, k); ++m)
	    {
	      modelled_out += gsl_matrix_get(NNI_rep_model_input,k - m, j) * gsl_vector_get(distribution_function, m);
<<<<<<< HEAD
=======
	      // gsl_matrix_set(modelled_events,
	      // 		     k,
	      // 		     j,
	      // 		     gsl_matrix_get(modelled_events, k, j) + (gsl_matrix_get(NNI_rep_model_input,k - m, j) * gsl_vector_get(distribution_function, m))); // severity_ratio[k - m, j] HAS BEEN USED IN THE PAST - THIS GIVES THE PROPORTION OF SYMPTOMATIC CASES WHO 
>>>>>>> d57401911579f16086976ff69cbcc365b7cbd93d
	    }// FOR
	  gsl_matrix_set(modelled_events, k, j, modelled_out);
	}// FOR
    }// FOR

  // aggregate these counts at intervals of 1/(l_reporting_time_steps_per_day) to counts per day
  output_per_selected_period(in_gmip.l_reporting_time_steps_per_day, modelled_events, expected_counts);

  gsl_matrix_free(NNI_rep_model_input);
  gsl_matrix_free(modelled_events);

}


void fn_background_model(gsl_matrix* total_expected, gsl_matrix* viropositivity, // top row, vector whose values are set.
			 const gsl_matrix* H1N1_expected, const gsl_matrix* background_counts, const global_model_instance_parameters& in_gmip, const bool& Viro_flag, const double& sensitivity, const double& specificity)
{

  gsl_matrix_memcpy(total_expected, H1N1_expected);
  gsl_matrix_add(total_expected, background_counts); // TE = S + B

  if(gsl_matrix_min(total_expected) < 0)
    gsl_matrix_set_all(viropositivity, -1.0);
  else {
    if(Viro_flag)
      {

	gsl_matrix* true_vironegativity = gsl_matrix_alloc(in_gmip.l_duration_of_runs_in_days, NUM_AGE_GROUPS);

	gsl_matrix_memcpy(viropositivity, H1N1_expected); // VP = S
	gsl_matrix_memcpy(true_vironegativity, background_counts); // VN = B

	// divide the H1N1_expected matrix by the matrix sum to get the naive positivity
	gsl_matrix_div_elements(viropositivity, total_expected); // VP = S / (S + B)
	gsl_matrix_div_elements(true_vironegativity, total_expected); // VN = B / (S + B)

	// scale by the sensitivity and the specificity
	gsl_matrix_scale(viropositivity, sensitivity);
	gsl_matrix_scale(true_vironegativity, 1 - specificity);

	// add them to get the virological positivity
	gsl_matrix_add(viropositivity, true_vironegativity);

	gsl_matrix_free(true_vironegativity);

      }

  }

}

// LOG-LIKELIHOOD OF BINOMIAL "POSITIVITY" DATA
double fn_log_lik_positivity(const gsl_matrix* mat_nsample,
			     const gsl_matrix* mat_npositive,
			     const gsl_matrix* mat_model_positivity)
{

  double lfx = 0.0;
  double x, n, p; // THESE REALLY SHOULD BE INTEGERS, BUT CURRENTLY THERE ARE NO GUARANTEES OF INTEGER COUNTS DUE TO THE NON-SPLITTING OF THE UNDER 4S.

  // gsl_matrix_const_view mat_prob_check = gsl_matrix_const_submatrix(mat_model_positivity, lbounds.lower - 1, 0, lbounds.upper - lbounds.lower + 1, mat_model_positivity->size2);

  if(gsl_matrix_min(mat_model_positivity) < 0 || gsl_matrix_max(mat_model_positivity) > 1)
    lfx += GSL_NEGINF;
  else {

    for(int inti = 0; inti < mat_nsample->size1; inti++)
      {

	for(int intj = 0; intj < mat_nsample->size2; intj++)
	  {
	    x = gsl_matrix_get(mat_npositive, inti, intj);
	    n = gsl_matrix_get(mat_nsample, inti, intj);
	    p = gsl_matrix_get(mat_model_positivity, inti, intj);

	    lfx += (p == 0.0 || p == 1.0 ? (p == 0 ? (x == 0 ? 0 : GSL_NEGINF) : (x == n ? 0 : GSL_NEGINF)) : (x * gsl_sf_log(p)) + ((n - x) * gsl_sf_log(1 - p)));

	  }

      }

  }

  return lfx;
}

// LOG-LIKELIHOOD OF POISSON "COUNT" DATA
double fn_log_lik_countdata(const gsl_matrix* mat_counts,
			    const gsl_matrix* mat_expected_counts)
{

  double lfx = 0.0;
  double x, mu;

  // Check there are no negative expected counts.
  // gsl_matrix_const_view mat_count_check = gsl_matrix_const_submatrix(mat_expected_counts, lbounds.lower - 1, 0, lbounds.upper - lbounds.lower + 1, mat_expected_counts->size2);

  if(gsl_matrix_min(mat_expected_counts) < 0)
    lfx += GSL_NEGINF;
  else{

    for(int inti = 0; inti < mat_counts->size1; inti++)
      {
	for(int intj = 0; intj < mat_counts->size2; intj++)
	  {
	    mu = gsl_matrix_get(mat_expected_counts, inti, intj);
	    x = gsl_matrix_get(mat_counts, inti, intj);
	    if((mu == 0) && (x == 0)){
	      lfx += 0;
	    } else if((mu == 0) && (x > 0)){
	      lfx = GSL_NEGINF;
	    } else {
	      lfx += ((x * gsl_sf_log(mu)) - mu);
	    }
	  }
      }

  }
  return lfx;
}


double fn_log_lik_negbindata(const gsl_matrix* mat_counts,
			     const gsl_matrix* mat_expected_counts,
			     const gsl_matrix* mat_dispersion_params)
{

  double lfx = 0.0;
  double mu, eta;
  int x;

  // Check there are no negative expected counts.
  // gsl_matrix_const_view mat_count_check = gsl_matrix_const_submatrix(mat_expected_counts, lbounds.lower - 1, 0, lbounds.upper - lbounds.lower + 1, mat_expected_counts->size2);

  if(gsl_matrix_min(mat_expected_counts) < 0)
    lfx += GSL_NEGINF;
  else{
    for(int inti = 0; inti < mat_counts->size1; inti++)
      {
	for(int intj = 0; intj < mat_counts->size2; intj++)
	  {
	    mu = gsl_matrix_get(mat_expected_counts, inti, intj);
	    x = (int) round(gsl_matrix_get(mat_counts, inti, intj));
	    if((mu == 0) && (x == 0)){
	      lfx += 0;
	    } else if((mu == 0) && (x != 0)){
	      lfx += GSL_NEGINF;
	    } else {
	      eta = gsl_matrix_get(mat_dispersion_params, inti, intj);
	      if(eta > 1.5e-08){
		double r = mu / eta;
		lfx += gsl_sf_lngamma(x + r) - gsl_sf_lngamma(r);
		double p = 1 - (1 / (eta + 1));
		lfx += ((r * gsl_sf_log(1 - p)) + (x * gsl_sf_log(p)));
	      } else 
		lfx += (x * gsl_sf_log(mu)) - mu; // Dispersion is so small, likelihood is practically Poisson. Any differences in likelihood would be small in comparison to prior/proposal differences and would lead to acceptance anyway.
	    }
	  }
      }
  }
  
  return lfx;
}

double fn_log_lik_loggaussian_fixedsd(const gsl_matrix* mat_data,
				      const gsl_matrix* mat_emeans,
				      const gsl_matrix* mat_sds)
{
  double lfx = 0.0;
  for(int inti = 0; inti < mat_data->size1; inti++)
    {
      for(int intj = 0; intj < mat_data->size2; intj++)
	{
	  double sd = gsl_matrix_get(mat_sds, inti, intj);
	  if(sd > DBL_EPSILON){
	    double mu = gsl_sf_log(gsl_matrix_get(mat_emeans, inti, intj));
	    double tmp = (gsl_matrix_get(mat_data, inti, intj) - mu) / sd;
	    lfx -= (tmp * tmp / 2);
	  }
	}
    }
  return lfx;
}



void fn_log_likelihood_global(glikelihood& llhood,
		       Region* meta_region,
		       int region,
		       bool flag_update_transmission_model,
		       bool flag_update_reporting_model,
		       bool flag_update_GP_likelihood,
		       bool flag_update_Hosp_likelihood,
		       bool flag_update_Viro_likelihood,
		       bool flag_update_Sero_likelihood,
		       bool flag_update_Prev_likelihood,
		       const global_model_instance_parameters &gmip,
		       gslVector& gp_distribution_function,
		       gslVector& hosp_distribution_function) {


  
  // PASSING A VALUE OF 0 FOR THE ARGUMENT region EVALUATES THE
  // SPECIFIED COMPONENTS OF THE LIKELIHOOD OVER ALL REGIONS
  // int low_region = (region == -1) ? 0 : region;
  // int hi_region = (region == -1) ? gmip.l_num_regions : region + 1;
  // double temp_log_likelihood;
  // double lfx_increment = 0.0;

  // #ifdef USE_THREADS
  //  int num_parallel_regions = FN_MIN(omp_get_num_procs(), hi_region - low_region);
  //  int num_subthread_teams = ceil(((double) omp_get_num_procs()) / ((double) num_parallel_regions));
  // #else
  //  int num_subthread_teams = 1;
  // #endif

  

  // #pragma omp parallel for private(temp_log_likelihood) default(shared) num_threads(num_parallel_regions) schedule(static) reduction(+:lfx_increment)

  double lfx_sum = 0;
  
#pragma omp parallel for default(shared) schedule(static) reduction(+:lfx_sum)
  for(int int_region = 0; int_region < gmip.l_num_regions; int_region++) {
    
    fn_log_likelihood_region(llhood.rlik[int_region], meta_region, int_region, flag_update_transmission_model, flag_update_reporting_model, flag_update_GP_likelihood, flag_update_Hosp_likelihood, flag_update_Viro_likelihood, flag_update_Sero_likelihood, flag_update_Prev_likelihood, gmip, gp_distribution_function, hosp_distribution_function);
    
    lfx_sum += llhood.rlik[int_region].region_lfx;
  }
  
  llhood.total_lfx += lfx_sum;
}


void fn_log_likelihood_region(rlikelihood& llhood,
			      Region* meta_region,
			      int int_region,
			      bool flag_update_transmission_model,
			      bool flag_update_reporting_model,
			      bool flag_update_GP_likelihood,
			      bool flag_update_Hosp_likelihood,
			      bool flag_update_Viro_likelihood,
			      bool flag_update_Sero_likelihood,
			      bool flag_update_Prev_likelihood,
			      const global_model_instance_parameters &gmip,
			      gslVector& gp_distribution_function,
			      gslVector& hosp_distribution_function)  {

  double lfx_region = 0;
  double temp_log_likelihood = 0;
  //double lfx_increment = 0;
      
  // Does the transmission model need to be re-evaluated?
  // (Not necessary when updating parameters of the reporting model)
  if(flag_update_transmission_model)
    fn_transmission_model(meta_region[int_region].det_model_params,
			  gmip,
			  meta_region[int_region].population,
			  meta_region[int_region].region_modstats,
			  meta_region[int_region].Vaccination_data,
			  meta_region[int_region].VBoosting_data,
			  meta_region[int_region].VFourth_data);
      
  // Having evaluated the transmission model, do we need to evaluate the seropositivity likelihood
  if(gmip.l_Sero_data_flag &&
     (flag_update_transmission_model || flag_update_Sero_likelihood)
     ){ // HERE!!! NEED TO ADD A CONDITION INTO HERE && (flag_update_transmission_model || new_flag_for_updating_serology)
    // Get the seropositivity at the HI>32 level - subtract the initial portion who are positive at HI>8 but not HI>32
    // This proportion is age, but not time dependent.
    gsl_matrix* test_positivity = gsl_matrix_alloc(meta_region[int_region].region_modstats.d_internal_AR->size1, meta_region[int_region].region_modstats.d_internal_AR->size2);
    gsl_matrix* test_sens_scaling = gsl_matrix_alloc(test_positivity->size1, test_positivity->size2);
    gsl_matrix_memcpy(test_positivity, meta_region[int_region].region_modstats.d_internal_AR);
    gsl_matrix_memcpy(test_sens_scaling, meta_region[int_region].det_model_params.l_sero_sensitivity);
	
    gsl_vector* prop_immune_baseline_nonseropositive = gsl_vector_alloc(meta_region[int_region].det_model_params.l_init_prop_sus->size);
    gsl_vector* temp_vec = gsl_vector_alloc(meta_region[int_region].det_model_params.l_init_prop_sus->size);
    gsl_vector_memcpy(prop_immune_baseline_nonseropositive, meta_region[int_region].det_model_params.l_init_prop_sus);
    gsl_vector_memcpy(temp_vec, meta_region[int_region].det_model_params.l_init_prop_sus_HI_geq_32);
    gsl_vector_add_constant(prop_immune_baseline_nonseropositive, -1.0);
    gsl_vector_add_constant(temp_vec, -1.0);
    gsl_vector_mul(prop_immune_baseline_nonseropositive, temp_vec);
    gsl_vector_free(temp_vec);

    for(int int_t = 0; int_t < test_positivity->size1; int_t++)
      {
	gsl_vector_view seropos_row = gsl_matrix_row(test_positivity, int_t);
	gsl_vector_sub(&seropos_row.vector, prop_immune_baseline_nonseropositive);
      }
    gsl_vector_free(prop_immune_baseline_nonseropositive);

    // ** Some account for test sensitivity and specificity. Will have to move from here
    // ** if these two quantities are ever to be allowed to vary by time, region or age.
    gsl_matrix_add(test_sens_scaling, meta_region[int_region].det_model_params.l_sero_specificity);
    gsl_matrix_add_constant(test_sens_scaling, -1.0);
    gsl_matrix_mul_elements(test_positivity, test_sens_scaling);
    gsl_matrix_sub(test_positivity, meta_region[int_region].det_model_params.l_sero_specificity);
    gsl_matrix_add_constant(test_positivity, 1.0);
    gsl_matrix_free(test_sens_scaling);
	
    // ** Is there any missing data - if dataset is of dimension less than the number of strata
    if(test_positivity->size2 != meta_region[int_region].Serology_data->getDim2()){
      // ** Yes: Aggregate seropositivities using a weighted mean
      gsl_matrix* weighted_positivity = gsl_matrix_calloc(meta_region[int_region].Serology_data->getDim1(),
							  meta_region[int_region].Serology_data->getDim2());
      for(int int_t = 0; int_t < test_positivity->size1; int_t++)
	{
	  gsl_vector_view seropos_row = gsl_matrix_row(test_positivity, int_t);
	  gsl_vector_mul(&seropos_row.vector, meta_region[int_region].Serology_data->access_weights());
	  gsl_vector_view seropos_out_row = gsl_matrix_row(weighted_positivity, int_t);
	  R_by_sum_gsl_vector_mono_idx(&seropos_out_row.vector,
				       NULL,
				       &seropos_row.vector,
				       meta_region[int_region].Serology_data->access_groups());
	}
      temp_log_likelihood = meta_region[int_region].Serology_data->lfx(weighted_positivity, NULL);
      gsl_matrix_free(weighted_positivity);  
    } else // ** Just calculate the likelihood based on the already computed positivities
      temp_log_likelihood = meta_region[int_region].Serology_data->lfx(test_positivity, NULL);
    // CCS
    lfx_region += (temp_log_likelihood - llhood.Sero_lfx);
    llhood.Sero_lfx = temp_log_likelihood;
    gsl_matrix_free(test_positivity);
  }

  if(gmip.l_Prev_data_flag && (flag_update_transmission_model || flag_update_Prev_likelihood))
    {
      gsl_matrix* summed_prevalence = gsl_matrix_calloc(meta_region[int_region].region_modstats.d_prevalence->size1,
							meta_region[int_region].Prevalence_data->getDim2());
      if(gsl_matrix_min(meta_region[int_region].region_modstats.d_prevalence) > DBL_EPSILON)
	{
	  // ** Is there any missing data - if dataset is of dimension less than the number of strata
	  if(meta_region[int_region].region_modstats.d_prevalence->size2 != meta_region[int_region].Prevalence_data->getDim2()){
	    // ** Aggregate total prevalence (we work with numbers prevalent, not a proportion)
	    for(int int_k = 0; int_k < summed_prevalence->size1; int_k++)
	      {
		gsl_vector_view full_strata_prev = gsl_matrix_row(meta_region[int_region].region_modstats.d_prevalence, int_k);
		gsl_vector_view agg_strata_prev = gsl_matrix_row(summed_prevalence, int_k);
		R_by_sum_gsl_vector_mono_idx(&agg_strata_prev.vector,
					     NULL,
					     &full_strata_prev.vector,
					     meta_region[int_region].Prevalence_data->access_groups());
	      }
	  } else gsl_matrix_memcpy(summed_prevalence, meta_region[int_region].region_modstats.d_prevalence);
	  temp_log_likelihood = meta_region[int_region].Prevalence_data->meld_lfx(summed_prevalence);
	  // CCS
	  lfx_region += (temp_log_likelihood - llhood.Prev_lfx);
	  //if (debug) cout << int_region << " Prev prev " << gsl_vector_get(*llhood.Prev_lfx, int_region) << " curr " << temp_log_likelihood << endl;
	  llhood.Prev_lfx = temp_log_likelihood;
	}
      gsl_matrix_free(summed_prevalence);
    }
      
  double lfx_sub_increment = 0.0;
  // #pragma omp parallel default(shared) num_threads(2) reduction(+:lfx_sub_increment)
  //       {

  //        	// Have assigned two processors to this scope - first processor deals with GP related data
  //        	if(omp_get_thread_num() == 0)
  //        	  {
  if((flag_update_GP_likelihood || flag_update_Viro_likelihood) && (((bool) gmip.l_GP_consultation_flag) || ((bool) gmip.l_Viro_data_flag)))
    {
      if(flag_update_reporting_model)
	{
	  fn_reporting_model(meta_region[int_region].region_modstats.d_H1N1_GP_Consultations,
			     meta_region[int_region].region_modstats.d_Delta_Dis,
			     meta_region[int_region].det_model_params.l_pr_symp,
			     meta_region[int_region].det_model_params.l_pr_onset_to_GP,
			     gmip,
			     gmip.l_GP_likelihood.upper,
			     meta_region[int_region].population,
			     *gp_distribution_function /*, num_subthread_teams */ );
	}
	
      // ADD THE BACKGROUND AND CALCULATE THE UPDATED POSITIVITY
      fn_background_model(meta_region[int_region].region_modstats.d_Reported_GP_Consultations,
			  meta_region[int_region].region_modstats.d_viropositivity,
			  meta_region[int_region].region_modstats.d_H1N1_GP_Consultations,
			  meta_region[int_region].det_model_params.l_background_gps_counts,
			  gmip,
			  true,
			  meta_region[int_region].det_model_params.l_sensitivity,
			  meta_region[int_region].det_model_params.l_specificity);
	  
      // Q-SURVEILLANCE DATA DOESN'T HAVE 100% COVERAGE OF THE POPULATION.
      // SCALE THE EXPECTED COUNTS BY THE DAILY % COVERAGE BY AGE (the latter is done inside likelihood functions
      if(gmip.l_GP_consultation_flag || flag_update_GP_likelihood)
	{

	  gsl_matrix_mul_elements(meta_region[int_region].region_modstats.d_Reported_GP_Consultations,
				  meta_region[int_region].det_model_params.l_day_of_week_effect);

	  if(gmip.l_GP_consultation_flag && flag_update_GP_likelihood)
	    {
	      if(gsl_matrix_min(meta_region[int_region].region_modstats.d_Reported_GP_Consultations) >= 0)
		{
		  gsl_matrix* mu_gp_counts = gsl_matrix_alloc(meta_region[int_region].region_modstats.d_Reported_GP_Consultations->size1,
							      meta_region[int_region].GP_data->getDim2());
		  if(meta_region[int_region].region_modstats.d_Reported_GP_Consultations->size2 != meta_region[int_region].GP_data->getDim2()){
		    // Yes: we have missing data
		    for(int int_t = 0; int_t < mu_gp_counts->size1; int_t++)
		      {
			gsl_vector_view full_strata_counts = gsl_matrix_row(meta_region[int_region].region_modstats.d_Reported_GP_Consultations, int_t);
			gsl_vector_view agg_strata_counts = gsl_matrix_row(mu_gp_counts, int_t);
			R_by_sum_gsl_vector_mono_idx(&agg_strata_counts.vector,
						     NULL,
						     &full_strata_counts.vector,
						     meta_region[int_region].GP_data->access_groups());
		      }
		  } else gsl_matrix_memcpy(mu_gp_counts, meta_region[int_region].region_modstats.d_Reported_GP_Consultations);
		  data_type dlfx = meta_region[int_region].GP_data->get_likelihood_type();
		  if(dlfx == cPOISSON_LIK)
		    temp_log_likelihood = meta_region[int_region].GP_data->lfx(mu_gp_counts, NULL);
		  else if(dlfx == cNEGBIN_LIK)
		    temp_log_likelihood = meta_region[int_region].GP_data->lfx(mu_gp_counts,
									       meta_region[int_region].det_model_params.l_gp_negbin_overdispersion);
		  else {
		    perror("Unrecognised likelihood selected\n");
		    exit(2);
		  }
		  lfx_sub_increment += (temp_log_likelihood - llhood.GP_lfx);
		  llhood.GP_lfx = temp_log_likelihood;
		  gsl_matrix_free(mu_gp_counts);
		}
	      else {
		lfx_sub_increment += GSL_NEGINF;
		llhood.GP_lfx = GSL_NEGINF;
	      }
	    }
	}
      if(gmip.l_Viro_data_flag && flag_update_Viro_likelihood)
	{
	  // Viropositivity already calculated by each modelled strata.. if there is missing data we need to take a weighted average of the positivities.
	  // Where the weights are proportional to the number of consultations in each strata.

	  // If missing data...
	  if(meta_region[int_region].region_modstats.d_viropositivity->size2 != meta_region[int_region].Virology_data->getDim2()){

	    // Aggregate viropositivities using a weighted mean
	    gsl_matrix* weighted_positivity = gsl_matrix_alloc(meta_region[int_region].Virology_data->getDim1(),
							       meta_region[int_region].Virology_data->getDim2());

	    for(int int_t = 0; int_t < meta_region[int_region].Virology_data->getDim1(); int_t++)
	      {
		gsl_vector_view viropos_out_row = gsl_matrix_row(weighted_positivity, int_t);
		gsl_vector_view viropos_row = gsl_matrix_row(meta_region[int_region].region_modstats.d_viropositivity, int_t);
		gsl_vector_view gp_row = gsl_matrix_row(meta_region[int_region].region_modstats.d_Reported_GP_Consultations, int_t);
		meta_region[int_region].Virology_data->data_population_sizes(&gp_row.vector);
		gsl_vector_mul(&viropos_row.vector, meta_region[int_region].Virology_data->access_weights());
		R_by_sum_gsl_vector_mono_idx(&viropos_out_row.vector,
					     NULL,
					     &viropos_row.vector,
					     meta_region[int_region].Virology_data->access_groups());
	      }
	    temp_log_likelihood = meta_region[int_region].Virology_data->lfx(weighted_positivity, NULL);
	    gsl_matrix_free(weighted_positivity);
	  } else
	    temp_log_likelihood = meta_region[int_region].Virology_data->lfx(meta_region[int_region].region_modstats.d_viropositivity, NULL);

	  lfx_sub_increment += (temp_log_likelihood - llhood.Viro_lfx);
	  llhood.Viro_lfx = temp_log_likelihood;
	      
	}
	  
    }
  // 	  }
      
      
  // Have second processor assigned to hospitalisations
  //        	if(omp_get_thread_num() == 1)
  //        	  {
  if(flag_update_Hosp_likelihood && ((bool) gmip.l_Hospitalisation_flag))
    {
      fn_reporting_model(meta_region[int_region].region_modstats.d_Reported_Hospitalisations,
			 meta_region[int_region].region_modstats.d_Delta_Dis,
			 meta_region[int_region].det_model_params.l_pr_symp,
			 meta_region[int_region].det_model_params.l_pr_onset_to_Hosp,
			 gmip,
			 gmip.l_Hosp_likelihood.upper,
			 meta_region[int_region].population,
			 *hosp_distribution_function);
	  
      gsl_matrix* mu_hosp_counts = gsl_matrix_alloc(meta_region[int_region].region_modstats.d_Reported_Hospitalisations->size1,
						    meta_region[int_region].Hospitalisation_data->getDim2());

      if(gsl_matrix_min(meta_region[int_region].region_modstats.d_Reported_Hospitalisations) >= 0)
	{

	  if(meta_region[int_region].region_modstats.d_Reported_Hospitalisations->size2 != meta_region[int_region].Hospitalisation_data->getDim2()){
	    // Yes: again we have missing data
	    for(int int_t = 0; int_t < mu_hosp_counts->size1; int_t++)
	      {
		gsl_vector_view full_strata_counts = gsl_matrix_row(meta_region[int_region].region_modstats.d_Reported_Hospitalisations, int_t);
		gsl_vector_view agg_strata_counts = gsl_matrix_row(mu_hosp_counts, int_t);
		R_by_sum_gsl_vector_mono_idx(&agg_strata_counts.vector,
					     NULL,
					     &full_strata_counts.vector,
					     meta_region[int_region].Hospitalisation_data->access_groups());
	      }
	  } else gsl_matrix_memcpy(mu_hosp_counts, meta_region[int_region].region_modstats.d_Reported_Hospitalisations);
	  data_type dlfx = meta_region[int_region].Hospitalisation_data->get_likelihood_type();
	  if(dlfx == cPOISSON_LIK)
	    temp_log_likelihood = meta_region[int_region].Hospitalisation_data->lfx(mu_hosp_counts,
										    NULL);
	  else if(dlfx == cNEGBIN_LIK)
	    temp_log_likelihood = meta_region[int_region].Hospitalisation_data->lfx(mu_hosp_counts,
										    meta_region[int_region].det_model_params.l_hosp_negbin_overdispersion);
	  else {
	    perror("Unrecognised likelihood selected\n");
	    exit(2);
	  }
	} else temp_log_likelihood = GSL_NEGINF;
	  
      lfx_sub_increment += (temp_log_likelihood - llhood.Hosp_lfx);

      //if (debug) cout << int_region << " Hosp prev " << gsl_vector_get(*llhood.Hosp_lfx, int_region) << " curr " << temp_log_likelihood << endl;

      llhood.Hosp_lfx = temp_log_likelihood;

	  
      gsl_matrix_free(mu_hosp_counts);	  
    }

  lfx_region += lfx_sub_increment;

  // Save the regional likelihood component for AMGS blocks
  llhood.region_lfx = lfx_region;

  // cout << "Set region " << int_region << ": " << llhood.region_lfx[int_region] << endl;

  // lfx_increment += lfx_region;
  // if (debug) cout << int_region << " end incr " << lfx_increment << endl;

  // }
  
  // llhood.total_lfx += lfx_increment;
  
}
