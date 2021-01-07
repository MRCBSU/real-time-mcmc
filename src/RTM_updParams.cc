#include "RTM_updParams.h"

// Note: Specifying bitmasks with 0b is a C++14 feature, but is a compiler
// extension in both gcc (4.3 or later) and clang (since at least 2013).
updParam::updParam(string name, double value, int flags)
  : param_name(name),
    param_value(value),
    flag_transmission_model (flags & 0b1000000),
    flag_reporting_model    (flags & 0b0100000),
    flag_GP_likelihood      (flags & 0b0010000),
    flag_Hosp_likelihood    (flags & 0b0001000),
    flag_Sero_likelihood    (flags & 0b0000100),
    flag_Viro_likelihood    (flags & 0b0000010),
    flag_Prev_likelihood    (flags & 0b0000001)
{
}


updParam::updParam(updateable_model_parameter &in, const int num_instances)
  : param_name (in.param_name),
    param_value (in.param_value),
    proposal_value (in.proposal_value),
    prior_distribution (in.prior_distribution),
    flag_hyperprior (in.flag_hyperprior),
    log_prior_dens (in.log_prior_dens),
    proposal_log_prior_dens (in.proposal_log_prior_dens), 
    flag_update (in.flag_update),
    // Length is the same as length of param_value
    prior_params (in.param_value->size),
    proposal_variances (in.proposal_variances),
    prior_multivariate_norm (in.prior_multivariate_norm),
    map_to_regional (in.map_to_regional),
    posterior_mean (in.posterior_mean),
    posterior_sumsq (in.posterior_sumsq),
    flag_transmission_model (in.flag_transmission_model), 
    flag_reporting_model (in.flag_reporting_model), 
    flag_GP_likelihood (in.flag_GP_likelihood),
    flag_Hosp_likelihood (in.flag_Hosp_likelihood),
    flag_Sero_likelihood (in.flag_Sero_likelihood), 
    flag_Viro_likelihood (in.flag_Viro_likelihood),
    flag_Prev_likelihood (in.flag_Prev_likelihood),
    flag_any_child_nodes (in.flag_any_child_nodes), 
    flag_child_nodes (num_instances) {
  
  // Some vectors in the input param are only allocated if flag_update is true
  if (flag_update) {
    // prior_params - gsl_vector**
    for (int i = 0; i < prior_params.size(); i++)
      prior_params[i] = gslVector(in.prior_params[i]);
    
  }
  // flag_child_nodes
  for (int i = 0; i < num_instances; i++)
    flag_child_nodes[i] = in.flag_child_nodes[i];
}




const std::map<std::string, int> updParamSet::nameMap = {
    { "exponential_growth_rate_hyper", EGR_HYPER },
    { "l_p_lambda_0_hyper", LPL0_HYPER },
    { "prop_susceptible_hyper", PROP_SUS_HYPER },
    { "log_beta_rw_sd", LBETA_RW_SD  },
    { "gp_negbin_overdispersion", GP_OVERDISP },
    { "hosp_negbin_overdispersion", HOSP_OVERDISP },
    { "latent_period", ALP },
    { "infectious_period", AIP },
    { "r1_period", AR1 },
    { "relative_infectiousness", REL_INFECT },
    { "prop_symptomatic", PROP_SYMP },
    { "contact_parameters", CONTACT },
    { "log_beta_rw", LBETA_RW },
    { "R0_amplitude_kA", R0_AMP },
    { "R0_seasonal_peakday", R0_PEAKDAY },
    { "exponential_growth_rate", EGR },
    { "log_p_lambda_0", LPL0 },
    { "prop_susceptible", PROP_SUS },
    { "prop_HI_32_to_HI_8", PROP_HI_GEQ_32 },
    { "prop_case_to_GP_consultation", PROP_GP },
    { "prop_case_to_hosp", PROP_HOSP },
    { "prop_case_to_death", PROP_DEATH },
    { "importation_rates", IMPORTATION },
    { "background_GP", BGR },
    { "test_sensitivity", SENS },
    { "test_specificity", SPEC },
    { "day_of_week_effects", DOW_EFFECTS },
    { "sero_test_sensitivity", SSENS },
    { "sero_test_specificity", SSPEC }
  };



void updParamSet::insertAndRegister(updateable_model_parameter& inPar, int num_instances, bool regional) {
  
  updParam newPar(inPar, num_instances);
  int index;
  
  if (regional) {
    localParams.push_back(std::move(newPar));
    index = localParams.size() * -1;
    
  } else {
    globalParams.push_back(std::move(newPar));
    index = globalParams.size();  
  }
  
  lookupVec[nameMap.at(inPar.param_name)] = index;

}


void updParamSet::populateVectors(int regions) {
    // Total number of local params
    int count = 0;
    for (auto &par : localParams)
      count += par.param_value.size();

    // TODO: Error check the two left devision calculations

    // Integer division
    count = count / regions;

    // Init the vectors for each region
    for (auto &vec : localParamValues)
      vec.create(count);

    // Copy

    // Insertion indices for each of the regional vectors
    std::vector<int> index(regions, 0);
    for (auto &par : localParams) {
      // Number of elts for each region
      // Integer division
      int regionSize = par.param_value.size() / regions;

      // Counter for reading from param_value
      int readIndex = 0;
      for (int r = 0; r < regions; r++)
	for (int i = 0; i < regionSize; i++)
	  localParamValues[r][index[r]++] = par.param_value[readIndex++];
    }
    
    // Global params
    count = 0;
    for (auto &par : globalParams)
      count += par.param_value.size();
    
    globalParamValues.create(count);
    
    int outInd = 0;
    for (auto &par : globalParams)
      for (int i = 0; i < par.param_value.size(); i++)
	globalParamValues[outInd++] = par.param_value[i];
    
  }


std::vector<int> updParamSet::lookupVec;

std::vector<updParam> updParamSet::globalParams;
std::vector<updParam> updParamSet::localParams;

gslVector updParamSet::globalParamValues;
std::vector<gslVector> updParamSet::localParamValues;
