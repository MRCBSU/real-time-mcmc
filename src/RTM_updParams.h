#ifndef HEADER_updParams_
#define HEADER_updParams_

#include <exception>
#include <map>
#include <string>
#include <vector>

#include "gslWrapper.h"
#include "RTM_StructDefs.h"

// Class contains a single model parameter
// updParam = Updateable model parameter

class updParam {

public:
  
  // Constructor for specifying default values
  // Flag is a binary 7-bit value consisting of the following 7 flags:
  // transmission, reporting, gp, hosp, sero, viro, prev data likelihood
  updParam(string name, double value, int flags);

  // Constructor for current init by copying from an existing param
  // TODO: Remove this when input code improved
  updParam(updateable_model_parameter &in, const int num_instances);


  std::string param_name;
  // CCS: TODO REMOVE
  // This is used to allow easy construction of the main param vectors
  // Will be set to null once input is complete, and not used throughout M-H
  gslVector param_value;
  unsigned int valuePtr;
  unsigned int param_index;
  gslVector proposal_value;
  gslVectorInt prior_distribution; // THE DISTRIBUTION OF THE INDIVIDUAL COMPONENTS
  bool flag_hyperprior;
  double log_prior_dens;
  double proposal_log_prior_dens;
  bool flag_update; // TRUE: IF PRIOR_DISTRIBUTION[i] != cDETERM FOR ALL i
  std::vector<gslVector> prior_params;
//  gsl_vector **prior_params; // WHAT THESE PARAMETERS REPRESENT DEPENDS ON THE prior_distribution MEMEMBER OF THE STRUCTURE. EACH ROW CORRESPONDS TO EACH OF THE DISTRIBUTIONS IN THE prior_distribution VECTOR
  gslVector proposal_variances;
  mvnorm *prior_multivariate_norm;
  regression_def map_to_regional;
  gslVector posterior_mean;
  gslVector posterior_sumsq;
  bool flag_transmission_model; // DOES THE NUMBER OF NEW INFECTEDS NEED TO RECALCULATED WHEN UPDATING THE PARAMETER THIS PARAMETER
  bool flag_reporting_model; // DO THE REPORTING DELAYS NEED TO BE RECALCULATED WHEN UPDATING WHEN UPDATING THIS PARAMETER
  bool flag_GP_likelihood; // DOES THE LIKELIHOOD FOR G.P. CONSULTATION NEED TO BE RECALCULATED WHEN UPDATING THIS PARAMETER
  bool flag_Hosp_likelihood; // DOES THE LIKELIHOOD FOR THE HOSPITALISATIONS NEED TO BE RECALCULATED WHEN UPDATING THIS PARAMETER
  bool flag_Sero_likelihood; // DOES THE LIKELIHOOD FOR THE SEROEPIDEMIOLOGY DATA NEED TO BE RECALCULATED WHEN UPDATING THIS PARAMETER
  bool flag_Viro_likelihood; // DOES THE LIKELIHOOD FOR THE VIROLOGY (POSITIVITY) DATA NEED TO BE RECALCULATED WHEN UPDATING THIS PARAMETER
  bool flag_Prev_likelihood; // DOES THE LIKELIHOOD FOR THE PREVALENCE DATA/ESTIMATES NEED TO BE RECALCULATED WHEN UPDATING THIS PARAMETER
  bool flag_any_child_nodes; // TRUE IF ANY FLAG_CHILD_NODES ARE TRUE. FALSE OTHERWISE
  std::vector<bool> flag_child_nodes; // FLAG FOR EACH OF THE OTHER PARAMETERS OF THE MODEL INDICATING WHETHER THEY ARE CHILD NODES OF THE CURRENT NODE
};



// Class for managing the set of params
// This is a 'static' class, i.e. we don't need an instance object
class updParamSet {

public:
  
  // Define parameter names as used by code, particularly the region params update
  // Need to be specified by eg. "ump::PROP_SUS_HYPER"
  enum paramIndex {
    EGR_HYPER, LPL0_HYPER, PROP_SUS_HYPER,
    LBETA_RW_SD, GP_OVERDISP, HOSP_OVERDISP,
    ALP, AIP, AR1,
    REL_INFECT, PROP_SYMP, CONTACT,
    LBETA_RW, R0_AMP, R0_PEAKDAY,
    EGR, LPL0, PROP_SUS,
    PROP_HI_GEQ_32, PROP_GP, PROP_HOSP,
    PROP_DEATH, IMPORTATION, BGR,
    SENS, SPEC, DOW_EFFECTS,
    SSENS, SSPEC,
    // This value must always be last. It gives the total number of params.
    UMP_FINAL_NUM_PARAMS
  };

  static const std::map<std::string, int> nameMap;

  // Lookup vector. For i = one of the paramIndices, lookup[i] = index:
  // - if (i > 0) param is globalParams[index - 1]
  // - if (i < 0) param is localParams[-(index - 1)
  // - if (i == 0) error (param doesn't exist or hasn't been created yet)
  
  static std::vector<int> lookupVec;
  
  static std::vector<updParam> globalParams;
  static std::vector<updParam> localParams;

  static gslVector globalParamValues;
  static std::vector<gslVector> localParamValues;

  // 'Constructor'
  static void init(int numRegions) {
    lookupVec.resize(UMP_FINAL_NUM_PARAMS);
    localParamValues.resize(numRegions);
  }

  // regional: true = local param, false = global
  static void insertAndRegister(updateable_model_parameter& inPar, int num_instances, bool regional);


  static void populateVectors(int regions);
  

  // Return a reference to the given parameter
  static updParam& lookup(paramIndex index) {
    int i = lookupVec[index];
    if (i > 0)
      return globalParams[i-1];
    else if (i < 0)
      return localParams[-1 * (i-1)];
    else // i == 0
      throw std::invalid_argument("lookup index is not a valid parameter index");
  }


  // Return a gsl_vector_view for the given parameter, i.e. a sub-vector
  // containing only this parameter's values
  static gsl_vector_const_view value(paramIndex index, int region = 0) {

    int i = lookupVec[index];
    updParam &par = lookup(index);

    // TODO: Better way of getting length of subvector
    if (i > 0)
      return gsl_vector_const_subvector(*globalParamValues, par.valuePtr, par.param_value.size());
    else if (i < 0) {
      return gsl_vector_const_subvector(*(localParamValues[region]), par.valuePtr, par.param_value.size());
    }
    else {
      throw std::invalid_argument("invalid index");
    }
  }
  
  // - - - - - - - - - - - - - - - - -
  
  // Future work: The mod_pars.txt file should be read from here. The existing
  // input code is difficult to read and edit. But it works, so leave it for now.
  // At present, these two methods are unused.
  // Current defaults in RTM_StructAssign
  updParam defaultParam(std::string &name) {
    if (name == "exponential_growth_rate_hyper")
      return updParam(name, 0.15, 0b0000000);
    else if (name == "l_p_lambda_0_hyper")
      return updParam(name, -15.0, 0b0000000);
    // else if (...)
    else
      throw std::invalid_argument("");
    }
  
  // Takes a parameter name and returns a parameter containing the relevant 
  // defaults. Setting defaults like this makes them easier to read and modify.
  void readInputFile();

  
};


// Equivalent of a typedef to shorten name
// Does this simplify or increase confusion?
using upd = updParamSet;








/*

class updateParamSet {
  
  std::map<std::string, int> nameMap = {
    { "exponential_growth_rate_hyper", 0 },
    { "l_p_lambda_0_hyper" 1 },
    { "prop_susceptible_hyper", 2 },
    { "log_beta_rw_sd", 3 },
    { "gp_negbin_overdispersion", 4 },
    { "hosp_negbin_overdispersion", 5 },
    { "latent_period", 6 },
    { "infectious_period", 7 },
    { "r1_period", 8 },
    { "relative_infectiousness", 9 },
    { "prop_symptomatic", 10 },
    { "contact_parameters", 11 },
    { "log_beta_rw", 12 },
    { "R0_amplitude_kA", 13 },
    { "R0_seasonal_peakday", 14 },
    { "exponential_growth_rate", 15 },
    { "log_p_lambda_0", 16 },
    { "prop_susceptible", 17 },
    { "prop_HI_32_to_HI_8", 18 },
    { "prop_case_to_GP_consultation", 19 },
    { "prop_case_to_hosp", 20 },
    { "prop_case_to_death", 21 },
    { "importation_rates", 22 },
    { "background_GP", 23 },
    { "test_sensitivity", 24 },
    { "test_specificity", 25 },
    { "day_of_week_effects", 26 },
    { "sero_test_sensitivity", 27 },
    { "sero_test_specificity", 28 }
  };

  
  
  updateParam &operator[](const std::string &name) {
    
  }

  updateParam &operator[](const int &id) {

  }
}

*/
#endif // HEADER_updParams_
