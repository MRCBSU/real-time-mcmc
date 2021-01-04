#ifndef HEADER_Params_
#define HEADER_Params_

#include <vector>
#include <string>

#include "gslWrapper.h"
#include "RTM_StructDefs.h"


class updateableParam {

public:  

  // TODO: Remove num_proposed/accepted_moves and any other unused vals
  
  // TODO: Remove this constructor
updateableParam(updateable_model_parameter &in, const int num_instances)
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
    //number_accepted_moves (in.number_accepted_moves->size, 0),
    //number_proposed_moves (in.number_proposed_moves->size, 0), 
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
	number_accepted_moves = vector<int>(in.number_accepted_moves->size, 0);
	number_proposed_moves = vector<int>(in.number_proposed_moves->size, 0);
      
	// prior_params - gsl_vector**
	for (int i = 0; i < prior_params.size(); i++)
	  prior_params[i] = gslVector(in.prior_params[i]);

      }
      // flag_child_nodes
      for (int i = 0; i < num_instances; i++)
	flag_child_nodes[i] = in.flag_child_nodes[i];
    }


  
  std::string param_name;
  // CCS: TODO REMOVE
  // This is used to allow easy construction of the main param vectors
  // Will be set to null once input is complete, and not used throughout M-H
  gslVector param_value;
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
  std::vector<int> number_accepted_moves;
  std::vector<int> number_proposed_moves;
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


extern gslVector globalParamValues;
extern std::vector<gslVector> localParamValues;

extern std::vector<updateableParam> globalParams;
extern std::vector<updateableParam> localParams;


#endif // HEADER_Params_
