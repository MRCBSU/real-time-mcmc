#ifndef HEADER_updParams_
#define HEADER_updParams_

#include <exception>
#include <map>
#include <string>
#include <vector>

#include <gsl/gsl_rng.h>

#include "gslWrapper.h"
#include "RTM_StructDefs.h"


// Forward declaration so that updParam can use upd::paramIndex
class updParam;
class updParamSet;

class updParamBlock {
public:
  gslVector vals;
  gslVector proposal;
  gslVector mu;
  double beta;
  gslMatrix sigma;
  //double logPosterior;
  double laccept;
  int acceptLastMove;  // Int for use in adaptive update
  int numAccept;

  void calcProposal(gsl_rng *rng);
  void calcAccept(updParamSet &paramSet, bool global);
  void doAccept(gsl_rng *rng);
  void adaptiveUpdate(int iter);
};

// Class for managing the set of params
class updParamSet {
public:

  // List of updateable parameter names
  // Need to be specified by class name, eg "ump::PROP_SUS_HYPER"
  enum paramIndex {
    INVALID=-1, // For parent vector when there are no parents
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
    // This must always be last. It gives total number of params.
    UMP_FINAL_NUM_PARAMS
  };
  
  updParamSet() :
    params(UMP_FINAL_NUM_PARAMS) { }

  // Mapping of strings to names above
  // Needed for copying param from existing input code; remove when tidied up
  static const std::map<std::string, paramIndex> nameMap;

  // Is this necessary? So far nameMap only used inside this class
  paramIndex getID(std::string name) {
    // TODO: Throw error if string not there
    return nameMap.at(name);
  }
    
  // Lookup vector. For i = one of the paramIndices, lookup[i] = index:
  // - if (i > 0) param is globalParams[index - 1]
  // - if (i <std 0) param is localParams[-(index - 1)
  // - if (i == 0) error (param doesn't exist or hasn't been created yet)
  //std::vector<int> lookupVec;

  // Vectors containing the parameters
  //std::vector<updParam> globalParams;
  //std::vector<updParam> localParams;

  std::vector<updParam> params;
  
  // TODO: For speed, return const?
  int numPars() const {
    return params.size(); //globalParams.size() + localParams.size();
  }

  updParam& operator[](paramIndex index) {
    return params[index];
  }
  const updParam& operator[](paramIndex index) const {
    return params[index];
  }
  updParam& operator[](int index) {
    // TODO BOUNDS CHECK
    return params[(paramIndex) index];
  }
  const updParam& operator[](int index) const {
    return params[(paramIndex) index];
  }
  
  // Contain the parameter values for each block
  updParamBlock globalParBlock;
  std::vector<updParamBlock> localParBlock;

  // Initialise: Called once to set up blocks
  void init(int regions);

  // Insert a new parameter (copy from old-style parameter)
  // regional: true = local par, false = global
  void insertAndRegister(updateable_model_parameter& inPar, int num_instances, bool regional, int numRegions);

  // Lookup parameter by param name from enum
  //const updParam& lookup(paramIndex index) const;

  // Perform M-H step
  void updateParams(gsl_rng *rng);

  // Return a gsl_vector for the given parameter, i.e. a sub-vector
  // containing only this parameter's values
  // WARNING: This is a pointer to a gsl_vector_const_view, so the pointer
  // does not own the data and should not modify it in any way; it should not
  // be passed 'upwards' to calling functions as the view may be stored on the stack
  const gsl_vector_const_view lookupValue(paramIndex index, int region = 0) const;
  const gsl_vector_const_view lookupValue(int index, int region = 0) const {
    return lookupValue((paramIndex) index, region);
  }

  // Return 1st element of value
  const double lookupValue0(paramIndex index, int region = 0) const;
  

  infection_to_data_delay gp_delay;
  infection_to_data_delay hosp_delay;
  infection_to_data_delay death_delay;
  
  // - - - - - - - - - - - - - - - - -
  
  // Future work: The mod_pars.txt file should be read from here. The existing
  // input code is difficult to read and edit. But it works, so leave it for now.
  // At present, these two methods are unused.
  // Current defaults in RTM_StructAssign
  updParam defaultParam(std::string &name);
  
  // Takes a parameter name and returns a parameter containing the relevant 
  // defaults. Setting defaults like this makes them easier to read and modify.
  void readInputFile();

  void calcProposals(gsl_rng *rng);

  void calcAccept();
  void doAccept(gsl_rng *rng);
  void adaptiveUpdate(int iter);


};


// Equivalent of a typedef to shorten name
// Does this simplify or increase confusion?
using upd = updParamSet;


// Class contains a single model parameter (upd = updateable)
class updParam {
public:

  std::string param_name;
  upd::paramIndex index;	// Index of param in names enum
  unsigned int size;    // Number of components in param
  int regionSize;    // For local params only, size of each region (= size/numRegions)
  
  gslVector init_value;   // Initial value of parameter. Used only in setup.
  unsigned int value_index;  // Starting index in either global or local value vector
  bool local;    // True if local param (differs over region); false if global

  // Prior params: Vector of length number of components.
  // Each element is the vector containing the prior for that component, or
  // if size 0, the prior is the parameters of the parameter in parents[i].
  std::vector<int> parents;
  std::vector<gslVector> prior_params;
  gslVectorInt prior_distribution; // Flag giving distribution of individual components
  bool flag_hyperprior;
  double log_prior_dens;
  double proposal_log_prior_dens;
  bool flag_update; // True if prior_distribution[i] != cDETERM for all i

  // gslVector proposal_variances;
  mvnorm *prior_multivariate_norm;
  regression_def map_to_regional;
  gslVector posterior_mean;
  gslVector posterior_sumsq;
  bool flag_transmission_model; // True if number of new infected needs to be recalculated
  bool flag_reporting_model; // If reporting delays need to be recalculated
  bool flag_GP_likelihood; // GP consultation likelihood needs to be recalculated
  bool flag_Hosp_likelihood; // Hospitalisation likelihood
  bool flag_Sero_likelihood; // Seroepidemiology data likelihood
  bool flag_Viro_likelihood; // Virology (positivity) data recalculated
  bool flag_Prev_likelihood; // Likelihood for prevalence data need to be recalculated
  bool flag_any_child_nodes; // True if any flag child nodes are true
  std::vector<bool> flag_child_nodes; // Flag for each of the other parameters of the mode, indicating whether they are child nodes of the current node

  // Constructor for specifying default values
  // Flag is a binary 7-bit value consisting of the following 7 flags:
  // transmission, reporting, gp, hosp, sero, viro, prev data likelihood
  // Not currently used.
  updParam(string name, double value, int flags);
  
  // Constructor for current init by copying from an existing param
  // TODO: Remove this when input code improved
  updParam(updateable_model_parameter &in, const int num_instances, bool regional);

  // Need default constructor for constructing params vec in updParamSet
  updParam() { }

  
};





#endif // HEADER_updParams_
