#ifndef HEADER_updParams_
#define HEADER_updParams_

#include <cassert>
#include <exception>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include <gsl/gsl_rng.h>

#include "gslWrapper.h"
#include "RTM_StructDefs.h"

// upd = Updateable
class updParam;
class updParamBlock;

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

  // The parameters
  std::vector<updParam> params;

  // Vector of parameter blocks. First block (index 0) is global.
  std::vector<updParamBlock> blocks;  

  // Mapping of string names to paramIndex above
  static const std::map<std::string, paramIndex> nameMap;

  infection_to_data_delay gp_delay;
  infection_to_data_delay hosp_delay;
  infection_to_data_delay death_delay;

  // Useful to have
  int numRegions;
  
  updParamSet() :
    params(UMP_FINAL_NUM_PARAMS) { }

  int numPars() const {
    return UMP_FINAL_NUM_PARAMS;
  }
  
  // Initialise blocks. Called once after all params created.
  void init(int numRegions_, const string& outdir);

  // Insert a new parameter (copy from old-style parameter)
  void insertAndRegister(updateable_model_parameter& inPar, int num_instances, bool local, int numRegions);

  // Return a view for the given parameter, i.e. a sub-vector
  // Warning: View is only valid as long as parameter exists, and is not moved etc
  const gsl_vector_const_view lookupValue(paramIndex index, int region = 0) const;
  const gsl_vector_const_view lookupValue(int index, int region = 0) const;

  // Return 1st element of value
  double lookupValue0(paramIndex index, int region = 0) const;

  // - - - - - - -
  // Key M-H code
  // Calculate new proposal vectors (sequential)
  void calcProposals(gsl_rng *rng);
  // Calculate acceptance ratios (in parallel over blocks)
  void calcAccept(Region* country, const global_model_instance_parameters& gmip, const mixing_model& base_mix);
  // Implement acceptance
  void doAccept(gsl_rng *rng, Region* country, const global_model_instance_parameters& gmip);
  // Adapt MH distribution over time
  void adaptiveUpdate(int iter);
  // DEBUG output
  void outputPars();

  void printAcceptRates(int numIters);

  // - - - - - - -
  // Future work: Initialise pars directly from input file rather than copying
  // from old-style param. Also, store defaults in this class.
  updParam defaultParam(std::string &name);
  void readInputFile();

  // Helper:
  // Use [] as shortcut to get individual param
  updParam& operator[](paramIndex index) {
    return params[index];
  }
  const updParam& operator[](paramIndex index) const {
    return params[index];
  }
  updParam& operator[](int index) {
    assert(index >= 0 && index < UMP_FINAL_NUM_PARAMS);
    return params[(paramIndex) index];
  }
  const updParam& operator[](int index) const {
    assert(index >= 0 && index < UMP_FINAL_NUM_PARAMS);
    return params[(paramIndex) index];
  }
};

// Short typedef for specifying paramIndex constants
using upd = updParamSet;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

// Block of parameter values
// Either all global param values, or local values for single region
class updParamBlock {
public:
  bool global;	// True for the global block
  gslVector vals;
  gslVector proposal;
  gslVector dist; // Distribution for each of the vals; copy of vals in updParam
  gslVector mu;
  double beta;
  gslMatrix sigma;
  double laccept;
  int acceptLastMove;  // Int for use in adaptive update
  int numAccept;
  // Region index corresponding to this block, zero-based
  // Also zero in the global block (for use in lookup index maths)
  int regionNum;
  
  likelihood lfx;
  likelihood prop_lfx;

  double childProposalDensity;
  double childCurrentDensity;
  
  // TODO: Breaks default copy/assign/destruct
  Region* propCountry;
  void copyCountry(Region* inCountry, int numRegions);
  void setLlhood(likelihood l) {
    lfx = l;
    // Do this to set sizes
    prop_lfx = l;
  }
  
  // Match the key M-H methods from updParamSet
  void calcProposal(gsl_rng *rng);
  void calcAccept(updParamSet &paramSet, Region* country, const global_model_instance_parameters& gmip, const mixing_model& base_mix);
  void doAccept(gsl_rng *rng, updParamSet& paramSet, Region* country, int numRegions, const global_model_instance_parameters& gmip);
  void adaptiveUpdate(int iter);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

// Class for a single model parameter
class updParam {
public:

  std::string param_name;
  upd::paramIndex index;    // Index of param in names enum
  unsigned int size;	    // Number of components in param
  int regionSize;	    // For local params, num components per region
  
  gslVector init_value;     // Initial value of parameter. Used only in setup.
  unsigned int value_index; // Starting index in either global or local value vector
  bool global;		    // True if global param, false if local

  // Prior params: Vector of length number of components.
  std::vector<gslVector> prior_params;
  // If prior is separate parameter, prior_params[i].size() == 0 and the index of
  // the prior is stored in parents[i].
  std::vector<int> parents;
  std::vector<distribution_type> prior_distribution; // Flag giving distribution of individual components
  bool flag_hyperprior;
  double log_prior_dens;
  double proposal_log_prior_dens;
  bool flag_update; // True if prior_distribution[i] != cDETERM for all i

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

  std::ofstream outfile;

  // Default constructor for constructing params vec in updParamSet
  updParam() { }

  // Constructor for current init by copying from an existing param
  // TODO: Remove this when input code improved
  updParam(updateable_model_parameter &in, const int num_instances, bool regional);

  // Constructor for specifying default values
  // Flag is a binary 7-bit value consisting of the following 7 flags:
  // transmission, reporting, gp, hosp, sero, viro, prev data likelihood
  // Not currently used.
  updParam(string name, double value, int flags);
};


#endif // HEADER_updParams_
