#ifndef HEADER_updParams_
#define HEADER_updParams_

#include <cassert>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include <gsl/gsl_rng.h>

#include "gslWrapper.h"
#include "RTM_StructDefs.h"

// Debugging flags
//#define USE_OLD_CODE
//#define SKIP_NEW_CODE
extern bool debug;  // Set in RTM_updParams.cc


// upd = Updateable
class updParam;
class updParamBlock;

// The new block structure needs an array of design matrices for each param
class updRegrDef {
public:
  gslVectorInt region_breakpoints;
  gslVectorInt age_breakpoints;
  gslVectorInt time_breakpoints;
  link_function regression_link;
  std::vector<gslMatrix> design_matrix;

  updRegrDef(regression_def& in)
    : region_breakpoints(in.region_breakpoints),
      age_breakpoints(in.age_breakpoints),
      time_breakpoints(in.time_breakpoints),
      regression_link(in.regression_link) { }

  updRegrDef() { }
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

  // The parameters
  std::vector<updParam> params;

  // Parameter component blocks. First block (index 0) is global, others local.
  std::vector<updParamBlock> blocks;  

  // Mapping of string names to paramIndex above
  static const std::map<std::string, paramIndex> nameMap;

  infection_to_data_delay gp_delay;
  infection_to_data_delay hosp_delay;
  infection_to_data_delay death_delay;

  int numRegions;
  
  updParamSet() :
    params(UMP_FINAL_NUM_PARAMS) { }

  int numPars() const {
    assert(params.size() == UMP_FINAL_NUM_PARAMS);
    return UMP_FINAL_NUM_PARAMS;
  }
  
  // Insert a new parameter (copy from old-style parameter)
  void insertAndRegister(updateable_model_parameter& inPar, int num_instances, bool local, int numRegions_);

  // Initialise blocks. Called once after all params created.
  void init(const string& outdir);

  // Return a view for the given parameter, i.e. a sub-vector
  // Warning: View is only valid as long as parameter exists, and is not moved etc
  const gsl_vector_const_view lookup(paramIndex index, int region) const;
  const gsl_vector_const_view lookup(int index, int region) const;

  // Return 1st element of value
  double lookup0(paramIndex index, int region) const;

  // - - - - - - -
  // Key M-H code
  // Calculate new proposal vectors (sequential)
  void calcProposals(updParamSet& paramSet, gsl_rng *rng);
  // Calculate acceptance ratios (in parallel over blocks)
  void calcAccept(Region* country, const global_model_instance_parameters& gmip, const mixing_model& base_mix);
  // Accept or reject
  void doAccept(gsl_rng *rng, Region* country, const global_model_instance_parameters& gmip);
  // Adapt MH distribution over time
  void adaptiveUpdate(int iter);
  
  // Debugging output
  void outputPars();
  void outputProposals();
  void printAcceptRates(int numIters);

  // - - - - - - -
  // Future work: Initialise pars directly from input file rather than copying
  // from old-style param. Also, store defaults in this class.
  updParam defaultParam(std::string &name);
  void readInputFile();

  // Helper:
  // Enable use of [] as shortcut to get individual param
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

  gslVector values;     // Values: Copied from parameters at each update
  gslVector parIndex;   // Index of param each component comes from
  gslVector parOffset;  // Index of component in that param
  gslVector proposal;   // Proposal value
  std::vector<distribution_type> dist;       // Prior distribution for each component
			// Copied from param object, but fixed during execution

  bool global;	// True for the global block
  // Zero-based region number for block. Also zero in global block.
  int regionNum;
 
  // Proposal distribution
  gslVector mu;
  gslMatrix sigma;

  // For adapting the proposal distribution
  // Static vars are a hack for the values shared among all local blocks
  double beta;
  static double regbeta;
  int acceptLastMove;
  static int regacceptLastMove;  

  double laccept;   // Acceptance score
  int numAccept;    // Acceptance rate
  int numProposed;
  
  likelihood lfx;
  likelihood prop_lfx;

  // Initialise likelihood
  void setLlhood(likelihood& l) {
    lfx = l;
    prop_lfx = l;  // Sets vector sizes of prop_lfx
  }
  double childProposalDensity;
  double childCurrentDensity;
  
  // TODO: Region ptr breaks default copy/assign/destruct
  // Safe for now as blocks aren't moved/copied once added to blocks vector
  Region* propCountry;
  void copyCountry(Region* inCountry, int numRegions);

  // M-H methods
  void calcProposal(updParamSet& paramSet, gsl_rng *rng);
  void calcAccept(updParamSet &paramSet, Region* country, const global_model_instance_parameters& gmip, const mixing_model& base_mix);
  void doAccept(gsl_rng *rng, updParamSet& paramSet, Region* country, int numRegions, const global_model_instance_parameters& gmip);
  void adaptiveUpdate(int iter);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

// Class for a single model parameter
class updParam {
public:

  std::string param_name;
  gslVector values;
  
  upd::paramIndex index;  // Index of param in names enum
  unsigned int size;  // Number of components in param
  int regionSize;     // For local params, num components per region
	              // Warning: *Not* number of components per proposal block
  bool global;	// True if global param, false if local

  // Prior params: Vector of length number of components.
  std::vector<gslVector> prior_params;
  // If prior is separate parameter, prior_params[i].size() == 0 and the index of
  // the prior is stored in parents[i].
  std::vector<int> parents;
  // Distribution for each indvidiual component
  std::vector<distribution_type> prior_distribution;
  bool flag_hyperprior;
  double log_prior_dens;
  double proposal_log_prior_dens;
  bool flag_update; // True if prior_distribution[i] != cCONSTANT for all i

  mvnorm *prior_multivariate_norm;
  updRegrDef map_to_regional;
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
