#include "RTM_updParams.h"
#include <cassert>


// Note: Specifying bitmasks with 0b is a C++14 feature, but is a compiler
// extension in both gcc (4.3 or later) and clang (since at least 2013).
// TODO: Confirm whether default should be single element
// TODO: Resizing when reading from input file?
updParam::updParam(string name, double value, int flags)
  : param_name(name),
    init_value(value),
    flag_transmission_model (flags & 0b1000000),
    flag_reporting_model    (flags & 0b0100000),
    flag_GP_likelihood      (flags & 0b0010000),
    flag_Hosp_likelihood    (flags & 0b0001000),
    flag_Sero_likelihood    (flags & 0b0000100),
    flag_Viro_likelihood    (flags & 0b0000010),
    flag_Prev_likelihood    (flags & 0b0000001)
{
}


updParam::updParam(updateable_model_parameter &in, const int num_instances, bool regional)
  : param_name (in.param_name),
    index(upd::INVALID), // -ve indicates invalid; set this during 'insert()'
    init_value (in.param_value),
    size (in.param_value->size),
    regionSize(-1),
    local(regional),
    //proposal_value (in.proposal_value),
    prior_distribution (in.prior_distribution),
    flag_hyperprior (in.flag_hyperprior),
    log_prior_dens (in.log_prior_dens),
    proposal_log_prior_dens (in.proposal_log_prior_dens), 
    flag_update (in.flag_update),
    // Length is the same as length of param_value
    // -1 indicates 'no prior'
    parents (in.param_value->size, upd::INVALID),
    prior_params (in.param_value->size),
    //proposal_variances (in.proposal_variances),
    prior_multivariate_norm(0),
    //map_to_regional(),
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




const std::map<std::string, upd::paramIndex> updParamSet::nameMap = {
  { "exponential_growth_rate_hyper", EGR_HYPER },  // 0
  { "l_p_lambda_0_hyper", LPL0_HYPER },
  { "prop_susceptible_hyper", PROP_SUS_HYPER },
  { "log_beta_rw_sd", LBETA_RW_SD  },
  { "gp_negbin_overdispersion", GP_OVERDISP }, 
  { "hosp_negbin_overdispersion", HOSP_OVERDISP }, // 5
  { "latent_period", ALP },
  { "infectious_period", AIP },
  { "r1_period", AR1 },
  { "relative_infectiousness", REL_INFECT },
  { "prop_symptomatic", PROP_SYMP },  // 10
  { "contact_parameters", CONTACT },
  { "log_beta_rw", LBETA_RW },
  { "R0_amplitude_kA", R0_AMP },
  { "R0_seasonal_peakday", R0_PEAKDAY },
  { "exponential_growth_rate", EGR },  // 15
  { "log_p_lambda_0", LPL0 },
  { "prop_susceptible", PROP_SUS },
  { "prop_HI_32_to_HI_8", PROP_HI_GEQ_32 },
  { "prop_case_to_GP_consultation", PROP_GP },
  { "prop_case_to_hosp", PROP_HOSP },  // 20
  { "prop_case_to_death", PROP_DEATH },
  { "importation_rates", IMPORTATION },
  { "background_GP", BGR },
  { "test_sensitivity", SENS },
  { "test_specificity", SPEC },  // 25
  { "day_of_week_effects", DOW_EFFECTS },
  { "sero_test_sensitivity", SSENS },
  { "sero_test_specificity", SSPEC }
};



void updParamSet::insertAndRegister(updateable_model_parameter& inPar, int num_instances, bool regional, int numRegions) {
  
  updParam newPar(inPar, num_instances, regional);
  newPar.index = nameMap.at(newPar.param_name);

  if (regional) {

    if (newPar.size % numRegions != 0) {
      std::cerr << "Error: parameter " << newPar.param_name << " is marked as regional but number of regions does not evenly divide number of components\n";
      std::cerr << "Param size: " << newPar.size << " , numRegions: " << numRegions << endl;
      exit(1);
    }

    // Integer division
    newPar.regionSize = newPar.size / numRegions;
  }
  
  params[newPar.index] = std::move(newPar);

  // We don't have deep copy/move for the regression_def struct, so we do this
  // after putting the param in the vector.
  
  updParam &par = params[newPar.index];

  // Create a new regression_def object
  par.map_to_regional = regression_def();
  
  if (inPar.map_to_regional.region_breakpoints != 0) {
    par.map_to_regional.region_breakpoints = gsl_vector_int_alloc(inPar.map_to_regional.region_breakpoints->size);
    gsl_vector_int_memcpy(par.map_to_regional.region_breakpoints, inPar.map_to_regional.region_breakpoints);
  } else {
    par.map_to_regional.region_breakpoints = 0;
  }

  if (inPar.map_to_regional.age_breakpoints != 0) {
    par.map_to_regional.age_breakpoints = gsl_vector_int_alloc(inPar.map_to_regional.age_breakpoints->size);
    gsl_vector_int_memcpy(par.map_to_regional.age_breakpoints, inPar.map_to_regional.age_breakpoints);
  } else {
    par.map_to_regional.age_breakpoints = 0;
  }

  if (par.map_to_regional.time_breakpoints != 0) {
    par.map_to_regional.time_breakpoints = gsl_vector_int_alloc(inPar.map_to_regional.time_breakpoints->size);
    gsl_vector_int_memcpy(par.map_to_regional.time_breakpoints, inPar.map_to_regional.time_breakpoints);
  } else {
    par.map_to_regional.time_breakpoints = 0;
  }
  
  par.map_to_regional.regression_link = inPar.map_to_regional.regression_link;

  // WARNING WARNING
  // INCOMPLETE: Should have vector of design matrices for regions

  if (inPar.map_to_regional.design_matrix != 0) {
    if (par.local) {
      // Take submatrices for each param
      // Assume main design matrix is square
      assert(inPar.map_to_regional.design_matrix->size1 == inPar.map_to_regional.design_matrix->size2);
      
      int size = inPar.map_to_regional.design_matrix->size1 / numRegions;
      par.map_to_regional.design_matrix = gsl_matrix_alloc(inPar.map_to_regional.design_matrix->size1, size);
      gsl_matrix_view submat = gsl_matrix_submatrix(inPar.map_to_regional.design_matrix, 0, 0, inPar.map_to_regional.design_matrix->size1, size);
      gsl_matrix_memcpy(par.map_to_regional.design_matrix, &submat.matrix);
      
    } else {
      // Global. Simply copy.
      par.map_to_regional.design_matrix = gsl_matrix_alloc(inPar.map_to_regional.design_matrix->size1, inPar.map_to_regional.design_matrix->size2);
      gsl_matrix_memcpy(par.map_to_regional.design_matrix, inPar.map_to_regional.design_matrix);
    }
  } else {
    inPar.map_to_regional.design_matrix = 0;
  }
  
    
  /*
  int index;
  
  if (regional) {
    localParams.push_back(std::move(newPar));
    index = localParams.size() * -1;
    
  } else {
    globalParams.push_back(std::move(newPar));
    index = globalParams.size();  
  }

  int i = nameMap.at(inPar.param_name);
  lookupVec[i] = index;
  */
  
  }


updParam updParamSet::defaultParam(std::string &name) {
    if (name == "exponential_growth_rate_hyper")
      return updParam(name, 0.15, 0b0000000);
    else if (name == "l_p_lambda_0_hyper")
      return updParam(name, -15.0, 0b0000000);
    // else if (...)
    else
      throw std::invalid_argument("");
    }

void updParamSet::init(int numRegions) {

  int globalSize = 0, localSize = 0;
  for (auto &par : params) {
    if (par.local)
      localSize += par.size;
    else
      globalSize += par.size;
  }
  
  // Split local params by region. Integer division
  localSize /= numRegions;
  
  globalParBlock.vals.alloc(globalSize);
  localParBlock.resize(numRegions);
  for (auto &block : localParBlock)
    block.vals.alloc(localSize);
  
  int globalIndex = 0;
  std::vector<int> localIndex(numRegions, 0);
  for (auto &par : params) {
    par.value_index = localIndex[0];
    if (par.local) {
      int readIndex = 0;
      int regionSize = par.size / numRegions; // Integer division
      for (int r = 0; r < numRegions; r++) {
	for (int i = 0; i < regionSize; i++)
	  localParBlock[r].vals[localIndex[r]++] = par.init_value[readIndex++];
      }
      
    } else { // par.global)
      par.value_index = globalIndex;
      for (int i = 0; i < par.size; i++)
	globalParBlock.vals[globalIndex++] = par.init_value[i];
    }
  }	

      

    /*
  
    // Total number of local params
    int count = 0;
    for (auto &par : localParams)
      count += par.param_value.size();

    // TODO: Error check the two left devision calculations

    // Integer division
    count = count / regions;

    // Init the vectors for each region
    localParBlock.resize(regions);
    for (auto &block : localParBlock)
      block.vals.alloc(count);

    // Copy

    // Insertion indices for each of the regional vectors
    std::vector<int> index(regions, 0);
    for (auto &par : localParams) {
      // Number of elts for each region
      // Integer division
      int regionSize = par.param_value.size() / regions;

      // Counter for reading from param_value
      int readIndex = 0;
      for (int r = 0; r < regions; r++) {
	for (int i = 0; i < regionSize; i++)
	  localParBlock[r].vals[index[r]++] = par.param_value[readIndex++];
      }
    }
    
    // Global params
    count = 0;
    for (auto &par : globalParams)
      count += par.param_value.size();
    
    globalParBlock.vals.alloc(count);
    
    int outInd = 0;
    for (auto &par : globalParams)
      for (int i = 0; i < par.param_value.size(); i++)
	globalParBlock.vals[outInd++] = par.param_value[i];
*/
  
    // Init the paramBlocks
  globalParBlock.mu = globalParBlock.vals;
  
  // Set the diagonal of sigma matrix
  globalParBlock.sigma.alloc(globalParBlock.vals.size(), globalParBlock.vals.size());
  for (int i = 0; i < globalParBlock.vals.size(); i++)
    globalParBlock.sigma[i][i] = fabs(globalParBlock.mu[i]) * 0.01;
  
  globalParBlock.beta = 0;
  globalParBlock.proposal.alloc(globalParBlock.vals.size());
  
  for (auto &block : localParBlock) { 
    block.mu = block.vals;
    
    block.sigma.alloc(block.vals.size(), block.vals.size());
    for (int i = 0; i < block.vals.size(); i++)
      block.sigma[i][i] = fabs(block.mu[i]) * 0.01;
    
    block.beta = 0;
    block.proposal.alloc(block.vals.size());
  }
    
  // TODO: Get initial likelihoods for each block
  
}


/*
const updParam& updParamSet::lookup(paramIndex index) const {
  int i = lookupVec[index];
  if (i > 0)
    return globalParams[i-1];
  else if (i < 0)
    return localParams[-1 * (i-1)];
  else // i == 0
    return nullptr;
    //throw std::invalid_argument("lookup index is not a valid parameter index");
}
*/

//const gsl_vector* updParamSet::lookupValue(paramIndex index, int region) const {
const gsl_vector_const_view updParamSet::lookupValue(upd::paramIndex index, int region) const {
//  int pind = lookupVec[index];
  const updParam &par = params[index];
  if (par.local) {
    //gsl_vector_const_view view = ...
    return gsl_vector_const_subvector(*localParBlock[region].vals, par.value_index, par.regionSize);
    //return &(view.vector);
  } else {
    return gsl_vector_const_subvector(*globalParBlock.vals, par.value_index, par.size);
    //return &(view.vector);
  }
}

const double updParamSet::lookupValue0(upd::paramIndex index, int region) const {
  const updParam &par = params[index];
  if (par.local) {
    return localParBlock[region].vals[0];
  } else {
    return globalParBlock.vals[0];
  }
}



void updParamSet::calcProposals(gsl_rng *rng) {
  globalParBlock.calcProposal(rng);
  
  // #pragma omp parallel for
  for (auto &block : localParBlock)
    block.calcProposal(rng);

}

// Main block update method

void updParamBlock::calcProposal(gsl_rng *rng) {

  // Calculate proposal
  
  // Multivar normal. mu = param values, not the variable 'mu'.
  gslMatrix covar = sigma * exp(beta);
  
  gsl_ran_multivariate_gaussian(rng, *vals, *covar, *proposal);

  std::cout << "vals: " << vals << std::endl;
  std::cout << "prop: " << proposal << std::endl;
  
}

void updParamBlock::adaptiveUpdate(int iter) {
  double eta = pow(iter - 199 + 1, -0.6);
  beta = beta + eta * (acceptLastMove - 0.234);

  // delta = proposed vals - curr vals
  gslVector delta = proposal - vals;
  mu = (1 - eta) * mu + eta * delta;
  gslMatrix deltaProd(delta.size(), delta.size());

  for (int i = 0; i < delta.size(); i++)
    for (int j = 0; j < delta.size(); j++)
      deltaProd[i][j] = delta[i] * delta[j];
  
  sigma = (1 - eta) * sigma + eta * deltaProd;
}


void updParamSet::adaptiveUpdate(int iter) {
  globalParBlock.adaptiveUpdate(iter);

  for (auto &block : localParBlock)
    block.adaptiveUpdate(iter);
}

void updParamSet::calcAccept() {
  // Parallelise
  globalParBlock.calcAccept(*this, false);

  // #pragma omp parallel for
  for (auto &block : localParBlock)
    block.calcAccept(*this, true);
}

void updParamBlock::calcAccept(updParamSet &paramSet, bool local) {

  laccept = 0;
  acceptLastMove = 0;
  
  for (auto &par : paramSet.params) {
    if (local && par.local) {
      // Param counts for this block
      //laccept += prior_density;
    }
  }
}

void updParamSet::doAccept(gsl_rng *rng) {
  globalParBlock.doAccept(rng);
  for (auto & block : localParBlock)
    block.doAccept(rng);
}

void updParamBlock::doAccept(gsl_rng *rng) {
  double acceptTest = gsl_sf_log(gsl_rng_uniform(rng));
  if (laccept > acceptTest) {
    acceptLastMove = 1;
  } else {
    // Undo necessary stuff
  }
}
