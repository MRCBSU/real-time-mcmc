#include "RTM_updParams.h"
#include "RTM_FunctDefs.h"


bool debug = false;


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


// Initialise from defaults in updParamSet. Method currently unused.
updParam::updParam(string name, double value, int flags)
  : param_name(name),
    values(value),
    // Note: Specifying bitmasks with 0b is a C++14 feature, but is a compiler
    // extension in both gcc (4.3 or later) and clang (since at least 2013).
    flag_transmission_model (flags & 0b1000000),
    flag_reporting_model    (flags & 0b0100000),
    flag_GP_likelihood      (flags & 0b0010000),
    flag_Hosp_likelihood    (flags & 0b0001000),
    flag_Sero_likelihood    (flags & 0b0000100),
    flag_Viro_likelihood    (flags & 0b0000010),
    flag_Prev_likelihood    (flags & 0b0000001)
{ }


updParam::updParam(updateable_model_parameter &in, const int num_instances, bool regional)
  : param_name (in.param_name),
    index(upd::INVALID), // -ve indicates invalid; set this during 'insert()'
    size (in.param_value->size),
    regionSize(in.param_value->size),  // For local pars this will be reset
    values (in.param_value),
    // value_index: set in paramSet initialisation
    global(!regional),
    prior_params (in.param_value->size),    
    parents (in.param_value->size, upd::INVALID),
    
    //prior_distribution: Change of type
    flag_hyperprior (in.flag_hyperprior),
    log_prior_dens (in.log_prior_dens),
    proposal_log_prior_dens (in.proposal_log_prior_dens), 
    flag_update (in.flag_update),

    prior_multivariate_norm(0),
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
  
  // Convert from gsl_vector_int* to std::vector<distribution_type>
  prior_distribution.resize(in.prior_distribution->size);
  for (int i = 0; i < prior_distribution.size(); i++)
    prior_distribution[i] = (distribution_type) gsl_vector_int_get(in.prior_distribution, i);
  
  if (flag_update) {
    // Copy from gsl_vector**
    for (int i = 0; i < prior_params.size(); i++)
      prior_params[i] = gslVector(in.prior_params[i]);
  }
  
  // Copy from bool*
  for (int i = 0; i < num_instances; i++)
    flag_child_nodes[i] = in.flag_child_nodes[i];
}


// Defaults. Incomplete
updParam updParamSet::defaultParam(std::string &name) {
  if (name == "exponential_growth_rate_hyper")
    return updParam(name, 0.15, 0b0000000);
  else if (name == "l_p_lambda_0_hyper")
    return updParam(name, -15.0, 0b0000000);
  // else if (...)
  else
    throw std::invalid_argument("");
}


// Insert parameter into parameter set
void updParamSet::insertAndRegister(updateable_model_parameter& inPar, int num_instances, bool regional, int numRegions_) {

  numRegions = numRegions_;
  
  updParam newPar(inPar, num_instances, regional);
  newPar.index = nameMap.at(newPar.param_name);
  newPar.map_to_regional.design_matrix.resize(numRegions);

  if (! regional) {
    // Global par. Make n copies of the design matrix (each region has its own copy)
    for (int r = 0; r < numRegions; r++)
      newPar.map_to_regional.design_matrix[r] = inPar.map_to_regional.design_matrix;
    
  } else {
    // Local
    if (newPar.size % numRegions != 0) {
      std::cerr << "Error: Number of components of parameter " << newPar.param_name << " is not evenly divided by number of regions\n";
      exit(1);
    }
    
    newPar.regionSize = newPar.size / numRegions;     // Integer division

    // Number of columns = num param components
    // Number of rows = total of (region, time, age) group combinations,
    //   OR num param components, whichever is larger
    int nrows = inPar.map_to_regional.design_matrix.nrows() / numRegions;
    int ncols = newPar.regionSize;

    for (int r = 0; r < numRegions; r++) {
      newPar.map_to_regional.design_matrix[r].eraseRealloc(nrows, ncols);
      gsl_matrix_view submat = gsl_matrix_submatrix(*inPar.map_to_regional.design_matrix, r*nrows, r*ncols, nrows, ncols);
      newPar.map_to_regional.design_matrix[r] = &submat.matrix;
    }
  }

  params[newPar.index] = std::move(newPar);
}

// low/high only used for uniform dist
double transform(double x, distribution_type dist, double low = 0, double high = 1) {
  
  switch(dist) {
  case cCONSTANT:
  case cNORMAL:
  case cMVNORMAL:
    return x;
  case cGAMMA:
  case cHALFNORMAL:
    // double min is ~= 2e-308
    return gsl_sf_log(x - std::numeric_limits<double>::min());
  case cBETA:
    return gsl_sf_log(x) - gsl_sf_log(1-x);
  case cUNIFORM:
    return gsl_sf_log(x - low) - gsl_sf_log(high - x);
  default:
    throw invalid_argument("transform: invalid distribution");
  }
}

double invTransform(double x, distribution_type dist, double low = 0, double high = 1) {

  switch(dist) {
  case cCONSTANT:
  case cNORMAL:
  case cMVNORMAL:
    return x;
  case cGAMMA:
  case cHALFNORMAL:
    return std::numeric_limits<double>::min() + gsl_sf_exp(x);
  case cBETA:
    return gsl_sf_exp(x) / (1 + gsl_sf_exp(x));
  case cUNIFORM:
    return (high * gsl_sf_exp(x) + low) / (1 + gsl_sf_exp(x));
  default:
    throw invalid_argument("invTransform: invalid distribution");
  }
}

double jacobianFactor(double prop, double curr, distribution_type dist, double low = 0, double high = 1) {
  
  switch(dist) {
  case cCONSTANT:
  case cNORMAL:
  case cMVNORMAL:
    return 0;
  case cGAMMA:
  case cHALFNORMAL:
    return gsl_sf_log(prop/curr);
  case cBETA:
    // a=0, b=1
    return gsl_sf_log(prop/curr) + gsl_sf_log((1-prop)/(1-curr));
  case cUNIFORM:
    return gsl_sf_log(((prop-low) * (high-prop)) / ((curr-low) * (high-curr)));
  default:
    throw invalid_argument("jacobianFactor: invalid distribution");
  }
}

void updParamSet::init(const string& dir) {

  // Debug output: parameter output files. Overwrites existing files
  // Requires that 'dir' already exists. Default param value is "pars"
  if (debug) {
    for (updParam& par : params) {
      if (par.flag_update) {
	string filename = dir + "/" + par.param_name + ".txt";
	par.outfile.open(filename, std::ofstream::out);
	if (! par.outfile.is_open()) {
	  std::cerr << "Error: Unable to open output file " << filename << "\n";
	  exit(1);
	}
      }
    }
  }

  // Count non-constant components
  int globalSize = 0, localSize = 0;
  for (updParam &par : params) {
    if (par.flag_update)
      for (int i = 0; i < par.size; i++)
	if (par.prior_distribution[i] != cCONSTANT) {
	  if (par.global)
	    globalSize++;
	  else
	    localSize++;
	}
  }
  localSize /= numRegions;

  // Alloc space
  blocks.resize(numRegions + 1);
  
  blocks[0].values.alloc(globalSize);
  blocks[0].dist.resize(globalSize);
  blocks[0].parIndex.alloc(globalSize);
  blocks[0].parOffset.alloc(globalSize);
  blocks[0].regionNum = 0;
  blocks[0].global = true;

  for (int i = 1; i < blocks.size(); i++) {
    blocks[i].values.alloc(localSize);
    blocks[i].dist.resize(localSize);
    blocks[i].parIndex.alloc(localSize);
    blocks[i].parOffset.alloc(localSize);
    blocks[i].regionNum = i-1;
    blocks[i].global = false;
  }

  // Insert components to blocks
  // Counters give current insertion index for each region
  int globalIndex = 0;
  std::vector<int> localIndex(numRegions, 0);
  for (updParam &par : params) {
    if (par.flag_update) {
      if (par.global) {
	for (int i = 0; i < par.size; i++) {
	  if (par.prior_distribution[i] != cCONSTANT) {
	    blocks[0].values[globalIndex] = par.values[i];
	    blocks[0].parIndex[globalIndex] = par.index;
	    blocks[0].parOffset[globalIndex] = i;
	    blocks[0].dist[globalIndex] = par.prior_distribution[i];
	    globalIndex++;
	  }
	}
      } else {
	// local block

	// Count number of non-const components
	int nonconst = 0;
	for (int i = 0; i < par.size; i++)
	  if (par.prior_distribution[i] != cCONSTANT)
	    nonconst++;

	if (nonconst % numRegions != 0) {
	  cerr << "Error: Number of non-constant components of param " << par.param_name << " is not evenly divisible by number of regions\n";
	  exit(1);
	}
	int perRegion = nonconst / numRegions; // Integer division
	
	int r = 0;   // Region number
	int rcount = 0; // Count components per region
	for (int i = 0; i < par.size; i++) {
	  if (par.prior_distribution[i] != cCONSTANT) {
	    // localIndex holds insertion points for each block
	    blocks[r+1].values[localIndex[r]] = par.values[i];
	    blocks[r+1].parIndex[localIndex[r]] = par.index;
	    blocks[r+1].parOffset[localIndex[r]] = i;
	    blocks[r+1].dist[localIndex[r]] = par.prior_distribution[i];

	    localIndex[r]++;
	    rcount++;
	    if (rcount >= perRegion) {
	      // We have filled this region
	      r++;
	      rcount = 0;
	    }
	  }
	}
      }
    }	
  }

  // Initialise mu and sigma for each block. Use transformed param values
  for (updParamBlock &block : blocks) { 

    block.proposal.alloc(block.values.size());

    block.mu.alloc(block.values.size());
    for (int i = 0; i < block.values.size(); i++)
      block.mu[i] = transform(block.values[i], block.dist[i]);

    // Set the diagonal of sigma matrix
    block.sigma.allocZero(block.values.size(), block.values.size());
    for (int i = 0; i < block.values.size(); i++) {
      block.sigma[i][i] = fabs(block.mu[i]);
      if (block.sigma[i][i] == 0)
	block.sigma[i][i] = 1;
    }

    block.sigma *= 0.01;   // Scale
    block.beta = 0;
  }
}


// Region is a number between 0 and r-1. This is ignored if the parameter is global.
const gsl_vector_const_view updParamSet::lookup(upd::paramIndex index, int region) const {
  const updParam& par = params[index];
  if (par.global) {
    // Return the whole vector
    return gsl_vector_const_subvector(*par.values, 0, par.size);
  } else {
    int start = par.regionSize * region;
    return gsl_vector_const_subvector(*par.values, start, par.regionSize);
  }
}

const gsl_vector_const_view updParamSet::lookup(int index, int region) const {
  return lookup((paramIndex) index, region);
}

double updParamSet::lookup0(upd::paramIndex index, int region) const {
  const updParam &par = params[index];
  if (par.global)
    return par.values[0];
  else {
    int start = par.regionSize * region;
    return par.values[start];
  }
}

void updParamBlock::copyCountry(Region* inCountry, int numRegions) {
  propCountry = new Region[numRegions];
  flagclass all_pos;
  for (int r = 0; r < numRegions; r++) {
    Region_alloc(propCountry[r], inCountry[r]);
    Region_memcpy(propCountry[r], inCountry[r], all_pos);
  }
}
void updParamSet::calcProposals(updParamSet& paramSet, gsl_rng *rng) {  
  for (auto &block : blocks)
    block.calcProposal(paramSet, rng);
}

void updParamBlock::calcProposal(updParamSet& paramSet, gsl_rng *rng) {

  // Calculate proposal
  // Assumes that the values vector is already up to date.
  // It is set on init, and set whenever proposal is accepted
  
  // Covariance matrix. Sigma changes every iter, so no benefit to caching
  gslMatrix covar;
  covar = sigma * exp(beta);

  // Transform
  gslVector transformed(values.size());
  for (int i = 0; i < values.size(); i++)
    transformed[i] = transform(values[i], dist[i]);
  
  // Random sample
  gsl_linalg_cholesky_decomp1(*covar);
  gsl_ran_multivariate_gaussian(rng, *transformed, *covar, *proposal);
  
  for (int i = 0; i < values.size(); i++)
    proposal[i] = invTransform(proposal[i], dist[i]);
}

double updParamBlock::regbeta = 0;
int updParamBlock::regacceptLastMove = 0;


void updParamSet::calcAccept(Region* country, const global_model_instance_parameters& gmip, const mixing_model& base_mix) {
  // #pragma omp parallel for
  for (auto &block : blocks)
    block.calcAccept(*this, country, gmip, base_mix);
}


void updParamBlock::calcAccept(updParamSet& paramSet, Region* country, const global_model_instance_parameters& gmip, const mixing_model& base_mix) {

  laccept = 0;
  if (global)
    acceptLastMove = 0;
  else
    regacceptLastMove = 0;

  // Get region
  // TODO: Can we only copy one region for local blocks?
  flagclass update_flags;
  for (int r = 0; r < paramSet.numRegions; r++) {
    regional_model_params_memcpy(propCountry[r].det_model_params, country[r].det_model_params, update_flags);
  }
  
  // Prior densities for components in the block
  for (int i = 0; i < values.size(); i++) {
    laccept += jacobianFactor(proposal[i], values[i], dist[i]);

    if (dist[i] != cMVNORMAL) {

      // Get prior from node
      // If the prior is a parent node, get the parent value.
      // (Original code used a pointer here so no need for extra logic)
      updParam& par = paramSet[parIndex[i]];
      gsl_vector* prior = *(par.prior_params[parOffset[i]]);
      if (prior->size == 0) {
	// Use parent
	int parent = par.parents[parOffset[i]];
	prior = *paramSet[parent].values;
      }
      
      laccept += R_univariate_prior_log_density_ratio(
	proposal[i], values[i], dist[i], prior);
      
    } else {
      // TODO: Fix
      cout << "ERROR: Multivariate normal code path not implemented\n";
      exit(1);
    }
  }
  
  if (laccept == GSL_NEGINF)
    return;

  for (updParam& par : paramSet.params) {

    // Skip over constant params
    if (! par.flag_update)
      continue;

    // If the flags match, the parameter is part of this block
    if (global == par.global) {
      
      if (par.flag_any_child_nodes) {
	// Find the child node(s)
	for (int ch = 0; ch < par.flag_child_nodes.size(); ch++) {
	  if (par.flag_child_nodes[ch]) {
	    
	    updParam &child = paramSet[ch];
	    	    
	    if (child.global) {
	      // TODO TODO TODO
	      // Sum R_univariate_prior_log_density for each component
	      std::cerr << "calcAccept(): child node is global; not implemented yet\n";
	      exit(1);
	    } else {

	      child.proposal_log_prior_dens = 0;
		
	      // We assume parent is 2 params, as it is acting as a prior
	      assert(par.size == 2);

	      // We haven't accepted, so the prior is the relevant proposal components

	      // Difference in child density between proposal and current value
	      // Caching no longer easy to do, and is inefficient if updating all
	      // regional blocks per iter
	      
	      // WARNING WARNING TOFIX
	      // Proper way of setting the prior.
	      // Mix of const and non-const, so can't just subset the proposal vector

	      // Current value
	      gslVector prevprior(2);
	      prevprior[0] = par.values[0];
	      prevprior[1] = values[0];
	      
	      child.log_prior_dens = 0;
	      for (int i = 0; i < child.size; i++) {
		child.log_prior_dens += R_univariate_prior_log_density(
		    child.values[i],
		    child.prior_distribution[i],
		    *prevprior);
	      }

	      // Proposed value
	      gslVector prior(2);
	      prior[0] = par.values[0];
	      prior[1] = proposal[0];

	      for (int i = 0; i < child.size; i++) {
		child.proposal_log_prior_dens += R_univariate_prior_log_density(
		    child.values[i],
		    child.prior_distribution[i],
		    *prior);
	      }

	      laccept += child.proposal_log_prior_dens - child.log_prior_dens;
	    }
	  }
	}
      }
    }
  }

  // Do the region update and likelihood for the block as a whole, rather than
  // param by param

  // TEMPORARY HACK
  // Global parameters: hosp_negbin_overdispersion, infectious_period, prop_case_to_hosp, sero_test_sensitivity, sero_test_specificity
  // Correspond to regional params: l_hosp_negbin_overdispersion, l_average_infectious_period, l_R0_init, l_I0, l_R0_Amplitude, l_pr_onset_to_Hosp, l_sero_sensitivity, l_sero_specificity
  
  // Local: contact_parameters, log_beta_rw, exponential_growth_rate, log_p_lambda_0
  // Corresponding regional: l_MIXMOD, l_lbeta_rw, l_EGR, l_R0_init, l_I0, l_R0_Amplitude

  flagclass blockflags;
  // Default is all flags set to true
  for (int i = 0; i < 27; i++)
    blockflags.regional_update_flags[i] = false;
  
  if (global) {
    blockflags.regional_update_flags[2] = true;
    blockflags.regional_update_flags[8] = true;
    blockflags.regional_update_flags[10] = true;
    blockflags.regional_update_flags[11] = true;
    blockflags.regional_update_flags[14] = true;
    blockflags.regional_update_flags[23] = true;
    blockflags.regional_update_flags[25] = true;
    blockflags.regional_update_flags[26] = true;
  } else {
    blockflags.regional_update_flags[6] = true;
    blockflags.regional_update_flags[7] = true;
    blockflags.regional_update_flags[8] = true;
    blockflags.regional_update_flags[10] = true;
    blockflags.regional_update_flags[11] = true;
    blockflags.regional_update_flags[18] = true;
  }
  
  // Copy proposal to paramSet to evaluate region
  for (int i = 0; i < values.size(); i++)
    paramSet[parIndex[i]].values[parOffset[i]] = proposal[i];
  
  if (global) {
    // Evaluate all regions
    for (int reg = 0; reg < paramSet.numRegions; reg++) {
      block_regional_parameters(propCountry[reg].det_model_params, paramSet, gmip, reg, propCountry[reg].population, propCountry[reg].total_population, base_mix, blockflags);
      
    }
  } else {
    // Local. Evaluate only region of this block
    int reg = regionNum;
    block_regional_parameters(propCountry[reg].det_model_params, paramSet, gmip, reg, propCountry[reg].population, propCountry[reg].total_population, base_mix, blockflags);
  }
  
  // Need to restore original values so that other blocks can run
  for (int i = 0; i < values.size(); i++)
    paramSet[parIndex[i]].values[parOffset[i]] = values[i];

  // For at least one of the updated parameters, all of the flags are set to true
  fn_log_likelihood(prop_lfx, propCountry, 0,
		    true, true, true, true, true, true, true, 
		    gmip,
		    paramSet.gp_delay.distribution_function,
		    paramSet.hosp_delay.distribution_function
    );
	
  laccept += prop_lfx.total_lfx - paramSet.lfx.total_lfx;
}

void updParamSet::doAccept(gsl_rng *rng, Region* country, const global_model_instance_parameters& gmip) {
  for (auto & block : blocks)
    block.doAccept(rng, *this, country, numRegions, gmip);
}

void updParamBlock::doAccept(gsl_rng *rng, updParamSet& paramSet, Region* country, int numRegions, const global_model_instance_parameters& gmip) {
  
  double acceptTest = gsl_sf_log(gsl_rng_uniform(rng));

  flagclass update_flags; // Set flags to update whole region

  numProposed++;

  if (laccept > acceptTest) {
    numAccept++;
   
    if (global)
      acceptLastMove = 1;
    else
      regacceptLastMove = 1;

    // Update the parameters
    for (int i = 0; i < values.size(); i++)
      paramSet[parIndex[i]].values[parOffset[i]] = proposal[i];
    
    // To calculate the adaptive delta, we swap these two, so proposal
    // temporarily holds the previous iteration's values.
    values.swapWith(proposal);

    // TODO: Can only global nodes have children?
    // TODO: Value no longer cached. Decide whether caching or recalculating is most efficient
    /*    
    if (global) {
      // Save child log prior density
      for (updParam& par : paramSet.params)
      if (par.flag_any_child_nodes)
	for (int i = 0; i < par.flag_child_nodes.size(); i++)
	  if (par.flag_child_nodes[i]) {
	    updParam& child = paramSet[i];
	    child.log_prior_dens = child.proposal_log_prior_dens;
	  }
    }
    */
    
    paramSet.lfx = prop_lfx;

    if (global) {
      for (int r = 0; r < numRegions; r++) {
	regional_model_params_memcpy(country[r].det_model_params, propCountry[r].det_model_params, update_flags);

	model_statistics_memcpy(
	  country[r].region_modstats, propCountry[r].region_modstats,
	  true, gmip.l_GP_consultation_flag, gmip.l_Hospitalisation_flag,
	  true, gmip.l_Viro_data_flag, gmip.l_Prev_data_flag);
      }
    } else {
      // local
      int r = regionNum;
      
      regional_model_params_memcpy(country[r].det_model_params, propCountry[r].det_model_params, update_flags);
      model_statistics_memcpy(
	country[r].region_modstats, propCountry[r].region_modstats,
	true, gmip.l_GP_consultation_flag, gmip.l_Hospitalisation_flag,
	true, gmip.l_Viro_data_flag, gmip.l_Prev_data_flag);
    }
  } else {
    // Reject proposal

    // We need to reset the proposal so that adpative delta is calculated correctly 
    proposal = values;
    
    // Undo region changes 
    if (global) {
      for (int r = 0; r < numRegions; r++) {
	// Don't need to restore region - we do this at the start of the next iter
	// regional_model_params_memcpy(propCountry[r].det_model_params, country[r].det_model_params, update_flags);
	
	model_statistics_memcpy(
	  propCountry[r].region_modstats, country[r].region_modstats,
	  true,  gmip.l_GP_consultation_flag, gmip.l_Hospitalisation_flag,
	  true, gmip.l_Viro_data_flag, gmip.l_Prev_data_flag);
      }
    } else {
      // local
      int r = regionNum;
      model_statistics_memcpy(
	propCountry[r].region_modstats, country[r].region_modstats,
	true, gmip.l_GP_consultation_flag, gmip.l_Hospitalisation_flag,
	true, gmip.l_Viro_data_flag, gmip.l_Prev_data_flag);
    }

    // TODO do we need to do this?
    prop_lfx = paramSet.lfx;
  }
}



void updParamBlock::adaptiveUpdate(int iter) {

  gslVector delta(values.size());
  for (int i = 0; i < values.size(); i++) {
    delta[i] = transform(values[i], dist[i]) - mu[i];
  }

  double betastart = beta;
  
  double eta = pow(iter - 199 + 1, -0.6);
  if (global)
    beta = beta + eta * (acceptLastMove - 0.234);
  else
    beta = beta + eta * (updParamBlock::regacceptLastMove - 0.234);
  
  for (int i = 0; i < mu.size(); i++)
    mu[i] = (1 - eta) * mu[i] + eta * transform(values[i], dist[i]);
  
  gslMatrix deltaProd(delta.size(), delta.size());
    
  for (int i = 0; i < delta.size(); i++)
    for (int j = 0; j < delta.size(); j++)
      deltaProd[i][j] = delta[i] * delta[j];
    
  sigma = (1 - eta) * sigma + eta * deltaProd;
}


void updParamSet::adaptiveUpdate(int iter) {
  for (updParamBlock& block : blocks)
    block.adaptiveUpdate(iter);
}


void updParamSet::outputPars() {
  for (auto& par : params) {
    if (par.flag_update) {
      par.outfile << "  ";
      for (int i = 0; i < par.size; i++)
	par.outfile << par.values[i] << " ";
      par.outfile << std::endl;
    }
  }
}

void updParamSet::outputProposals() {
  cerr << "WARNING: outputProposals() currently disabled\n";
  for (auto& par : params) {
    if (! par.flag_update)
      continue;
    par.outfile << "P ";
    if (par.global) {
      //for (int i = 0; i < par.size; i++)
      //par.outfile << blocks[0].proposal[par.value_index + i] << " ";
      //par.outfile << std::endl;
    } else {
      // local
      //for (int r = 0; r < numRegions; r++)
      //for (int i = 0; i < par.regionSize; i++)
      //par.outfile << blocks[r+1].proposal[par.value_index + i] << " ";
      //par.outfile << std::endl;
    }
  }
}


void updParamSet::printAcceptRates(int numIters) {
  cout << "Acceptance Rates: ";
  for (auto& block : blocks)
    cout << setprecision(4) << block.numAccept / (double) block.numProposed << " ";
  cout << endl;
}
