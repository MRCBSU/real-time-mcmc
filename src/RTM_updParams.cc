#include "RTM_updParams.h"
#include "RTM_FunctDefs.h"


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
    init_value(value),
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
    init_value (in.param_value),
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
  
  // Prior distribution is a vector, not gsl_vector_int*, so need copy and cast
  prior_distribution.resize(in.prior_distribution->size);
  for (int i = 0; i < prior_distribution.size(); i++)
    prior_distribution[i] = (distribution_type) gsl_vector_int_get(in.prior_distribution, i);
  
  if (flag_update) {
    // prior_params - copy from gsl_vector**
    for (int i = 0; i < prior_params.size(); i++)
      prior_params[i] = gslVector(in.prior_params[i]);
  }
  // flag_child_nodes - copy from bool*
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
void updParamSet::insertAndRegister(updateable_model_parameter& inPar, int num_instances, bool regional, int numRegions) {
  
  updParam newPar(inPar, num_instances, regional);
  newPar.index = nameMap.at(newPar.param_name);

  if (regional) {
    if (newPar.size % numRegions != 0) {
      std::cerr << "Error: parameter " << newPar.param_name << " is marked as regional but number of regions does not evenly divide number of components\n";
      std::cerr << "Param size: " << newPar.size << " , numRegions: " << numRegions << endl;
      exit(1);
    }

    newPar.regionSize = newPar.size / numRegions;     // Integer division
  }

  // For local parameters, we need to resize design matrix
  if (inPar.map_to_regional.design_matrix.nrows() > 0) {
    if (! newPar.global) {
      // Take submatrices for each param
      // Assume main design matrix is square
      assert(inPar.map_to_regional.design_matrix.nrows() == inPar.map_to_regional.design_matrix.ncols());
      
      int size = inPar.map_to_regional.design_matrix.nrows() / numRegions;
      newPar.map_to_regional.design_matrix.eraseRealloc(size, size);
      gsl_matrix_view submat = gsl_matrix_submatrix(*inPar.map_to_regional.design_matrix, 0, 0, size, size);
      newPar.map_to_regional.design_matrix = &submat.matrix;
    }
  }
  
  params[newPar.index] = std::move(newPar);
}

// low/high only used for uniform
double transform(double x, int dist, double low = 0, double high = 1) {
  switch(dist) {
  case cCONSTANT:
    return x;
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

// newparam, trunc_flag, low, high
double invTransform(double x, int dist, double low = 0, double high = 1) {
  switch(dist) {
  case cCONSTANT:
    return x;
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
    throw invalid_argument("transform: invalid distribution");
  }
}

void updParamSet::init(int numRegions_, const string& dir) {

  numRegions = numRegions_;

  // Open output files
  // Will write over existing contents
  for (auto& par : params) {
    if (! par.flag_update)
      continue;
    string filename = dir + "/" + par.param_name + ".txt";
    par.outfile.open(filename, std::ofstream::out);
    if (! par.outfile.is_open()) {
      std::cerr << "Error: Unable to open output file " << filename << "\n";
      exit(1);
    }
  }
  
  int globalSize = 0, localSize = 0;
  for (auto &par : params) {
    if (par.flag_update)
      if (par.global)
	globalSize += par.size;
      else
	localSize += par.size;
  }
  
  // Split local params by region. Integer division
  localSize /= numRegions;

  blocks.resize(numRegions + 1);
  
  blocks[0].vals.alloc(globalSize);
  blocks[0].dist.alloc(globalSize);
  blocks[0].regionNum = 0;
  blocks[0].global = true;

  for (int i = 1; i < blocks.size(); i++) {
    blocks[i].vals.alloc(localSize);
    blocks[i].dist.alloc(localSize);
    blocks[i].regionNum = i-1;
    blocks[i].global = false;
  }

  
  int globalIndex = 0;
  std::vector<int> localIndex(numRegions, 0);
  for (auto &par : params) {
    if (par.flag_update) {
      par.value_index = localIndex[0];
      if (par.global) {
	par.value_index = globalIndex;
	for (int i = 0; i < par.size; i++) {
	  blocks[0].vals[globalIndex] = par.init_value[i];
	  blocks[0].dist[globalIndex] = par.prior_distribution[i];
	  globalIndex++;
	}
      } else {
	// local block
	int readIndex = 0;
	int regionSize = par.size / numRegions; // Integer division
	for (int r = 0; r < numRegions; r++) {
	  for (int i = 0; i < regionSize; i++) {
	    blocks[r+1].vals[localIndex[r]] = par.init_value[readIndex];
	    blocks[r+1].dist[localIndex[r]] = par.prior_distribution[readIndex];
	    localIndex[r]++;
	    readIndex++;
	  }
	}
      }
    }	
  }
  
  // Initialise mu and sigma for each block.
  // Need to transform params accordingly
      
  for (auto &block : blocks) { 
    // Init the paramBlocks
    block.mu.alloc(block.vals.size());
    for (int i = 0; i < block.vals.size(); i++)
      block.mu[i] = transform(block.vals[i], block.dist[i]);

    // Set the diagonal of sigma matrix
    block.sigma.allocZero(block.vals.size(), block.vals.size());
    for (int i = 0; i < block.vals.size(); i++) {
      block.sigma[i][i] = fabs(block.mu[i]) * 0.01;
      if (block.sigma[i][i] == 0)
	block.sigma[i][i] = 0.01;
    }
    
    block.beta = 0;
    block.proposal.alloc(block.vals.size());
  }
}


const gsl_vector_const_view updParamSet::lookupValue(upd::paramIndex index, int region) const {
  const updParam &par = params[index];
  if (! par.flag_update) {
    // Constant param, return initial value
    return gsl_vector_const_subvector(*par.init_value, 0, par.size);
  }
  if (par.global) {
    return gsl_vector_const_subvector(*blocks[0].vals, par.value_index, par.size);
  } else {
    // local
    return gsl_vector_const_subvector(*blocks[region+1].vals, par.value_index, par.regionSize);
  }
}

const gsl_vector_const_view updParamSet::lookupValue(int index, int region) const {
  // TODO: CHECK REGION BOUNDS
  return lookupValue((paramIndex) index, region);
}


double updParamSet::lookupValue0(upd::paramIndex index, int region) const {
  const updParam &par = params[index];
  int i = par.value_index;
  if (! par.flag_update)
    return par.init_value[0];
  if (par.global) {
    return blocks[0].vals[i];
  } else {
    return blocks[region+1].vals[i];
  }
}


void updParamSet::calcProposals(gsl_rng *rng) {  
  // #pragma omp parallel for
  for (auto &block : blocks)
    block.calcProposal(rng);
}

// Main block update method
void updParamBlock::calcProposal(gsl_rng *rng) {

  // Calculate proposal
  
  // Multivar normal. mu = param values, not the variable 'mu'.
  // After 1st 200 iters, sigma/beta change every iter, so no benefit to caching
  gslMatrix covar = sigma * exp(beta);

  gslVector transformed(vals.size());
  for (int i = 0; i < vals.size(); i++)
    transformed[i] = transform(vals[i], dist[i]);
  
  gsl_linalg_cholesky_decomp1(*covar);
  gsl_ran_multivariate_gaussian(rng, *transformed, *covar, *proposal);

  for (int i = 0; i < vals.size(); i++) {
    if (dist[i] == cCONSTANT)
      proposal[i] = vals[i];
    else
      proposal[i] = invTransform(proposal[i], dist[i]);
  }
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
  for (auto &block : blocks)
    block.adaptiveUpdate(iter);
}

void updParamSet::calcAccept(Region* country, const global_model_instance_parameters& gmip, const mixing_model& base_mix) {
  // #pragma omp parallel for
  for (auto &block : blocks)
    block.calcAccept(*this, country, gmip, base_mix);
}

void updParamBlock::copyCountry(Region* inCountry, int numRegions) {
  propCountry = new Region[numRegions];
  flagclass all_pos;
  for (int r = 0; r < numRegions; r++) {
    Region_alloc(propCountry[r], inCountry[r]);
    Region_memcpy(propCountry[r], inCountry[r], all_pos);
  }
}

void updParamBlock::calcAccept(updParamSet &paramSet, Region* country, const global_model_instance_parameters& gmip, const mixing_model& base_mix) {

  laccept = 0;
  acceptLastMove = 0;
  
  for (auto &par : paramSet.params) {

    // Skip over constant params
    if (! par.flag_update)
      continue;
    
    // Both block and par are either both global or both local, so consider par
    if (global == par.global) {

      // par.proposal_log_prior_dens = 0;

      // Iterate over components in block
      // For each region, we only consider that region's components

      flagclass update_flags(par.index);

      // Iterate over components in the current block
      // (for global blocks, regionSize is the same as the parameter size)
      for (int i = 0; i < par.regionSize; i++) {
	// value = block.vals[par.value_index + i]
	// Region component index in par: regionNum * regionSize + i;

	// Index of component in parameter and block
	int parIndex = regionNum * par.regionSize + i;
	int blockIndex = par.value_index + i;

	if (par.prior_distribution[parIndex] != cCONSTANT) {
	  if (par.prior_distribution[parIndex] != cMVNORMAL) {
	  
	    // Code save proposal_log_prior_dens at this point
	    // laccept += univariate_prior_ratio(par, a_component = component index)

	    // dist_component = flag_hyperprior ? 1 : component
	    const gsl_vector* prior;
	    int distIndex;
	    if (par.flag_hyperprior) {
	      assert(par.prior_params[parIndex].size() == 0);
	      // Use parent distribution
	      gsl_vector_const_view parentView = paramSet.lookupValue(par.parents[parIndex]);
	      prior = &parentView.vector;
	      distIndex = 1;
	    } else {
	      prior = *par.prior_params[parIndex];
	      distIndex = parIndex;
	    }
	      
	    laccept += R_univariate_prior_log_density_ratio(
	      proposal[blockIndex],
	      vals[blockIndex],
	      par.prior_distribution[distIndex],
	      prior);
	    
	  }
	}
      }

      if (laccept == GSL_NEGINF)
	return;

      if (par.prior_distribution[0] == cMVNORMAL) {

	// TODO: Test this code path
	cout << "WARNING: Multivariate normal code path untested\n";
	gsl_vector_view prop = gsl_vector_subvector(*proposal, par.value_index, par.regionSize);
	gsl_vector_view value = gsl_vector_subvector(*vals, par.value_index, par.regionSize);
	
	laccept += par.prior_multivariate_norm->ld_mvnorm_ratio(
	  &prop.vector, &value.vector);
      }
	
      if (laccept == GSL_NEGINF)
	return;
      
      if (par.flag_any_child_nodes) {
	// Loop over other parameters
	
	for (int ch = 0; ch < par.flag_child_nodes.size(); ch++) {
	  if (par.flag_child_nodes[ch]) {
	    
	    updParam &child = paramSet[ch];
	    	    
	    if (child.global) {
	      // TODO TODO TODO
	      // Iterate over child and add up R_univariate_prior_log_density for each component
	      // For now, we have no global child nodes
	    } else {

	      child.proposal_log_prior_dens = 0;
		
	      // Child components are split over regions.
	      for (int r = 0; r < paramSet.numRegions; r++) {
		for (int i = 0; i < par.regionSize; i++) {

		  double value = paramSet.blocks[r+1].vals[i];
		  gsl_vector_view prop = gsl_vector_subvector(*proposal, par.value_index, par.regionSize);

		  // We assume parent is 2 params, as it is acting as a prior
		  assert(par.size == 2);
		  
		  child.proposal_log_prior_dens += R_univariate_prior_log_density(
		    value,
		    child.prior_distribution[0],
		    // parameter is theta_i->proposal_value, the parent proposal
		    &prop.vector);
		}
	      }

	      laccept += child.proposal_log_prior_dens - child.log_prior_dens;
	    }
	  }
	}
      } else {
	// No child nodes
	
	gslVector originalVals = vals;
	vals = proposal;
	
	if (par.global) {
	  // Evaluate all regions
	  for (int reg = 0; reg < paramSet.numRegions; reg++) {
	    block_regional_parameters(propCountry[reg].det_model_params, paramSet, gmip, reg, propCountry[reg].population, propCountry[reg].total_population, base_mix, update_flags);
	  }
	} else {
	  // Local. Evaluate only region of this block
	  int reg = regionNum;
	  block_regional_parameters(propCountry[reg].det_model_params, paramSet, gmip, reg, propCountry[reg].population, propCountry[reg].total_population, base_mix, update_flags);
	}
	
	vals = originalVals;
	
	fn_log_likelihood(prop_lfx, propCountry, 0,
			  par.flag_transmission_model,
			  par.flag_reporting_model,
			  par.flag_GP_likelihood,
			  par.flag_Hosp_likelihood,
			  par.flag_Viro_likelihood,
			  par.flag_Sero_likelihood,
			  par.flag_Prev_likelihood,
			  gmip,
			  paramSet.gp_delay.distribution_function,
			  paramSet.hosp_delay.distribution_function
	  );
	
	laccept += prop_lfx.total_lfx - lfx.total_lfx;
      }
    }
  }
}

void updParamSet::doAccept(gsl_rng *rng, Region* country, const global_model_instance_parameters& gmip) {
  for (auto & block : blocks)
    block.doAccept(rng, *this, country, numRegions, gmip);
}

void updParamBlock::doAccept(gsl_rng *rng, updParamSet& paramSet, Region* country, int numRegions, const global_model_instance_parameters& gmip) {
  double acceptTest = gsl_sf_log(gsl_rng_uniform(rng));

  // TODO: Initialise properly??
  flagclass update_flags;
    
  if (laccept > acceptTest) {
    //cout << "accept\n";
    
    numAccept++;
    acceptLastMove = 1;
    vals = proposal;

    // Log prior dens: Original code saves this, but appears not to use the
    // saved values anywhere. Except for child nodes:
    
    // TODO: Speed this up?
    for (auto& par : paramSet.params) {
      if (par.flag_any_child_nodes) {
	for (int i = 0; i < par.size; i++)
	  if (par.flag_child_nodes[i]) {
	    updParam& child = paramSet[i];
	    child.log_prior_dens = child.proposal_log_prior_dens;
	  }
      }
    }
    
    lfx = prop_lfx;

    if (global) {
      for (int r = 0; r < numRegions; r++) {
	regional_model_params_memcpy(country[r].det_model_params, propCountry[r].det_model_params, update_flags);

	// TODO: Existing code does this parameter by parameter, so passes a bunch of parameter-specific flags
	// For now, set all these flags to true. (Only downside is a slower copy.)
	model_statistics_memcpy(country[r].region_modstats, propCountry[r].region_modstats,
				true,
				(bool) gmip.l_GP_consultation_flag,
			        (bool) gmip.l_Hospitalisation_flag,
				true,
				(bool) gmip.l_Viro_data_flag,
				(bool) gmip.l_Prev_data_flag);

      }
      
    } else {
      // local
      int r = regionNum;
      
      regional_model_params_memcpy(country[r].det_model_params, propCountry[r].det_model_params, update_flags);
      model_statistics_memcpy(country[r].region_modstats, propCountry[r].region_modstats,
			      true,
			      (bool) gmip.l_GP_consultation_flag,
			      (bool) gmip.l_Hospitalisation_flag,
			      true,
			      (bool) gmip.l_Viro_data_flag,
			      (bool) gmip.l_Prev_data_flag);
    }
  } else {
    //cout << "reject\n";
    
    // Undo region changes stuff
    if (global) {
      for (int r = 0; r < numRegions; r++) {
	regional_model_params_memcpy(propCountry[r].det_model_params, country[r].det_model_params, update_flags);
	
	model_statistics_memcpy(propCountry[r].region_modstats, country[r].region_modstats,
				true,
				(bool) gmip.l_GP_consultation_flag,
				(bool) gmip.l_Hospitalisation_flag,
				true,
				(bool) gmip.l_Viro_data_flag,
				(bool) gmip.l_Prev_data_flag);
      }
    } else {
      // local
      int r = regionNum;
      
      regional_model_params_memcpy(propCountry[r].det_model_params, country[r].det_model_params, update_flags);
      
      model_statistics_memcpy(propCountry[r].region_modstats, country[r].region_modstats,
			      true,
			      (bool) gmip.l_GP_consultation_flag,
			      (bool) gmip.l_Hospitalisation_flag,
			      true,
			      (bool) gmip.l_Viro_data_flag,
			      (bool) gmip.l_Prev_data_flag);

    }
    // TODO:
    // if (!par.flag_any_child_nodes)
    prop_lfx = lfx;
  }
  
}


void updParamSet::outputPars() {
  for (auto& par : params) {
    if (! par.flag_update)
      continue;
    if (par.global) {
      for (int i = 0; i < par.size; i++)
	par.outfile << blocks[0].vals[par.value_index + i] << " ";
      par.outfile << std::endl;
    } else {
      // local
      for (int r = 0; r < numRegions; r++)
	for (int i = 0; i < par.regionSize; i++)
	  par.outfile << blocks[r+1].vals[par.value_index + i] << " ";
      par.outfile << std::endl;
    }
  }
}

void updParamSet::printAcceptRates(int numIters) {
  cout << "Accept: ";
  for (auto& block : blocks)
    cout << setprecision(2) << block.numAccept / (double) numIters << " ";
  cout << endl;
}
