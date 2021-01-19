#include "RTM_updParams.h"
#include "RTM_FunctDefs.h"

#include <fstream>

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


void updParamSet::init(int numRegions_) {

  numRegions = numRegions_;
  
  int globalSize = 0, localSize = 0;
  for (auto &par : params) {
    if (par.global)
      globalSize += par.size;
    else
      localSize += par.size;
  }
  
  // Split local params by region. Integer division
  localSize /= numRegions;

  blocks.resize(numRegions + 1);
  
  blocks[0].vals.alloc(globalSize);
  blocks[0].regionNum = -1;
  for (int i = 1; i < blocks.size(); i++) {
    blocks[i].vals.alloc(localSize);
    blocks[i].regionNum = i-1;
  }

  
  int globalIndex = 0;
  std::vector<int> localIndex(numRegions, 0);
  for (auto &par : params) {
    par.value_index = localIndex[0];
    if (par.global) {
      par.value_index = globalIndex;
      for (int i = 0; i < par.size; i++)
	blocks[0].vals[globalIndex++] = par.init_value[i];
    } else {
      // local block
      int readIndex = 0;
      int regionSize = par.size / numRegions; // Integer division
      for (int r = 0; r < numRegions; r++) {
	for (int i = 0; i < regionSize; i++)
	  blocks[r+1].vals[localIndex[r]++] = par.init_value[readIndex++];
      }
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
    
  for (auto &block : blocks) { 
    // Init the paramBlocks
    block.mu = block.vals;

    // Set the diagonal of sigma matrix
    block.sigma.allocZero(block.vals.size(), block.vals.size());
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
  if (par.global) {
    return gsl_vector_const_subvector(*blocks[0].vals, par.value_index, par.size);
  } else {
    // local
    //gsl_vector_const_view view = ...
    return gsl_vector_const_subvector(*blocks[region+1].vals, par.value_index, par.regionSize);
    //return &(view.vector);
  }
}

const gsl_vector_const_view updParamSet::lookupValue(int index, int region) const {
  // TODO: CHECK REGION BOUNDS
  return lookupValue((paramIndex) index, region);
}


const double updParamSet::lookupValue0(upd::paramIndex index, int region) const {
  const updParam &par = params[index];
  int i = par.value_index;
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
  gslMatrix covar = sigma * exp(beta);

  gsl_ran_multivariate_gaussian(rng, *vals, *covar, *proposal);

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
    // Either both flags are global or both are local
    if (global && par.global) {

      // par.proposal_log_prior_dens = 0;

      // Iterate over components in block
      // For each region, we only consider that region's components

      flagclass update_flags(par.index);
      
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
	    int distIndex = parIndex;
	    if (par.flag_hyperprior)
	      parIndex = 1;
	    
	    laccept += R_univariate_prior_log_density_ratio(
	      proposal[blockIndex],
	      vals[blockIndex],
	      par.prior_distribution[parIndex],
	      *par.prior_params[parIndex]);
	  }
	}
      }

      if (laccept == GSL_NEGINF)
	return;

      if (par.prior_distribution[0] == cMVNORMAL) {
	gsl_vector_view prop = gsl_vector_subvector(*proposal, par.value_index, par.regionSize);
	gsl_vector_view value = gsl_vector_subvector(*vals, par.value_index, par.regionSize);
	
	laccept += par.prior_multivariate_norm->ld_mvnorm_ratio(
	  &prop.vector, &value.vector);
	
	if (laccept == GSL_NEGINF)
	  return;
	
	if (par.flag_any_child_nodes) {
	  // Loop over other parameters
	  for (int ch = 0; ch < par.flag_child_nodes.size(); ch++) {
	    if (par.flag_child_nodes[ch]) {

	      updParam &child = paramSet[ch];
	      
	      cout << "Parent-child: parent global " << par.global << ", child global " << child.global << endl;
	      
	      for (int i = 0; i < child.size; i++) {
		// univariate_prior_log_dens:
		int childValueIndex = child.value_index + i;

		assert(prop.vector.size == 2);
		
		laccept += R_univariate_prior_log_density(
		  vals[childValueIndex],
		  child.prior_distribution[0],
		  // parameter is theta_i->proposal_value, the parent proposal
		  // It is assumed this is 2 params long
		  &prop.vector);
	      }
	    } else {
	      // No child nodes

	      // TODO TODO TODO - COPY REGIONS
	      
	      if (par.global) {
		// Evaluate all regions

		// TODO: Save param_value and set param_value = proposal_value

		for (int reg = 0; reg < paramSet.numRegions; reg++) {
		  block_regional_parameters(propCountry[reg].det_model_params, paramSet, gmip, reg, country[reg].population, propCountry[reg].total_population, base_mix, update_flags);
		}
	      } else {
		// Local. Evaluate only region of this block
		int reg = regionNum;
		block_regional_parameters(propCountry[reg].det_model_params, paramSet, gmip, reg, country[reg].population, propCountry[reg].total_population, base_mix, update_flags);
	      }
	      // TODO: Reset param_value
	      
	      block_log_likelihood(prop_lfx, propCountry, 0,
				   par.flag_transmission_model,
				   par.flag_reporting_model,
				   par.flag_GP_likelihood,
				   par.flag_Hosp_likelihood,
				   par.flag_Viro_likelihood,
				   par.flag_Sero_likelihood,
				   par.flag_Prev_likelihood,
				   gmip,
				   paramSet);

	      laccept += prop_lfx.total_lfx - lfx.total_lfx;
	      
	    }
	  }
	}
      }
    }
  }
}

void updParamSet::doAccept(gsl_rng *rng, Region* country) {
  for (auto & block : blocks)
    block.doAccept(rng, country, numRegions);
}

void updParamBlock::doAccept(gsl_rng *rng, Region* country, int numRegions) {
  double acceptTest = gsl_sf_log(gsl_rng_uniform(rng));

  // TODO: Initialise properly??
  flagclass update_flags;

  
  if (laccept > acceptTest) {
    numAccept++;
    acceptLastMove = 1;
    vals = proposal;

    // TODO: Save log prior dens?
    // if (par.flag_any_child_nodes)
    //   save child log prior dens


    lfx = prop_lfx;

    for (int r = 0; r < numRegions; r++)
      regional_model_params_memcpy(country[r].det_model_params, propCountry[r].det_model_params, update_flags);

    // TODO: model_statistics_memcpy
    
  } else {
    // Undo necessary stuff
    for (int r = 0; r < numRegions; r++)
      regional_model_params_memcpy(propCountry[r].det_model_params, country[r].det_model_params, update_flags);

    // TODO: model_statistics_memcpy

    // TODO:
    // if (par.flag_any_child_nodes)
    prop_lfx = lfx;
  }
}


// TODO: Hardcoded dir name
// TODO: If this is kept in production, open files at startup and close on exit
void updParamSet::outputPars() {
  for (auto par : params) {
    string filename = "pars/" + par.param_name + ".txt";
    std::ofstream outfile(filename, std::ofstream::app);
    if (par.global) {
      for (int i = 0; i < par.size; i++)
	outfile << blocks[0].vals[par.value_index + i] << " ";
      outfile << std::endl;
    } else {
      // local
      for (int r = 0; r < numRegions; r++)
	for (int i = 0; i < par.regionSize; i++)
	  outfile << blocks[r+1].vals[par.value_index + i] << " ";
      outfile << std::endl;
    }
    outfile.close();
  }
}
