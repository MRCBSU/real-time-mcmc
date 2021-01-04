#include "RTM_BlockUpdate.h"

void updateGlobalParams(const gslVector &globalParamValues, const std::vector<updateableParam> &globalParams) {

    // class param_block {
    //   gslVector paramValues
    //   vector beta
    //   matrix sigma
    //   log_posterior
    // }
    
    // Proposal
    
    // theta_star = sample_normal(globalParamValues, exp(beta) * sigma);

    // Evaluate log posterior of theta_star
    // ?? ld_mvnorm_ratio
    
    // accept_test = log_post(theta_star) - log_post(theta)
    //
    // # pragma omp barrier
    // log post accesses all params across all regions; cannot
    // update until all accept_tests are complete
    // 
    // if (accept_test) > threshold
    //   globalParamValues = theta_star
    //   save log_posterior
    //

    // if (t > 199)
    //   Calculate nu and update beta
    //   Calculate delta and mu and update sigma
    
}

// This will be called in parallel, so cannot modify any shared values in localParams
void updateRegionParams(int region, const gslVector &localParamValues, const std::vector<updateableParam> &localParams) {
    
}
