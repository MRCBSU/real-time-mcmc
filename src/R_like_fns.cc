#include "R_like_fns.h"
#include "gsl_mat_ext.h"

using std::numeric_limits;

void R_gl_fac_sublevel(gsl_vector_int *factor_sequence, const int start_level,
                       const int start_sublevel, const int num_sublevels,
                       const int sequence_length) {

    int counter1 = start_level, counter2 = start_sublevel;
    register int i;

    for (i = 0; i < sequence_length; i++) {

        gsl_vector_int_set(factor_sequence, i, counter1);

        if (++counter2 > num_sublevels) {
            ++counter1;
            counter2 = 1;
        }
    }
}

double R_univariate_prior_log_density_nonnorm(const double &x,
                                              const distribution_type &dist,
                                              const gsl_vector *params) {

    switch (dist) {
    case cCONSTANT:
        return 0;
    case cGAMMA: // FIRST PARAMETER IS THE SHAPE PARAMETER, SECOND PARAMETER IS
                 // THE RATE PARAMETER
        return ((gsl_vector_get(params, 0) - 1) * gsl_sf_log(x)) -
               (gsl_vector_get(params, 1) * x);
    case cBETA: // FIRST PARAMETER IS THE FIRST SHAPE PARAMETER 'a',
        // SECOND PARAMETER IS THE SECOND SHAPE PARAMETER, 'b'
        return (((gsl_vector_get(params, 0) - 1) * gsl_sf_log(x)) +
                ((gsl_vector_get(params, 1) - 1) * gsl_sf_log(1 - x)));
    case cNORMAL:     // FIRST PARAMETER IS THE MEAN, SECOND PARAMETER IS THE
                      // STANDARD DEVIATION
    case cHALFNORMAL: // AS THE NORMAL DISTRIBUTION, BUT CONSTRAINED TO BE
                      // GREATER THAN ZERO
        return gsl_pow_2((x - gsl_vector_get(params, 0)) /
                         gsl_vector_get(params, 1)) /
               (-2.0);
    case cMVNORMAL: // SHOULD BE EVALUATED BY A DIFFERENT FUNCTION
        MULTIVARIATE_DISTRIBUTION;
    case cUNIFORM: // FIRST PARAMETER IS THE LOWER BOUND, THE SECOND PARAMETER
                   // THE UPPER BOUND
        return (x < gsl_vector_get(params, 0) || x > gsl_vector_get(params, 1))
                   ? GSL_NEGINF
                   : 0;
    default:
        DEFAULT_NO_DISTRIBUTION;
    }
}

double R_univariate_prior_log_density(const double &x,
                                      const distribution_type &dist,
                                      const gsl_vector *params) {

    switch (dist) {
    case cCONSTANT:
        return 0;
    case cGAMMA:
        return gsl_sf_log(gsl_ran_gamma_pdf(x, gsl_vector_get(params, 0),
                                            1 / gsl_vector_get(params, 1)));
    case cBETA:
        return gsl_sf_log(gsl_ran_beta_pdf(x, gsl_vector_get(params, 0),
                                           gsl_vector_get(params, 1)));
    case cNORMAL:
        return gsl_sf_log(gsl_ran_gaussian_pdf(x - gsl_vector_get(params, 0),
                                               gsl_vector_get(params, 1)));
    case cHALFNORMAL:
        return gsl_sf_log(2 *
                          gsl_ran_gaussian_pdf(x - gsl_vector_get(params, 0),
                                               gsl_vector_get(params, 1)));
    case cMVNORMAL: // SHOULD BE EVALUATED BY A DIFFERENT FUNCTION
        MULTIVARIATE_DISTRIBUTION;
    case cUNIFORM: // FIRST PARAMETER IS THE LOWER BOUND, THE SECOND PARAMETER
                   // THE UPPER BOUND
        return gsl_sf_log(gsl_ran_flat_pdf(x, gsl_vector_get(params, 0),
                                           gsl_vector_get(params, 1)));
    default:
        DEFAULT_NO_DISTRIBUTION;
    }
}

double R_univariate_prior_log_density_ratio(const double &x1, const double &x2,
                                            const distribution_type &dist,
                                            const gsl_vector *params) {
    bool flag1, flag2;
    // CALCULATES THE RATIO OF DENSITIES PI(X1) / PI(X2)
    switch (dist) {
    case cCONSTANT:
        return 0;
    case cGAMMA: // FIRST PARAMETER IS THE SHAPE PARAMETER, SECOND PARAMETER IS
                 // THE RATE PARAMETER
        return ((gsl_vector_get(params, 0) - 1) * gsl_sf_log(x1 / x2)) -
               (gsl_vector_get(params, 1) * (x1 - x2));
    case cBETA: // FIRST PARAMETER IS THE FIRST SHAPE PARAMETER 'a',
        // SECOND PARAMETER IS THE SECOND SHAPE PARAMETER, 'b'
        return (((gsl_vector_get(params, 0) - 1) * gsl_sf_log(x1 / x2)) +
                ((gsl_vector_get(params, 1) - 1) *
                 gsl_sf_log((1 - x1) / (1 - x2))));
    case cNORMAL:     // FIRST PARAMETER IS THE MEAN, SECOND PARAMETER IS THE
                      // STANDARD DEVIATION
    case cHALFNORMAL: // SEE ABOVE
        return (-1 / (2 * gsl_pow_2(gsl_vector_get(params, 1)))) *
               (gsl_pow_2(x1 - gsl_vector_get(params, 0)) -
                gsl_pow_2(x2 - gsl_vector_get(params, 0)));
    case cMVNORMAL: // SHOULD BE EVALUATED BY A DIFFERENT FUNCTION
        MULTIVARIATE_DISTRIBUTION;
    case cUNIFORM:
        flag1 = (x1 >= gsl_vector_get(params, 0)) &&
                (x1 <= gsl_vector_get(params, 1));
        flag2 = (x2 >= gsl_vector_get(params, 0)) &&
                (x1 <= gsl_vector_get(params, 1));
        return (flag1) ? ((flag2) ? 0 : GSL_POSINF)
                       : ((flag2) ? GSL_NEGINF : 0); // CHECK GSL_INF....
    default:
        DEFAULT_NO_DISTRIBUTION;
    }
}

int num_parameters_by_distribution(const distribution_type &dist,
                                   const int &dimension) {
    switch (dist) {
    case cCONSTANT:
        return 0;
    case cGAMMA:
    case cBETA:
    case cNORMAL:
    case cHALFNORMAL:
    case cUNIFORM:
        return 2;
    case cMVNORMAL:
        return (dimension + 1);
    default:
        DEFAULT_NO_DISTRIBUTION;
    }
}

double inverse_rw_transformation(const double &in_y, const char &trunc_flag,
                                 const double &a, const double &b) {
    switch (trunc_flag) {
    case cNO_TRUNC:
        return in_y;
    case cTRUNC_LO:
        return a + gsl_sf_exp(in_y);
    case cTRUNC_UP:
        return b - gsl_sf_exp(-in_y);
    case cTRUNC:
        return ((b * gsl_sf_exp(in_y)) + a) / (1 + gsl_sf_exp(in_y));
    }
}
// IMPORTANT FUNCTION FOR GENERATING CANDIDATE PROPOSALS FROM A POSSIBLY
// TRUNCATED RANDOM WALK
double R_proposenew2(double &newparam, const double &oldparam,
                     const double &propvar, gsl_rng *r, const char &trunc_flag,
                     const double &trunc_lo = 0.0,
                     const double &trunc_up = 1.0) {
    double param_trans = 0, ratio = 1;
    double proposenew_sd = sqrt(propvar);

    // Forward transformation
    if (trunc_flag) {
        if ((trunc_flag == cTRUNC) || (trunc_flag == cTRUNC_LO))
            param_trans = gsl_sf_log(oldparam - trunc_lo);
        if ((trunc_flag == cTRUNC) || (trunc_flag == cTRUNC_UP))
            param_trans -= gsl_sf_log(trunc_up - oldparam);
    } else
        param_trans = oldparam;
    // Generate a new value in the transformed space
    double newparam_trans =
        param_trans + gsl_ran_gaussian_ziggurat(r, proposenew_sd);
    // Inverse transformation
    newparam = inverse_rw_transformation(newparam_trans, trunc_flag, trunc_lo,
                                         trunc_up);
    // Get acceptance contribution
    if (trunc_flag) {
        if ((trunc_flag == cTRUNC) || (trunc_flag == cTRUNC_LO))
            ratio = (newparam - trunc_lo) / (oldparam - trunc_lo);
        if ((trunc_flag == cTRUNC) || (trunc_flag == cTRUNC_UP))
            ratio *= (trunc_up - newparam) / (trunc_up - oldparam);
    }
    return ratio;
}

double R_proposenew(double &newparam, const double &oldparam,
                    const double &propvar, gsl_rng *r, const char &trunc_flag,
                    const double &trunc_lo = 0.0,
                    const double &trunc_up = 1.0) {
    // proposes a new parameter value newparam, from a random walk proposal with
    // mean oldparam and variance propvar. The function returns the proposal
    // acceptance ratio, equal to the ratio of the two normalisations where
    // these are applicable (in the case of truncation), else 1. Optional
    // arguments should indicate the lower and upper limits of the truncation
    // respectively
    double U, ratio, proposenew_sd;
    double proposenew_upper, proposenew_lower, deva, devb;

    proposenew_sd = sqrt(propvar);

    if (trunc_flag) {
        if ((trunc_flag == cTRUNC) || (trunc_flag == cTRUNC_LO)) {
            // REPARAMETRIZATION FOR UNIFORM GAUSSIAN
            proposenew_lower = (trunc_lo - oldparam) / proposenew_sd;
        }
        if ((trunc_flag == cTRUNC) || (trunc_flag == cTRUNC_UP)) {
            // REPARAMETRIZATION FOR UNIFORM GAUSSIAN
            proposenew_upper = (trunc_up - oldparam) / proposenew_sd;
        }
        // COMPUTATION OF CDF AT THE RIGHT AND LEFT BOUNDARIES OF THE TRUNCATED
        // DISTRIBUTION
        devb = (trunc_flag == cTRUNC_LO)
                   ? 1.0
                   : gsl_cdf_ugaussian_P(proposenew_upper);
        deva = (trunc_flag == cTRUNC_UP)
                   ? 0.0
                   : gsl_cdf_ugaussian_P(proposenew_lower);
    } else {
        devb = 1.0;
        deva = 0.0;
    }
    // propose the new value
    ratio = devb - deva;

    U = gsl_rng_uniform_pos(r);
    U = deva + (U * ratio);

    newparam = oldparam + (proposenew_sd * gsl_cdf_ugaussian_Pinv(U));
    // COMPUTATION OF THE RATIO OF THE PROPOSED GIVEN THE CURRENT OVER THE
    // CURRENT GIVEN THE PROPOSED PROPOSAL DISTRIBUTIONS
    if (trunc_flag) {
        if ((trunc_flag == cTRUNC) || (trunc_flag == cTRUNC_UP)) {
            if (newparam > trunc_up)
                newparam = trunc_up - DBL_EPSILON;
            proposenew_upper = (trunc_up - newparam) / proposenew_sd;
        }
        if ((trunc_flag == cTRUNC) || (trunc_flag == cTRUNC_LO)) {
            if (newparam < trunc_lo)
                newparam = trunc_lo + DBL_EPSILON;
            proposenew_lower = (trunc_lo - newparam) / proposenew_sd;
        }
        devb = (trunc_flag == cTRUNC_LO)
                   ? 1.0
                   : gsl_cdf_ugaussian_P(proposenew_upper);
        deva = (trunc_flag == cTRUNC_UP)
                   ? 0.0
                   : gsl_cdf_ugaussian_P(proposenew_lower);
    }
    ratio /= (devb - deva);
    return ratio;
}
//

// WRAPPER FUNCTION FOR R_proposenew, RETURNS LOG-ACCEPTANCE RATIO CONTRIBUTION
double random_walk_proposal(double &proposal, const double &currentval,
                            const distribution_type &dist,
                            const double &proposal_var, gsl_rng *r,
                            const double &trunclo, const double &trunchi) {
    int lower_bound, upper_bound;
    double out;
    switch (dist) {
    case cMVNORMAL:
    case cNORMAL:
        return gsl_sf_log(
            R_proposenew2(proposal, currentval, proposal_var, r, cNO_TRUNC));
    case cGAMMA:
    case cHALFNORMAL:
        out = R_proposenew2(proposal, currentval, proposal_var, r, cTRUNC_LO,
                            std::numeric_limits<double>::min());
        return ((out == 0) ? GSL_NEGINF : gsl_sf_log(out));
    case cBETA:
        out = R_proposenew2(proposal, currentval, proposal_var, r, cTRUNC);
        return ((out == 0) ? GSL_NEGINF : gsl_sf_log(out));
    case cUNIFORM:
        out = R_proposenew2(proposal, currentval, proposal_var, r, cTRUNC,
                            trunclo, trunchi);
        return ((out == 0) ? GSL_NEGINF : gsl_sf_log(out));
    case cCONSTANT:
        proposal = currentval;
        return 0;
    default:
        DEFAULT_NO_DISTRIBUTION;
    }
}

double variance_inflation_factor(const double &proposal_variance,
                                 const distribution_type &dist,
                                 const gsl_vector *prior_params) {

    switch (dist) {
    case cBETA:
        return FN_MIN(proposal_variance, MAX_ZERO_ONE_PROPOSAL_VARIANCE);
    case cUNIFORM:
        return FN_MIN(proposal_variance, MAX_ZERO_ONE_PROPOSAL_VARIANCE *
                                             (gsl_vector_get(prior_params, 1) -
                                              gsl_vector_get(prior_params, 0)));
    case cGAMMA:
    case cNORMAL:
    case cMVNORMAL:
    case cHALFNORMAL:
        return FN_MIN(proposal_variance, MAX_UNBOUNDED_PROPOSAL_VARIANCE);
    default:
        return proposal_variance;
    }
}

double R_inverse_link_function(double x, const link_function g) {

    switch (g) {
    case cIDENTITY:
        return x;
    case cLOG:
        return gsl_sf_exp(x);
    case cLOGIT:
        return gsl_sf_exp(x) / (1 + gsl_sf_exp(x));
    default:
        DEFAULT_NO_LINK;
    }
}

void R_generalised_linear_regression(gsl_vector *out_y,
                                     const gsl_matrix *design_X,
                                     const gsl_vector *param_beta,
                                     link_function g) {

    // multiply the parameter vector by the design matrix
    gsl_matrix_multiplication_vector_dbl(out_y, design_X, param_beta);

    // apply the inverse of the selected link function
    for (int int_i = 0; int_i < out_y->size; int_i++)
        gsl_vector_set(
            out_y, int_i,
            R_inverse_link_function(gsl_vector_get(out_y, int_i), g));
}

// use a vector of breakpoints to aggregate a gsl_vector into sums by factor
void R_by_sum_gsl_vector_mono_idx(
    gsl_vector *oVec, gsl_vector *oVec_wts, const gsl_vector *iVec,
    const vector<int>
        &idx_breaks) { // Aggregates the population size over observed strata
                       // (oVec) and outputs if not NULL the weights of each
                       // input strata within the
    // observed strata.
    vector<int>::const_iterator itIdx = idx_breaks.begin();
    int int_width = *itIdx, offset = 0;
    gsl_vector_view iVec_sub;
    gsl_vector *iVec_cpy;
    if (oVec_wts != NULL)
        if (oVec_wts->size != iVec->size) {
            gsl_vector_free(oVec_wts);
            oVec_wts = gsl_vector_calloc(iVec->size);
        }
    if (iVec->size < idx_breaks.size()) {
        cout
            << "R_by_sum_stl_vector_mono_idx: Indices exceed maximum length of "
               "input vector."
            << endl;
        cout << "Input vector padded with zeros." << endl;
        iVec_cpy = gsl_vector_alloc(idx_breaks.size());
        gsl_vector_set_zero(iVec_cpy);
        iVec_sub = gsl_vector_subvector(iVec_cpy, 0, iVec->size);
        gsl_vector_memcpy(&iVec_sub.vector, iVec);
    } else {
        iVec_cpy = gsl_vector_alloc(iVec->size);
        gsl_vector_memcpy(iVec_cpy, iVec);
    }
    for (int int_i = 0; int_i < oVec->size; itIdx++) {
        iVec_sub = gsl_vector_subvector(iVec_cpy, offset, int_width);
        gsl_vector_set(oVec, int_i++,
                       gsl_vector_sum_elements(&iVec_sub.vector));
        if (oVec_wts != NULL) {
            gsl_vector_view oVec_wts_sub =
                gsl_vector_subvector(oVec_wts, offset, int_width);
            gsl_vector_memcpy(&oVec_wts_sub.vector, &iVec_sub.vector);
            gsl_vector_scale(&oVec_wts_sub.vector,
                             1 / gsl_vector_get(oVec, int_i - 1));
        }
        offset += int_width;
        int_width = (((itIdx + 1) == idx_breaks.end()) ? iVec_cpy->size
                                                       : *(itIdx + 1)) -
                    *itIdx;
    }
    gsl_vector_free(iVec_cpy);
}
