#include "rtm_data_class.h"
#include "R_like_fns.h"
#include "RTM_FunctDefs.h"
#include "stl_vec_ext.h"

using namespace std;
using std::string;

#define ERROR_TEXT(str_err_msg, str_err_file)    printf(str_err_msg, str_err_file); \
  exit(2);

#define NULL_ALLOC      denoms = NULL;		\
  counts = NULL;				\
  data_population = NULL;			\
  popn_weights = NULL;

/////// ///////
rtmData::rtmData()
{
  lbounds.lower = 0;
  lbounds.upper = 100;
  likelihood_type = cPOISSON;
  NULL_ALLOC;
}
rtmData::rtmData(const likelihood_bounds& lfx_range, const data_type& lfx)
{
  lbounds = lfx_range;
  likelihood_type = lfx;
  NULL_ALLOC;
}
rtmData::~rtmData()
{
  if(counts != NULL)
    gsl_matrix_free(counts);
  if(denoms != NULL)
    gsl_matrix_free(denoms);
  if(popn_weights != NULL)
    gsl_vector_free(popn_weights);
}
data_type rtmData::get_likelihood_type()
{
  return likelihood_type;
}
void rtmData::alloc(int& t, int& a)
{
  counts = gsl_matrix_alloc(t, a);
  denoms = gsl_matrix_alloc(t, a);
  data_population = gsl_vector_alloc(a);
  if(a < NUM_AGE_GROUPS){
    if(likelihood_type == cBINOMIAL)
      popn_weights = gsl_vector_alloc(NUM_AGE_GROUPS);
  }
}
void rtmData::setDim(const bool& missing_data, const string& str_source, const string& str_property, const int& smax)
{
  int dim_s;
  col_grouping.clear();

  if(!missing_data){
    for(int int_i = 0; int_i < smax; col_grouping.push_back(int_i++)){};
    dim_s = smax;
  } else {
    string temp_string;
    if (!read_string_from_instruct(temp_string, str_property, str_source)) {
      std::cerr << "Could not read " << str_property << " from " << str_source << std::endl;
      exit(2);
    }
    stl_vector_sscanf<int>(temp_string, col_grouping);
    dim_s = col_grouping.size() + 1;
  }
  this->alloc(lbounds.upper, dim_s);
}
void rtmData::read(const string& denomfile, const string& countfile)
{
  // Are denominators provided?
  if(denomfile.compare("") != 0)
    // YES: read denominator input file.
    data_matrices_fscanf(denomfile.c_str(), denoms, 1);
  else { // NO: are denominators required
    if(likelihood_type == cBINOMIAL){
      // YES: error.
      ERROR_TEXT("Required denominator file %s not found\n", denomfile.c_str());
    } else gsl_matrix_set_all(denoms, 1.0); // Default value of 1.. everyone is observable.
  }
  
  // Are counts provided?
  if(countfile.compare("") != 0)
    //yes: read count input file.
    data_matrices_fscanf(countfile.c_str(), counts, 1);
  else { // ERROR!
    ERROR_TEXT("Required count file %s not found\n", countfile.c_str());
  }
}
void rtmData::data_population_sizes(const gsl_vector* ptr_vec_popn)
{
  // data_population_sizes *should* have length denoms->size2
  // ptr_vec_popn doesn't necessarily have the requisite length, might need to aggregate over col_grouping
  if(ptr_vec_popn->size != denoms->size2){
    if(ptr_vec_popn->size > denoms->size2){
      R_by_sum_gsl_vector_mono_idx(data_population, popn_weights, ptr_vec_popn, col_grouping);
    }
    else {
      cout << "Input vector shorter than output vector, cannot contract over indices." << endl;
      exit(2);
    }
  } else gsl_vector_memcpy(data_population, ptr_vec_popn);
}
void rtmData::normalise(const gsl_vector* ptr_vec_popn)
{
  // Don't do this if working with sample data
  if(likelihood_type != cBINOMIAL){
    data_population_sizes(ptr_vec_popn);
    // Loop through the rows of the denominator matrix, normalising each.
    for(int int_t = 0; int_t < denoms->size1; int_t++){
      gsl_vector_view denom_row = gsl_matrix_row(denoms, int_t);
      gsl_vector_div(&denom_row.vector, data_population);
    }
  }
}
double rtmData::lfx(gsl_matrix* mat_expected, const gsl_matrix* mat_dispersion)
{ // Note this function will alter the contents of mat_expected after scaling by denoms.
  gsl_matrix_view N = gsl_matrix_submatrix(denoms, lbounds.lower - 1, 0, lbounds.upper - lbounds.lower + 1, denoms->size2);
  gsl_matrix_view mu = gsl_matrix_submatrix(mat_expected, lbounds.lower - 1, 0, lbounds.upper - lbounds.lower + 1, mat_expected->size2);
  if(likelihood_type != cBINOMIAL)
    gsl_matrix_mul_elements(&mu.matrix, &N.matrix);
  gsl_matrix_view X = gsl_matrix_submatrix(counts, lbounds.lower - 1, 0, lbounds.upper - lbounds.lower + 1, counts->size2);
  switch(likelihood_type){
  case cPOISSON :
    return fn_log_lik_countdata(&X.matrix, &mu.matrix);
  case cNEGBIN :
    {
      gsl_matrix_const_view eta = gsl_matrix_const_submatrix(mat_dispersion, lbounds.lower - 1, 0, lbounds.upper - lbounds.lower + 1, mat_dispersion->size2);
      return fn_log_lik_negbindata(&X.matrix, &mu.matrix, &eta.matrix);
    }
  case cBINOMIAL :
    return fn_log_lik_positivity(&N.matrix, &X.matrix, &mu.matrix);
  default :
    DEFAULT_NO_LIKELIHOOD;
  }
}
int rtmData::getDim1(){
  return counts->size1;
}
int rtmData::getDim2(){
  return counts->size2;
}
const gsl_vector* rtmData::access_weights() const {
  return popn_weights;
}
const vector<int> rtmData::access_groups() const {
  return col_grouping;
}
