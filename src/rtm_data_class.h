#ifndef HEADER_dataclass_
#define HEADER_dataclass_

#include <vector>
#include <gsl/gsl_vector.h>
#include "RTM_Header.h"
#include "distributions.h"

using namespace std;
using std::string;

// MODEL FOR LIKELIHOOD BOUNDS
struct likelihood_bounds{
  int lower;
  int upper;
};

// DEFINE A CLASS FOR HOLDING DATA STRUCTURES

class rtmData
{
 private:
  likelihood_bounds lbounds;
  data_type likelihood_type;
  gsl_matrix* counts;
  gsl_matrix* denoms;
  gsl_vector* data_population;
  gsl_vector* popn_weights;
  vector<int> col_grouping;
 public:
  rtmData();
  rtmData(const likelihood_bounds&, const data_type&);
  ~rtmData();
  data_type get_likelihood_type();
  void alloc(int&, int&);
  void setDim(const bool&, const string&, const string&, const int&);
  void read(const string&, const string&);
  void data_population_sizes(const gsl_vector*);
  void normalise(const gsl_vector*);
  double lfx(gsl_matrix*, const gsl_matrix*);
  int getDim1();
  int getDim2();
  const gsl_vector* access_weights() const;
  const vector<int> access_groups() const;
};

#endif
