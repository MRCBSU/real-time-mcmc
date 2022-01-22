#ifndef HEADER_modelstate_
#define HEADER_modelstate_

#include <cstddef>
#include <vector>
#include "RTM_Header.h"

using namespace std;

// DEFINE A CLASS FOR MONITORING THE END STATE OF THE EPIDEMIC
// (FOR OUTPUT TO SUBSEQUENT SMC ANALYSES)

class model_state
{
 private:
  vector<double> S;
  vector<double> E_1;
  vector<double> E_2;
  vector<double> I_1;
  vector<double> I_2;
  vector<double> R_pos;
  vector<double> R_neg;
  vector<double> W;
  vector<double> p_lambda;
 public:
  model_state();
  model_state(const unsigned int);
  model_state(const model_state&);
  model_state& operator=(const model_state&);
  void fill(const gsl_vector*, const gsl_vector*, const gsl_vector*, const gsl_vector*, const gsl_vector*,
	    const gsl_vector*, const gsl_vector*, const gsl_vector*, const gsl_vector*);
  bool write(ofstream&);
};

#endif
