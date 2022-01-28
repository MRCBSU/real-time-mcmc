#include "RTM_modelstate.h"
#include "gsl_vec_ext.h"

using namespace std;

// public model_state functions
model_state::model_state() :
  S(STATE_DIM),
  E_1(STATE_DIM),
  E_2(STATE_DIM),
  I_1(STATE_DIM),
  I_2(STATE_DIM),
  R_pos(STATE_DIM),
  R_neg(STATE_DIM),
  W(STATE_DIM),
  WS(STATE_DIM),
  p_lambda(NUM_AGE_GROUPS)
{}
model_state::model_state(const unsigned int iSize) :
  S(iSize), E_1(iSize), E_2(iSize), I_1(iSize), I_2(iSize), R_pos(iSize), R_neg(iSize), W(iSize), WS(iSize), p_lambda(iSize)
{}
model_state::model_state(const model_state& oC) : 
  S(oC.S), E_1(oC.E_1), E_2(oC.E_2), I_1(oC.I_1), I_2(oC.I_2),
  R_pos(oC.R_pos), R_neg(oC.R_neg), W(oC.W), WS(oc.WS), p_lambda(oC.p_lambda)
{}
model_state& model_state::operator=(const model_state& oC)
{
  if (this != &oC) {
    S = oC.S;
    E_1 = oC.E_1;
    E_2 = oC.E_2;
    I_1 = oC.I_1;
    I_2 = oC.I_2;
    R_pos = oC.R_pos;
    R_neg = oC.R_neg;
    W = oC.W;
    WS = oC.WS;
    p_lambda = oC.p_lambda;
  };
  return *this;
}
void model_state::fill(const gsl_vector* inS, const gsl_vector* inE_1, const gsl_vector* inE_2, const gsl_vector* inI_1, const gsl_vector* inI_2, const gsl_vector* inR_pos, const gsl_vector* inR_neg, const gsl_vector* inW, const gsl_vector* inWS, const gsl_vector* inp_lambda)
{
  gsl_populate_vector(S, inS);
  gsl_populate_vector(E_1, inE_1);
  gsl_populate_vector(E_2, inE_2);
  gsl_populate_vector(I_1, inI_1);
  gsl_populate_vector(I_2, inI_2);
  gsl_populate_vector(R_pos, inR_pos);
  gsl_populate_vector(R_neg, inR_neg);
  gsl_populate_vector(W, inW);
  gsl_populate_vector(WS, inWS);
  gsl_populate_vector(p_lambda, inp_lambda);
}
bool model_state::write(ofstream& outFile)
{
  if(!outFile.is_open())
    {
      cout << "model_state::write unable to access open file outFile" << endl;
      return false;
    }
  else
    {
      outFile.write(reinterpret_cast<char const*>(&S[0]),
		    S.size() * sizeof(double));
      outFile.write(reinterpret_cast<char const*>(&E_1[0]),
		    E_1.size() * sizeof(double));
      outFile.write(reinterpret_cast<char const*>(&E_2[0]),
		    E_2.size() * sizeof(double));
      outFile.write(reinterpret_cast<char const*>(&I_1[0]),
		    I_1.size() * sizeof(double));
      outFile.write(reinterpret_cast<char const*>(&I_2[0]),
		    I_2.size() * sizeof(double));
      outFile.write(reinterpret_cast<char const*>(&R_pos[0]),
		    R_pos.size() * sizeof(double));
      outFile.write(reinterpret_cast<char const*>(&R_neg[0]),
		    R_neg.size() * sizeof(double));
      outFile.write(reinterpret_cast<char const*>(&W[0]),
		    W.size() * sizeof(double));
      outFile.write(reinterpret_cast<char const*>(&WS[0]),
		    WS.size() * sizeof(double));
      outFile.write(reinterpret_cast<char const*>(&p_lambda[0]),
		    p_lambda.size() * sizeof(double));
      return true;
    }
}
