// LOADING IN USER-DEFINED HEADER FILE, CONTAINING ALL STANDARD HEADER FILES
#include "RTM_StructDefs.h"
#include "gsl_mat_ext.h"

using namespace std;
using std::string;

void mixing_model_alloc(struct mixing_model &l_MIXMOD,
			const int nbreakpoints,
			const int num_scalants,
			const int num_strata)
{

  l_MIXMOD.num_breakpoints = nbreakpoints;
  l_MIXMOD.breakpoints = (nbreakpoints > 0) ? gsl_vector_int_alloc(nbreakpoints) : 0;
  l_MIXMOD.MIXMAT = new gsl_matrix* [nbreakpoints + 1];
  l_MIXMOD.MIXMAT_scaled = new gsl_matrix* [nbreakpoints + 1];
  l_MIXMOD.MIXMAT_param = new gsl_matrix_int* [nbreakpoints + 1];
  l_MIXMOD.scalants = gsl_vector_alloc(num_scalants);
  l_MIXMOD.eval_dominant = gsl_vector_alloc(nbreakpoints + 1);
  l_MIXMOD.evector_MIXMAT_normalised = new gsl_vector* [nbreakpoints + 1];
  for(int inti = 0; inti < nbreakpoints + 1; inti++){
    l_MIXMOD.MIXMAT[inti] = gsl_matrix_alloc(num_strata, num_strata);
    l_MIXMOD.MIXMAT_scaled[inti] = gsl_matrix_alloc(num_strata, num_strata);
    l_MIXMOD.MIXMAT_param[inti] = gsl_matrix_int_alloc(num_strata, num_strata);
    l_MIXMOD.evector_MIXMAT_normalised[inti] = gsl_vector_alloc(num_strata);
  }

}

void mixing_model_free(struct mixing_model &l_MIXMOD)
{

  for(int inti = 0; inti < l_MIXMOD.num_breakpoints + 1; inti++){
    gsl_matrix_free(l_MIXMOD.MIXMAT[inti]);
    gsl_matrix_free(l_MIXMOD.MIXMAT_scaled[inti]);
    gsl_matrix_int_free(l_MIXMOD.MIXMAT_param[inti]);
    if(l_MIXMOD.evector_MIXMAT_normalised != NULL)
      gsl_vector_free(l_MIXMOD.evector_MIXMAT_normalised[inti]);
  }
  if(l_MIXMOD.evector_MIXMAT_normalised != NULL)
    delete [] l_MIXMOD.evector_MIXMAT_normalised;
  if(l_MIXMOD.eval_dominant != NULL)
    gsl_vector_free(l_MIXMOD.eval_dominant);
  if(l_MIXMOD.scalants != NULL)
    gsl_vector_free(l_MIXMOD.scalants);
  delete [] l_MIXMOD.MIXMAT;
  delete [] l_MIXMOD.MIXMAT_scaled;
  delete [] l_MIXMOD.MIXMAT_param;
  if(l_MIXMOD.num_breakpoints > 0) gsl_vector_int_free(l_MIXMOD.breakpoints);
  
}

// SET ONE MIXING_MODEL STRUCTURE TO CONTAIN (PRE-ALLOCATED) POINTERS TO VALUES
// IDENTICAL TO THOSE POINTED TO BY MEMBERS OF ANOTHER MIXING_MODEL STRUCTURE
void mixing_model_memcpy(mixing_model &mix_dest, const mixing_model& mix_src)
{
  // function presumes all memory has been prior allocated
  if(mix_src.num_breakpoints > 0)
    gsl_vector_int_memcpy(mix_dest.breakpoints, mix_src.breakpoints);
  if(mix_src.eval_dominant != 0)
    gsl_vector_memcpy(mix_dest.eval_dominant, mix_src.eval_dominant);
  if(mix_src.scalants != 0)
    gsl_vector_memcpy(mix_dest.scalants, mix_src.scalants);
  for(register int i = 0; i <= mix_src.num_breakpoints; i++){
    gsl_matrix_memcpy(mix_dest.MIXMAT[i], mix_src.MIXMAT[i]);
    gsl_matrix_memcpy(mix_dest.MIXMAT_scaled[i], mix_src.MIXMAT_scaled[i]);
    gsl_matrix_int_memcpy(mix_dest.MIXMAT_param[i], mix_src.MIXMAT_param[i]);
    if(mix_src.evector_MIXMAT_normalised != 0)
      gsl_vector_memcpy(mix_dest.evector_MIXMAT_normalised[i], mix_src.evector_MIXMAT_normalised[i]);
  }
}


int maximal_index(const mixing_model src_mix)
{
  int max_index = 0;
  for(int int_i = 0; int_i < src_mix.num_breakpoints; int_i++)
    max_index = FN_MAX(max_index, gsl_matrix_int_max(src_mix.MIXMAT_param[int_i]));

  return max_index;
}

// NEED FOR string_functions.cc??
// COUNTS INSTANCES OF STRING sep_string IN FILE filename
int count_instances(const char *filename, const char *sep_string)
{

  string input;
  int cat_appearances = 0;

  fn_load_file(&input, filename);

  for(int i = input.find("breakpoint", 0); i != string::npos; i = input.find("breakpoint", i))
    {
      cat_appearances++;
      i++;
    }

  return cat_appearances;
}

int max_mm_param_int(const char *filename, const int nbreakpoints)
{

  FILE* fl_filename;
  register int i, max_int_param = 0, int_dummy;
  gsl_matrix_int* tempmat = gsl_matrix_int_alloc(NUM_AGE_GROUPS, NUM_AGE_GROUPS);
  char str_dummy[20];

  fl_filename = fopen(filename, "r");

  for(i = 0; i < nbreakpoints; i++){

    gsl_matrix_int_fscanf(fl_filename, tempmat);

    if(i != nbreakpoints)
      fscanf(fl_filename, "%s %d\n", str_dummy, &int_dummy);

    max_int_param = GSL_MAX(max_int_param, gsl_matrix_int_max(tempmat));

  }

  gsl_matrix_int_free(tempmat);
  fclose(fl_filename);

  return max_int_param;
}


void input_mixing_matrix_model(mixing_model &l_MIXMAT, const char* infile, const char* in_param_file)
{

  register FILE* age_mixing_matrix;
  register FILE* age_mixing_matrix_param;
  register char temp_string[25];
  register int inti, temp_int;
  int nbreakpoints;
  int num_scalants;
  double evalue; // dummy variable for storing unrequired dominant eigenvalues

  nbreakpoints = count_instances(infile, "breakpoint");

  num_scalants = max_mm_param_int(in_param_file, nbreakpoints) + 1;

  mixing_model_alloc(l_MIXMAT, nbreakpoints, num_scalants);

  age_mixing_matrix = fopen(infile, "r");
  age_mixing_matrix_param = fopen(in_param_file, "r");

  // READ IN CONTACT MATRICES

  for(int i = 0; i <= l_MIXMAT.num_breakpoints; i++){

    gsl_matrix_fscanf(age_mixing_matrix, l_MIXMAT.MIXMAT[i]);
    gsl_matrix_int_fscanf(age_mixing_matrix_param, l_MIXMAT.MIXMAT_param[i]);
    if(i != l_MIXMAT.num_breakpoints){
      fscanf(age_mixing_matrix, "%s", temp_string);
      fscanf(age_mixing_matrix, "%d", gsl_vector_int_ptr(l_MIXMAT.breakpoints, i));
      fscanf(age_mixing_matrix_param, "%s %d", temp_string, &temp_int);
    }

    // CALCULATE DOMINANT EIGENVALUE AND CALCULATE AND STORE DOMINANT EIGENVECTOR
    gsl_compute_max_evalue_evector2(gsl_vector_ptr(l_MIXMAT.eval_dominant, i),
				   l_MIXMAT.evector_MIXMAT_normalised[i],
				    l_MIXMAT.MIXMAT[i]);
  }

  fclose(age_mixing_matrix);
  fclose(age_mixing_matrix_param);

}

void mixing_matrix_parameterise(mixing_model& l_MIXMOD)
{
  int inti, intj, intk;
  gsl_matrix *matrix_i;

  for(inti = 0; inti <= l_MIXMOD.num_breakpoints; inti++){
    matrix_i = l_MIXMOD.MIXMAT_scaled[inti];
    for(intj = 0; intj < matrix_i->size1; intj++)
      for(intk = 0; intk < matrix_i->size2; intk++)
	gsl_matrix_set(matrix_i, intj, intk,
		       gsl_matrix_get(l_MIXMOD.MIXMAT[inti], intj, intk) * gsl_vector_get(l_MIXMOD.scalants, gsl_matrix_int_get(l_MIXMOD.MIXMAT_param[inti], intj, intk)));

  }

}
// MIXING_MATRIX_NORMALIZATION: SCALES THE TIME T MIXING MATRIX BY THE DOMINANT EIGENVALUE OF THE POPULATION SCALED MIXING MATRIX AT TIME 0
void mixing_matrix_scale_and_normalisation(
					   mixing_model& l_MIXMOD_ADJUSTED2, // FINAL OUTPUT
					   const gsl_vector* l_populationAG,
					   const bool evec_flag = true) // POPULATION PER AGE GROUP
{ 
 
  register gsl_matrix* MIXMAT_dummy = gsl_matrix_alloc(l_MIXMOD_ADJUSTED2.MIXMAT[0]->size1, l_MIXMOD_ADJUSTED2.MIXMAT[0]->size2);

  // for each mixing matrix used
  for(register int i = 0; i <= l_MIXMOD_ADJUSTED2.num_breakpoints; i++){

    if(i == 0){ // TIME-0 NEXT GENERATION MATRIX
      gsl_matrix_memcpy(MIXMAT_dummy, l_MIXMOD_ADJUSTED2.MIXMAT_scaled[i]);

      for (register int b = 0; b < l_MIXMOD_ADJUSTED2.MIXMAT[i]->size2; ++b)
	{
	  gsl_vector_view MIXMAT_dummy_col = gsl_matrix_column(MIXMAT_dummy, b);
	  gsl_vector_mul(&MIXMAT_dummy_col.vector,
			 l_populationAG);
	}

      // now get dominant eigenvalue and eigenvector of these matrices
      gsl_compute_max_evalue_evector2(gsl_vector_ptr(l_MIXMOD_ADJUSTED2.eval_dominant, i),
				      l_MIXMOD_ADJUSTED2.evector_MIXMAT_normalised[i],
				      MIXMAT_dummy,
				      evec_flag);
    }

    // final matrices are to be scaled by the dominant eigenvalue
    gsl_matrix_scale(l_MIXMOD_ADJUSTED2.MIXMAT_scaled[i], 1.0 / gsl_vector_get(l_MIXMOD_ADJUSTED2.eval_dominant, 0));

  }

  gsl_matrix_free(MIXMAT_dummy);

} 

// SHORT FUNCTION TO PICK OUT THE REQUIRED MIXING MATRIX
int mix_timecat(const int time_steps, const gsl_vector_int* breakpoints, const double time_steps_per_day){

  int out_cat = 0;

  for(; out_cat < breakpoints->size;){

    if(time_steps_per_day * ((double) gsl_vector_int_get(breakpoints, out_cat)) <= (double) time_steps){
      out_cat++;}
    else {
      break;}

  }

  return out_cat;
}
