// LOADING IN USER-DEFINED HEADER FILE, CONTAINING ALL STANDARD HEADER FILES REQUIRED
#include "string_fns.h"
#include "gsl_vec_ext.h"

void gsl_vector_int_set_1ton(gsl_vector_int* v, int start)
{
  for(int int_i = 0; int_i < v->size; int_i++)
    gsl_vector_int_set(v, int_i, start + int_i);
}

double gsl_vector_sum_elements(register const gsl_vector* input_vector) 
{
  register double sum_elements = 0;
  for (register int i=0; i<(input_vector->size); ++i)
    {
      sum_elements += gsl_vector_get(input_vector,i);
    }
  return sum_elements;
} 

int gsl_vector_int_sum_elements(register const gsl_vector_int* input_vector)
{
  register int sum_elements = 0;
  for (register int i=0; i<(input_vector->size); ++i)
    {
      sum_elements += gsl_vector_int_get(input_vector,i);
    }
  return sum_elements;

}

void gsl_vector_sscanf(const string in_string, gsl_vector* v, const char* chr_delim)
{

  int indx = 0;
  for(int i = 0; i < v->size; i++)
    gsl_vector_set(v, i, read_from_delim_string<double>(in_string, chr_delim, indx));

}

void gsl_vector_int_sscanf(const string in_string, gsl_vector_int* v, const char* chr_delim)
{

  int indx = 0;
  for(int i = 0; i < v->size; i++)
    gsl_vector_int_set(v, i, read_from_delim_string<int>(in_string, chr_delim, indx));

}

void gsl_vector_subvector_memcpy(gsl_vector* v_out, gsl_vector* v_in, const size_t offset)
{
  gsl_vector_view v_in_sub = gsl_vector_subvector(v_in, offset, v_out->size);
  gsl_vector_memcpy(v_out, &(v_in_sub.vector));
}

void gsl_vector_int_set_seq(gsl_vector_int* v, const int i_from, const int i_by)
{
  for(int i = 0; i < v->size; i++)
    gsl_vector_int_set(v, i, i_from + (i * i_by));
}

void gsl_vector_set_seq(gsl_vector* v, const double dbl_from, const double dbl_by)
{
  for(int i = 0; i < v->size; i++)
    gsl_vector_set(v, i, dbl_from + (i * dbl_by));
}

// GSL_DOUBLE_PRODUCT_OF_VECTOR_ELEMENTS RETURNS THE PRODUCT OF THE ELEMENTS OF A VECTOR
double gsl_double_product_of_vector_elements(const gsl_vector* vector) 
{
  double product = 1.0;
  for (register int i=0; i<(vector->size); ++i)
    {
      product *= gsl_vector_get(vector,i);
    }
  return product;
}

bool in_gsl_vector(gsl_vector* v, const double dbl_match)
{
  bool inflag = 0;
  for(int i = 0; (i < v->size) && !inflag; i++)
    {
      inflag = (gsl_vector_get(v, i) == dbl_match) ? 1 : 0;
    }
  return inflag;
}

bool in_gsl_vector_int(gsl_vector_int* v, const int int_match)
{
  bool inflag = 0;
  for(int i = 0; (i < v->size) && !inflag; i++)
    {
      inflag = (gsl_vector_int_get(v, i) == int_match) ? 1 : 0;
    }
  return inflag;
}
void gsl_populate_vector(vector<double>& v, const gsl_vector* gsl_v)
{
  if(v.size() == gsl_v->size)
    {
      for(int i = 0; i < gsl_v->size; i++)
	v[i] = gsl_vector_get(gsl_v, i);
    }
  else
    perror("Dimension mismatch in gsl_populate_vector");
}
