// LOADING IN USER-DEFINED HEADER FILE, CONTAINING ALL STANDARD HEADER FILES REQUIRED
#include "string_fns.h"
#include "gsl_mat_ext.h"
#include "gsl_vec_ext.h"

// NEED TO FIND THE SUBROW AND SUBCOLUMN FUNCTIONS


// SET A ROW OF A MATRIX EQUAL TO A DOUBLE
void gsl_matrix_set_row_dbl(gsl_matrix* m, const size_t i_row, const double dbl_input)
{
  gsl_vector_view target_row = gsl_matrix_row(m, i_row);
  gsl_vector_set_all(&(target_row.vector), dbl_input);
}
// SET A COLUMN OF A MATRIX EQUAL TO A DOUBL
void gsl_matrix_set_column_dbl(gsl_matrix* m, const size_t i_col, const double dbl_input)
{
  gsl_vector_view target_col = gsl_matrix_column(m, i_col);
  gsl_vector_set_all(&(target_col.vector), dbl_input);
}


void gsl_matrix_set_subrow_dbl_memcpy(gsl_matrix* m, const size_t i_row, const size_t offset, const size_t n, const gsl_vector* dbl_input)
{
  gsl_vector_view target_subrow = gsl_matrix_subrow(m, i_row, offset, n);
  gsl_vector_memcpy(&(target_subrow.vector), dbl_input);
}

void gsl_matrix_set_subcolumn_dbl_memcpy(gsl_matrix* m, const size_t i_col, const size_t offset, const size_t n, const gsl_vector* dbl_input)
{
  gsl_vector_view target_subcolumn = gsl_matrix_subcolumn(m, i_col, offset, n);
  gsl_vector_memcpy(&(target_subcolumn.vector), dbl_input);
}

void gsl_matrix_set_subrow_dbl(gsl_matrix* m, const size_t i_row, const size_t offset, const size_t n, const double dbl_input)
{
  register int i;
  if(i_row < 0 || offset < 0 || i_row >= m->size1 || offset - 1 + n >= m->size2)
    perror("Error; trying to assign outside matrix limits in gsl_matrix_set_subrow_dbl\n");
  if(n < 0)
    perror("Error: trying to assign a negative number of entries in gsl_matrix_set_subrow_dbl\n");
  for(i = 0; i < n; i++)
    gsl_matrix_set(m, i_row, offset + i, dbl_input);
}

void gsl_matrix_set_subcolumn_dbl(gsl_matrix* m, const size_t j_col, const size_t offset, const size_t n, const double dbl_input)
{
  register int j;
  if(j_col < 0 || offset < 0 || j_col >= m->size2 || offset - 1 + n >= m->size1)
    perror("Error; trying to assign outside matrix limits in gsl_matrix_set_subcolumn_dbl\n");
  if(n < 0)
    perror("Error: trying to assign a negative number of entries in gsl_matrix_set_subcolumn_dbl\n");
  for(j = 0; j < n; j++)
    gsl_matrix_set(m, offset + j, j_col, dbl_input);
}


// GSL_MATRIX_MULTIPLICATION_VECTOR_DOUBLE; MULTIPLYING A MATRIX WITH A VECTOR IN DOUBLE FORMAT
void gsl_matrix_multiplication_vector_dbl(
					     gsl_vector* result, // PASSED BY REFERENCE AND CHANGES WITHIN THE FUNCTION
					     const gsl_matrix* A, 
					     const gsl_vector* b)
{
  // ERROR IDENTIFYING
  if(A->size2 != b->size) perror("Dimension mismatch in gsl_matrix_multiplication_vector_double");

  gsl_vector_set_zero(result); // INITIALIZATION

  for (register int i=0; i < A->size1; ++i)
    {
      for (register int j=0; j < A->size2; ++j)
	{
	  gsl_vector_set(result,i,gsl_vector_get(result,i) + (gsl_matrix_get(A,i,j)*gsl_vector_get(b,j)));
	}// FOR
    }// FOR
}

// MULTIPLIES THE ROWS OF MATRIX A BY VECTOR b
void gsl_matrix_multiplication_vector_product_by_row_dbl(
							 gsl_matrix* A,
							 const gsl_vector* b)
{
  if(A->size2 != b->size) perror("Dimension mismatch in gsl_matrix_multiplication_vector_product_by_row_dbl");

  gsl_vector_view A_row;
  for(int int_i = 0; int_i < A->size1; int_i++)
    {
      A_row = gsl_matrix_row(A, int_i);
      gsl_vector_mul(&A_row.vector, b);
    }
}


// COMPUTE_MAX_EVALUE_EVECTOR: COMPUTES THE MAXIMUM EIGENVALUE (evalue) AND CORRESPONDING EIGENVECTOR (evector) OF ANY GIVEN SYMMETRIC
// MATRIX M-THIS IS ONLY USED READING IN MIXING MATRIX
void gsl_compute_max_evalue_evector(
				    register double* evalue, // EIGENVALUE OUTPUT
				    gsl_vector* evector, // EIGENVECTOR OUTPUT IS NORMALIZED!
				    const gsl_matrix* M, // M IS ALWAYS A MATRIX OF DIMENSIONS POP_RANGES X POP_RANGES
				    const int  nu_iter, 
				    const double  tolerance)
{
  gsl_vector_set_all(evector, 1.0 / evector->size); // SET evector TO AN INITIAL VALUE SATISFYING SUM TO 1

  // SET THE OLD NORM TO A VERY LARGE POSITIVE NUMBER
  register double oldnorm = 0; 
  register gsl_vector* temp_prod = gsl_vector_alloc(evector->size);
 
  gsl_vector_set_zero(temp_prod);

  // MAIN COMPUTATIONS
  for (register int i=0; i < nu_iter; ++i)
    {

      gsl_matrix_multiplication_vector_dbl(temp_prod,M, evector); // MULTIPLY THE EIGENVECTOR WITH THE CONTACT MATRIX
      gsl_vector_memcpy(evector,temp_prod); // CREATE A COPY OF THE EIGENVECTOR
      *evalue = gsl_vector_sum_elements(evector); // CALCULATE THE EIGENVALUE
      gsl_vector_scale(evector, 1.0 / *evalue); // NORMALISE THE NEW EIGENVECTOR

      if ( (fabs( *evalue - oldnorm ) < tolerance )  && i > 0 )
	break;
      else oldnorm = *evalue;

    }//FOR

  // FINALLY YOU CAN NORMALIZE THE EIGENVECTOR THAT CORRESPONDS TO THE MAXIMUM EIGENVALUE:

  gsl_vector_free(temp_prod);
}


void gsl_matrix_sscanf(const string in_string, gsl_matrix* m, const char* chr_delim)
{

  int indx = 0;
  for(int i = 0; i < m->size1; i++)
    for(int j = 0; j < m->size2; j++)
      gsl_matrix_set(m, i, j, read_from_delim_string<double>(in_string, chr_delim, indx));

}

void gsl_matrix_int_sscanf(const string in_string, gsl_matrix_int* m, const char* chr_delim)
{

  int indx = 0;
  for(int i = 0; i < m->size1; i++)
    for(int j = 0; j < m->size2; j++)
      gsl_matrix_int_set(m, i, j, read_from_delim_string<int>(in_string, chr_delim, indx));

}
