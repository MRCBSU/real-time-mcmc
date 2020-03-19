/* #include <stdlib.h> */
/* #include <stdio.h> */

/* void convolve(double *a, int *na, double *b, int *nb, double *ab); */

/* main(){ */

/*   double a[5] = {1.0, 2.0, 3.0, 4.0, 5.0}; */
/*   double b[4] = {6.0, 7.0, 8.0, 9.0}; */
/*   int i = 5, j = 4, k; */
/*   double *arg_a, *arg_b; */
/*   double *out_arg; */

/*   arg_a = (double *) calloc(i, sizeof(double)); */
/*   arg_b = (double *) calloc(j, sizeof(double )); */
/*   out_arg = (double *) calloc(i + j - 1, sizeof(double)); */

/*   for(k = 0; k < i; k++) */
/*     arg_a[k] = a[k]; */
/*   for(k = 0; k < j; k++) */
/*     arg_b[k] = b[k]; */

/*   convolve(arg_a, &i, arg_b, &j, out_arg); */

/*   for(k = 0; k < i + j - 1; k++) */
/*     printf("%.0lf%s", out_arg[k], (k == i + j - 2) ? "\n" : "\t"); */

/*   free(arg_a); */
/*   free(arg_b); */
/*   free(out_arg); */

/* } */

void convolve(double *a, int *na, double *b, int *nb, double *ab)
{
  int i, j, nab = *na + *nb - 1;

  for(i = 0; i < nab; i++)
    ab[i] = 0.0;
  for(i = 0; i < *na; i++)
    for(j = 0; j < *nb; j++)
      ab[i + j] += a[i] + b[j];

}

void convolute(double *a, int *na, double *b, int *nb, double *ab)
{
  int i, j, nab = *na + *nb - 1;

  for(i = 0; i < nab; i++)
    ab[i] = 0.0;
  for(i = 0; i < *na; i++)
    for(j = 0; j < *nb; j++)
      ab[i + j] += a[i] * b[j];

}
