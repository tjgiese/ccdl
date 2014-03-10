/*-------------------------------------------------------|
|  NIST SPARSE BLAS v. 0.9 (Sat Jul 6 14:27:21 EDT 1996) |
|                                                        |
|  Authors:                                              |
|     Karin A. Remington and Roldan Pozo                 |
|     National Institute of Standards and Technology     |
|                                                        |
|  Based on the interface standard proposed in:          | 
|   "A Revised Proposal for a Sparse BLAS Toolkit" by    |
|    S. Carney and K. Wu -- University of Minnesota      |
|    M. Heroux and G. Li -- Cray Research                |  
|    R. Pozo and K.A. Remington -- NIST                  |
|                                                        |
|  Contact:                                              |
|     Karin A. Remington, email: kremington@nist.gov     |
--------------------------------------------------------*/

#include "spblas.h"

void ScaleArray_double(int m, int n, double *c, int ldc, const double beta)
{
  /* C <- beta*C         */
  /* C is m x n          */
  /* beta is a scalar    */

  int i, j;

  if (beta == 1.0)
    return;

  if (beta == 0.0) {
    if (n == 1)
      for (j = 0; j < m; j++)
	c[j] = 0.0;
    else
      for (i = 0; i < n; i++){
	for (j = 0; j < m; j++)
	  c[j] = 0.0;
        c+=ldc;
      }
  } else {
    if (n == 1)
      for (j = 0; j < m; j++)
	c[j] *= beta;
    else
      for (i = 0; i < n; i++){
	for (j = 0; j < m; j++)
	  c[j] *= beta;
        c+=ldc;
      }
  }
}

/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/

void ScaleArray_float(int m, int n, float *c, int ldc, const float beta)
{
  /* C <- beta*C          */
  /* C is m x n           */
  /* beta is a scalar     */

  int i, j;

  if (beta == 1.0)
    return;

  if (beta == 0.0) {
    if (n == 1)
      for (j = 0; j < m; j++)
	c[j] = 0.0;
    else
      for (i = 0; i < n; i++){
	for (j = 0; j < m; j++)
	  c[j] = 0.0;
        c+=ldc;
      }
  } else {
    if (n == 1)
      for (j = 0; j < m; j++)
	c[j] *= beta;
    else
      for (i = 0; i < n; i++){
	for (j = 0; j < m; j++)
	  c[j] *= beta;
        c+=ldc;
      }
  }
}

/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/

void ScaleArray_complex(int m, int n, float *c, int ldc, const float *beta)
{
  /* C <- beta*C          */
  /* C is m x n           */
  /* beta is a scalar     */

  float real, imag;
  int i, j, indexr, cldc = 2*ldc;

  if (beta[0] == 1.0 && beta[1] == 0.0 )
    return;

  if (beta[0] == 0.0 && beta[1] == 0.0 ) {
      for (j = 0; j < 2*m*n; j++)
	c[j] = 0.0;
  } else {
      for (i = 0; i < n; i++) {
	for (j = 0; j < m; j++) {
          indexr = j*2;
          real = c[indexr]*beta[0] - c[indexr+1]*beta[1];
          imag = c[indexr]*beta[1] + c[indexr+1]*beta[0];
	  c[indexr] = real;
	  c[indexr+1] = imag;
        }
        c+=cldc;
      }
  }
}

/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/

void ScaleArray_Z(int m, int n, double *c, int ldc, const double *beta)
{
  /* C <- beta*C          */
  /* C is m x n           */
  /* beta is a scalar     */

  double real, imag;
  int i, j, indexr, cldc = 2*ldc;

  if (beta[0] == 1.0 && beta[1] == 0.0 )
    return;

  if (beta[0] == 0.0 && beta[1] == 0.0 ) {
      for (j = 0; j < 2*m*n; j++)
	c[j] = 0.0;
  } else {
      for (i = 0; i < n; i++){
	for (j = 0; j < m; j++) {
          indexr = j*2;
          real = c[indexr]*beta[0] - c[indexr+1]*beta[1];
          imag = c[indexr]*beta[1] + c[indexr+1]*beta[0];
	  c[indexr] = real;
	  c[indexr+1] = imag;
        }
        c+=cldc;
      }
  }
}
