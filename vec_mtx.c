#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>



/*--- Temporary matrix. ----*/
double   **T=NULL;

/*--------------------------------------------------------
 * Procedure to allocate an n X n matrix and return the 
 * pointer to the matrix.
 */
double **alloc_mtx(int n)
{
	double  **B;
	int     i;

  B = (double **) malloc(sizeof(double *)*n);
  for(i=0;i<n;i++) 
	B[i] = (double *) malloc(sizeof(double)*n);
  return (B);
}



/*----------------------------------------------------------------------
 * Procedure to make an identity matrix.
 *    A: the input matrix,
 *   n: dimension of the matrix.
 */
void make_identity_mtx(double **A, int n)
{
   int i, j;

   for(i=0;i<n;i++)
	   for(j=0;j<n;j++)
		   if(i!=j) A[i][j] = 0.0;
		   else A[i][j] = 1.0;
}

/*--------------------------------------------------------------------
 * Post-multiply a matrix by another matrix.
 *      A=  A*R,
 * Input:
 *      A: the destination matrix,
 *      R: the mtx to multiply with A
 *      n: dimension of the matrix.
 */
void post_mtx_mult(double **A, double **R, int n)
{
	int  i, j, k;

	if(T==NULL)  T=alloc_mtx(n);

	for(i=0;i<n;i++)
		for(j=0;j<n;j++){
			T[i][j] = 0.0;
			for(k=0;k<n;k++)
				T[i][j] += A[i][k]*R[k][j];
		}

	for(i=0;i<n;i++)
       for(j=0;j<n;j++) A[i][j] = T[i][j];
}

/*--------------------------------------------------------------------
 * Pre-multiply a matrix by another matrix.
 *      A=  R*A,
 * Input:
 *      A: the destination matrix,
 *      R: the mtx to multiply with A
 *      n: dimension of the matrix.
 */
void pre_mtx_mult(double **R, double **A, int n)
{
	int  i, j, k;

	if(T==NULL)  T=alloc_mtx(n);

	for(i=0;i<n;i++)
		for(j=0;j<n;j++){
			T[i][j] = 0.0;
			for(k=0;k<n;k++)
				T[i][j] += R[i][k]*A[k][j];
		}

	for(i=0;i<n;i++)
       for(j=0;j<n;j++) A[i][j] = T[i][j];
}


/*-------------------------------------------------------------------
 * Procedure to form a rotational mtx,
 *     R: the input matrix,
 *     p, q: indinces of the cos() and sin() values,
 *     c, s: cos and sine values,
 *     n: dimension of the matrix.
 */
void make_rotate_mtx(double **R,  int p, int q, double c, double s, int n)
{
    make_identity_mtx(R, n);
	R[p][p] = c;
	R[q][q] = c;
	R[p][q] = -s;
	R[q][p] = s;
}


/*-------------------------------------------------------------------
 * Procedure to transpose a matrix.
 *   R: the input matrix.
 *   n: dimension of the matrix.
 */
void transpose_mtx(double **R, int n)
{
	double  temp;
    int     i, j;

	for(i=0;i<n;i++)
		for(j=i+1;j<n;j++){
			temp = R[i][j];
			R[i][j] = R[j][i];
			R[j][i] = temp;
		}
}


/*-----------------------------------------------------------
 * Find the maximum off-diagonal entry, return the value and
 * the indices.
 *    A: the mtx,
 *    p, q: the indices,
 *    n: matrix dimension.
 */
double max_off_diag_entry(double **A, int *p, int *q, int n)
{
	int  i, j;
	double max;

	max = A[0][1];
	*p = 0;
	*q = 1;
	for(i=0;i<n;i++)
		for(j=i+1;j<n;j++)
			if(fabs(A[i][j])>fabs(max)){
				*p = i;
				*q = j;
				max = A[i][j];
			}
	return(max);
}

/*-----------------------------------------------------
 * Procedure to allocate space for an n-dimensional
 * vector. Return the pointer to the vector.
 */
double *alloc_vec(int n)
{
	double *t;

	t = (double *) malloc(sizeof(double)*n);
	return(t);
}


/*--------------------------------------------------
 * Compute inner product <a, b>, where a and b are 
 * n-dimensional vectors.
 */
double  inner_product(double *a, double *b, int n)
{
	int   i;
	double sum;

	sum = 0.0;
	for(i=0;i<n;i++)
		sum += a[i]*b[i];
	return (sum);
}

/*---------------------------------------------------------
 * Compute a = A*b, where a and b are vectors and A an n x n
 * matrix.
 */
void mtx_vec_mult(double *a, double **A, double *b, int n)
{
	int  i, j;
	
	for(i=0;i<n;i++){
		a[i] = 0.0;
		for(j=0;j<n;j++)
			a[i] += A[i][j]*b[j];
	}
}


/*--------------------------------------------------------
 * Procedure to compute the 2-norm of a vector.
 */
double vec_norm(double *a, int n)
{
	double   sum;
	int      i;

	sum = 0.0;
	for(i=0;i<n;i++)
		sum += a[i]*a[i];
	sum = sqrt(sum);
	return(sum);
}


/*--------------------------------------------------------------------
 *Procedure to normalize a vector
 */
void normalize_vec(double *x, int n)
{
  double   norm;
  int         i;

  norm = vec_norm(x, n);
  if(norm==0) return;
  for(i=0;i<n;i++)
     x[i] = x[i]/norm;
 }


/*-------------------------------------------------------------
 * Procedure to compute residual vector  r[] = x[] - u*y[].
 *    u: the eigen value.
 */
void comp_residual(double *r, double *x, double u, double *y, int n)
{
	int   i;

	for(i=0;i<n;i++)
		r[i] = x[i] - u*y[i];
}

/* --------------------------------------------------
 * free matrix
 */
void free_mtx(double** B, int n)
{
    int i;
    for (i = n-1; i >= 0; i--)
        free(B[i]);
}


/* -----------------------------
 * generate a N*N symmetric matrix. 
 */ 
void gen_sym_mtx(double **A, int n){
    int i, j;
    srand(0);
    for(i = 0; i < n; i++){
        for(j = i; j < n; j++)
            A[i][j] = A[j][i] = (double)(rand()%20);
    }
}


