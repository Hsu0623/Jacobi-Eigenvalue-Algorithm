/*--------------------------------------------------------
 * Procedure to allocate an n X n matrix and return the 
 * pointer to the matrix.
 */
double **alloc_mtx(int n);

/*----------------------------------------------------------------------
 * Procedure to make an identity matrix.
 *   A: the input matrix,
 *   n: dimension of the matrix.
 */
void make_identity_mtx(double **A, int n);


/*--------------------------------------------------------------------
 * Post-multiply a matrix by another matrix.
 *      A=  A*R,
 * Input:
 *      A: the destination matrix,
 *      R: the mtx to multiply with A
 *      n: dimension of the matrix.
 */
void post_mtx_mult(double **A, double **R, int n);


/*--------------------------------------------------------------------
 * Pre-multiply a matrix by another matrix.
 *      A=  R*A,
 * Input:
 *      A: the destination matrix,
 *      R: the mtx to multiply with A
 *      n: dimension of the matrix.
 */
void pre_mtx_mult(double **R, double **A, int n);

/*-------------------------------------------------------------------
 * Procedure to form a rotational mtx,
 *     R: the input matrix,
 *     p, q: indinces of the cos() and sin() values,
 *     c, s: cos and sine values,
 *     n: dimension of the matrix.
 */
void make_rotate_mtx(double **R,  int p, int q, double c, double s, int n);

/*-------------------------------------------------------------------
 * Procedure to transpose a matrix.
 *   R: the input matrix.
 *   n: dimension of the matrix.
 */
void transpose_mtx(double **R, int n);


/*-----------------------------------------------------------
 * Find the maximum off-diagonal entry, return the value and
 * the indices.
 *    A: the mtx,
 *    p, q: the indices,
 *    n: matrix dimension.
 */
double max_off_diag_entry(double **A, int *p, int *q, int n);

/*--------------------------------------------------
 * Compute inner product <a, b>, where a and b are 
 * n-dimensional vectors.
 */
double  inner_product(double *a, double *b, int n);

/*-----------------------------------------------------
 * Procedure to allocate space for an n-dimensional
 * vector. Return the pointer to the vector.
 */
double *alloc_vec(int n);

/*---------------------------------------------------------
 * Compute a = A*b, where a and b are vectors and A an n x n
 * matrix.
 */
void mtx_vec_mult(double *a, double **A, double *b, int n);

/*--------------------------------------------------------
 * Procedure to compute the 2-norm of a vector.
 */
double vec_norm(double *a, int n);


/*--------------------------------------------------------------------
 *Procedure to normalize a vector
 */
void normalize_vec(double *x, int n);


/*-------------------------------------------------------------
 * Procedure to compute residual vector  r[] = x[] - u*y[].
 *    u: the eigen value.
 */
void comp_residual(double *r, double *x, double u, double *y, int n);


/* ------------------------------------------------------------
 * free matrix
 */
void free_mtx(double** B, int n);


/* -----------------------------
 * generate a N*N symmetric matrix. 
 */ 
void gen_sym_mtx(double **A, int n);


