/* Copyright © 2019 Ben Hermans */
 
#include "mex.h"
#include <stdlib.h>
#include <math.h>
#include "lapack.h"

#define TOL 1e-5
#  define c_absval(x) (((x) < 0) ? -(x) : (x)) /* absolute value */

double vec_prod(const double *a, const double *b, size_t n) {
  double prod = 0.0;
  size_t i = 0; 

  if(n >= 4) {
      for (; i <= n-4; i+=4) {
        prod += (a[i]*b[i] + a[i+1]*b[i+1] + a[i+2]*b[i+2] + a[i+3]*b[i+3]);
      }
  }
  for (; i < n; i++) {
    prod += a[i] * b[i];
  }

  return prod;
}

double norm(const double *a, size_t n) {
    return sqrt(vec_prod(a, a, n));
}

double vec_norm_inf(const double *a, size_t n) {
    register size_t j = 0;
    register double s0 = 0.;
    register double s1 = 0.;
    register double s2 = 0.;
    register double s3 = 0.;
    register double max0 = 0.;
    register double max1 = 0.;
    register double max2 = 0.;
    register double max3 = 0.;
    const size_t block_size = 4;
    const size_t block_len = n >> 2;
    const size_t remaining = n % block_size; /*Initializing four blocks for 
                                                * efficient implementation of 
                                                * inner product using registers */

    while (j < block_len * block_size) {
      s0 = c_absval(a[j]); max0 = s0 > max0 ? s0 : max0;
      s1 = c_absval(a[j+1]); max1 = s1 > max1 ? s1 : max1;
      s2 = c_absval(a[j+2]); max2 = s2 > max2 ? s2 : max2;
      s3 = c_absval(a[j+3]); max3 = s3 > max3 ? s3 : max3;
      j+=4;
    }    

    max0 = max0 > max1 ? max0 : max1;
    max0 = max0 > max2 ? max0 : max2;
    max0 = max0 > max3 ? max0 : max3;
    j = block_size * block_len;
    switch (remaining) {
        case 3: max0 = max0 > c_absval(a[j+2]) ? max0 : c_absval(a[j+2]);
        case 2: max0 = max0 > c_absval(a[j+1]) ? max0 : c_absval(a[j+1]);
        case 1: max0 = max0 > c_absval(a[j+0]) ? max0 : c_absval(a[j]); /*Taking contribution from the last terms
                                    * that were not included in the block*/
        case 0:;
    }
    return max0;

}

void mat_vec(const double *Ax, const long *Ap, const long *Ai, size_t n, const double *x, double *y) {
    /* NB: Assume A is square */
    size_t Anzmax = Ap[n];
    size_t k, i, p;
    for (k = 0; k < n; k++) {
        y[k] = 0.0;
    }
    for (k = 0; k < n; k++) {
        for (p = Ap[k]; p < Ap[k+1]; p++) {
            i = Ai[p];
            y[i] += Ax[p]*x[k];
        }
    }

}

void vec_add_scaled(const double *a, const double *b, double *c, double sc, size_t n) {
  size_t i;

  for (i = 0; i < n; i++) {
    c[i] =  a[i] + sc * b[i];
  }
}

void vec_double_add_scaled(const double *a, const double *b1, const double *b2, double *c, double sc1, double sc2, size_t n) {
  size_t i;

  for (i = 0; i < n; i++) {
    c[i] =  a[i] + sc1*b1[i] + sc2*b2[i];
  }
}

void vec_mult_add_scaled(double *a, const double *b, double sc1, double sc2, size_t n) {
  size_t i;

  for (i = 0; i < n; i++) {
    a[i] = sc1*a[i] + sc2*b[i];
  }
}


void vec_mult_scalar(const double *a, double sc, double *b, size_t n) {
  size_t i;

  for (i = 0; i < n; i++) {
    b[i] = sc*a[i];
  }
}

void vec_self_mult_scalar(double *a, double sc, size_t n) {
  size_t i;

  for (i = 0; i < n; i++) {
    a[i] *= sc;
  }
}


void lobcpg_cleanup(double *x, double *Ax, double *w, double *Aw, double *p, double *Ap){
    mxFree(x);
    mxFree(Ax);
    mxFree(w);
    mxFree(Aw);
    mxFree(p);
    mxFree(Ap);
}

void mexFunction(int nlhs, mxArray * plhs [], int nrhs, const mxArray * prhs []) {

    if (nrhs != 1) {
        mexErrMsgTxt("Wrong number of input arguments.");
    }

    /* Get the matrix from mex */
    const mxArray* A = prhs[0];
    size_t n, m;
    m = mxGetM(A);
    n = mxGetN(A);
    if (n != m) {
        mexErrMsgTxt("Input needs to be a square matrix.");
    }
    long *Acol = (long *) mxGetJc (A);
    long *Arow = (long *) mxGetIr (A);
    double *Aval = mxGetPr(A);

    /* Output pointer */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *lambda_min_result; 
    lambda_min_result = mxGetPr(plhs[0]);

    /* Malloc for the x, w and p vectors, and their A* counterparts */
    double *x = mxCalloc(n, sizeof(double)); /*Current guess of the eigenvector */
    double *Ax = mxCalloc(n, sizeof(double));
    double *w = mxCalloc(n, sizeof(double)); /*Current residual, Ax - lambda*x */
    double *Aw = mxCalloc(n, sizeof(double));
    double *p = mxCalloc(n, sizeof(double)); /* Conjugate gradient direction */
    double *Ap = mxCalloc(n, sizeof(double));
    double p_norm_inv;

    double lambda_min;    
    double B[3][3]; /*Compressed A on the [x, w, p] basis */
    double C[3][3] = {1, 0, 0, 0, 1, 0, 0, 0, 1}; /* Takes into account that p is not orthonormal with x, w */
    double lambda_B[3]; /* The eigenvalues of B */
    double *y; /* Eigenvector corresponding to min(lambda_B) */
    double xAw, wAw, xAp, wAp, pAp, xp, wp;

    /* Initialize eigenvector randomly */
    size_t i;
    for (i = 0; i < n; i++) {
        x[i] = (double) rand()/RAND_MAX;
    }
    vec_self_mult_scalar(x, 1.0/norm(x, n), n);
    mat_vec(Aval, Acol, Arow, n, x, Ax);
    lambda_min = vec_prod(x, Ax, n);

    /* Compute residual and make it orthonormal wrt x */
    vec_add_scaled(Ax, x, w, -lambda_min, n);
    vec_add_scaled(w, x, w, -vec_prod(x, w, n), n);
    vec_self_mult_scalar(w, 1.0/norm(w, n), n);
    mat_vec(Aval, Acol, Arow, n, w, Aw);
    xAw = vec_prod(Aw, x, n);
    wAw = vec_prod(Aw, w, n);

    /* In the first compressed system, there is no p yet, so it is 2 by 2 */
    double B_init[2][2] = {lambda_min, xAw, xAw, wAw};
    double lambda_init[2];

    /* Lapack variables */
    long int info = 0, dim = 2, lwork = 10, itype = 1;
    double work[10];
    char jobz = 'V';
    char uplo = 'L';
 
    /* Solve eigenvalue problem */
    dsyev(&jobz, &uplo, &dim, *B_init, &dim, lambda_init, work, &lwork, &info);
    lambda_min = lambda_init[0];
    y = B_init[0];

    /* Compute first p */
    vec_mult_scalar(w, y[1], p, n);
    vec_mult_scalar(Aw, y[1], Ap, n);
    vec_add_scaled(p, x, x, y[0], n);
    vec_add_scaled(Ap, Ax, Ax, y[0], n);
    
    dim = 3; /* From now on, the dimension of the eigenproblem to solve will be 3 */
    size_t max_iter = 1000;
    for (i = 0; i < max_iter; i++) {

        /* Update w */
        vec_add_scaled(Ax, x, w, -lambda_min, n);
        if (vec_norm_inf(w, n) < TOL) {
            *lambda_min_result = lambda_min;
            lobcpg_cleanup(x, Ax, w, Aw, p, Ap);
            return;
        } 
        vec_add_scaled(w, x, w, -vec_prod(x, w, n), n);
        vec_self_mult_scalar(w, 1.0/norm(w, n), n);
        mat_vec(Aval, Acol, Arow, n, w, Aw);
        xAw = vec_prod(Ax, w, n);
        wAw = vec_prod(w, Aw, n);

        /* Normalize p */
        p_norm_inv = 1.0/norm(p, n);
        vec_self_mult_scalar(p, p_norm_inv, n);
        vec_self_mult_scalar(Ap, p_norm_inv, n);

        /* Compress the system */
        xAp = vec_prod(Ax, p, n);
        wAp = vec_prod(Aw, p, n);
        pAp = vec_prod(Ap, p, n);
        xp = vec_prod(x, p, n);
        wp = vec_prod(w, p, n);

        B[0][0] = lambda_min; B[0][1] = xAw; B[0][2] = xAp; 
        B[1][0] = xAw; B[1][1] = wAw; B[1][2] = wAp; 
        B[2][0] = xAp; B[2][1] = wAp; B[2][2] = pAp;

        C[0][2] = xp; C[1][2] = wp; C[2][0] = xp; C[2][1] = wp; 
        C[2][2] = 1.0; /* The dsygv routine might override this element, therefore we reset it here.*/

        /* Solve eigenproblem B*x = lambda*C*x */
        dsygv(&itype, &jobz, &uplo, &dim, *B, &dim, *C, &dim, lambda_B, work, &lwork, &info);
        lambda_min = lambda_B[0];
        y = B[0];

        /* Update p and x */
        vec_mult_add_scaled(p, w, y[2], y[1], n);
        vec_mult_add_scaled(Ap, Aw, y[2], y[1], n);
        vec_mult_add_scaled(x, p, y[0], 1, n);
        vec_mult_add_scaled(Ax, Ap, y[0], 1, n);

    }

    *lambda_min_result = lambda_min;
    mexPrintf("Not converged!\n");
    lobcpg_cleanup(x, Ax, w, Aw, p, Ap);
    return;
}