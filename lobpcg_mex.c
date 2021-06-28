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

void prea_vec_copy(const double *a, double *b, size_t n) {
  size_t i;

  for (i = 0; i < n; i++) {
    b[i] = a[i];
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

#define c_min(a, b) (((a) < (b)) ? (a) : (b)) /**< minimum of two values */
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

double min_root_third_order(double a, double b, double c, double d)
{
    double r[3] = {0};
    double di, di_sqrt;
    if (a == 0)
    {
        // Not a cubic polynomial, should not happen 
        // mexPrintf("Error: Not a cubic polynomial. This should not happen\n");
    }
    else if (d == 0)
    {
        di = b*b - 4*a*c;
        if (d < 0)
        {
            mexPrintf("Error: Imaginary roots. This should not happen\n");
        }
        di_sqrt = sqrt(di);
        r[0] = (-b-di_sqrt)/(2*a);
        // r[1] = (-b+d_sqrt)/(2*a); //Not relevant since we want the min root
    }
    else
    {
        double temp, q, p, re, an, r13;
        temp = 1/a;
        b = b*temp;
        c = c*temp;
        d = d*temp;
        q = (3.0*c - (b*b))/9.0;
        p = (-(27.0*d) + b*(9.0*c - 2.0*(b*b)))/54.0;
        di = q*q*q + p*p;
        re = b/3.0;
        if (di > 0)
            mexPrintf("Imaginary roots, should not happen\n"); 
        else 
        {
            q = -q;
            an = q*q*q;
            an = acos(p/sqrt(an));
            r13 = 2.0*sqrt(q);
            r[0] = -re + r13*cos(an/3.0);
            r[1] = -re + r13*cos((an + 2.0*M_PI)/3.0);
            r[2] = -re + r13*cos((an + 4.0*M_PI)/3.0);
        }
    }
    
    if (r[0] <= r[1] && r[0] <= r[2]) return r[0];
    else return c_min(r[1], r[2]);
}

#include <math.h>
#  define c_sqrt sqrt /**< square root */

#define RREF_TOL 1e-8
# ifndef c_absval
#  define c_absval(x) (((x) < 0) ? -(x) : (x)) /**< absolute value */
# endif /* ifndef c_absval */

# ifndef c_max
#  define c_max(a, b) (((a) > (b)) ? (a) : (b)) /**< maximum of two values */
# endif /* ifndef c_max */

# ifndef mod
#  define mod(a,b) ((((a)%(b))+(b))%(b)) /**< modulo operation (positive result for all values) */
#endif

int custom_rref(double D[3][3])
{
    double p, k, temp[3], a[3];

    // First column
    a[0] = c_absval(D[0][0]); a[1] = c_absval(D[0][1]); a[2] = c_absval(D[0][2]); 
    if (a[0] < a[1] || a[0] < a[2])
    {
        if (a[1] > a[2])
        {
            if (a[1] < RREF_TOL) return 0;
            // swap row 0 and 1
            temp[0] = D[0][0]; temp[1] = D[0][1]; temp[2] = D[0][2];
            D[0][0] = D[1][0]; D[0][1] = D[1][1]; D[0][2] = D[1][2];
            D[1][0] = temp[0]; D[1][1] = temp[1]; D[1][2] = temp[2];  
        }
        else
        {
            if (a[2] < RREF_TOL) return 0;
            // swap row 0 and 2
            temp[0] = D[0][0]; temp[1] = D[0][1]; temp[2] = D[0][2];
            D[0][0] = D[2][0]; D[0][1] = D[2][1]; D[0][2] = D[2][2];
            D[2][0] = temp[0]; D[2][1] = temp[1]; D[2][2] = temp[2];  
        }
    }
    else
    {
        if (a[0] < RREF_TOL) return 0;
    }
    
    p = 1.0/D[0][0];
    D[0][1] *= p; D[0][2] *= p; D[0][0] = 1.0;
    D[1][1] -= D[1][0]*D[0][1]; D[1][2] -= D[1][0]*D[0][2]; D[1][0] = 0;
    D[2][1] -= D[2][0]*D[0][1]; D[2][2] -= D[2][0]*D[0][2]; D[2][0] = 0;

    // Second column
    a[1] = c_absval(D[1][1]); a[2] = c_absval(D[2][1]);
    if (a[1] < a[2])
    {
        if (a[2] < RREF_TOL) return 1;
        temp[2] = D[1][2];
        D[1][1] = D[2][1]; D[1][2] = D[2][2];
        D[2][2] = temp[2]; D[2][1] = 0;
    }
    else
    {
        if (a[1] < RREF_TOL) return 1;
    }

    p = 1.0/D[1][1];
    D[1][2] *= p; D[1][1] = 1.0;
    D[0][2] -= D[0][1]*D[1][2]; D[0][1] = 0;
    D[2][2] -= D[2][1]*D[1][2]; D[2][1] = 0;

    return 2;
}

double custom_eig(const double B[3][3], const double C[3][3], double x[3])
{
    double a, b, c, d;
    // TODO: this unpacking is silly
    double xqx = B[0][0], xqw = B[0][1], xqp = B[0][2], wqw = B[1][1], wqp = B[1][2], pqp = B[2][2], xp = C[0][2], wp = C[1][2];
    a = wp*wp + xp*xp - 1;
    b = (-xqx*wp*wp + 2*xqw*wp*xp - 2*wqp*wp - wqw*xp*xp - 2*xqp*xp + pqp + wqw + xqx);
    c = (wqp*wqp - 2*xp*wqp*xqw + 2*wp*xqx*wqp + xqp*xqp - 2*wp*xqp*xqw + 2*wqw*xp*xqp + xqw*xqw - pqp*wqw - pqp*xqx - wqw*xqx);
    d = - xqx*wqp*wqp + 2*wqp*xqp*xqw - wqw*xqp*xqp - pqp*xqw*xqw + pqp*wqw*xqx;
    double lam = min_root_third_order(a, b, c, d);
    double D[3][3];
    D[0][0] = B[0][0] - lam*C[0][0];
    D[0][1] = B[0][1];
    D[0][2] = B[0][2] - lam*C[0][2];
    D[1][0] = B[1][0];
    D[1][1] = B[1][1] - lam*C[1][1];
    D[1][2] = B[1][2] - lam*C[1][2];
    D[2][0] = B[2][0] - lam*C[2][0];
    D[2][1] = B[2][1] - lam*C[2][1];
    D[2][2] = B[2][2] - lam*C[2][2];
    int ind = custom_rref(D);

    if (ind == 0)
    {
        x[0] = 1; x[1] = 0; x[2] = 0;
    }
    else 
    if (ind == 1)
    {
        x[0] = 0; x[1] = 1; x[2] = 0;
    }
    else
    {
        double temp = 1/c_sqrt(1 + D[0][2]*D[0][2] - 2*D[0][2]*C[0][2] + D[1][2]*D[1][2] - 2*D[1][2]*C[1][2]);

        x[0] = -D[0][2]*temp;
        x[1] = -D[1][2]*temp;
        x[2] = temp;
    }

    return lam;
    
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
    
    double xAw, wAw, xAp, wAp, pAp, xp, wp;

    double *iter; 
    if (nlhs > 1)
    {
        /* Output pointer */
        plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
        iter = mxGetPr(plhs[1]);
    }

    double *eig_vec;
    if (nlhs > 2)
    {
        /* Output pointer */
        plhs[2] = mxCreateDoubleMatrix(n, 1, mxREAL);
        eig_vec = mxGetPr(plhs[2]);
    }

    double *lam_vec;
    if (nlhs > 3)
    {
        /* Output pointer */
        plhs[3] = mxCreateDoubleMatrix(10001, 1, mxREAL);
        lam_vec = mxGetPr(plhs[3]);
    }

    

    /* Initialize eigenvector randomly */
    // srand(2);
    size_t i;
    for (i = 0; i < n; i++) {
        x[i] = (double) rand()/RAND_MAX;
        // x[i] = 1.0/(i+1);
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

    #ifdef DLAPACK
    double *y; /* Eigenvector corresponding to min(lambda_B) */
    /* Lapack variables */
    long int info = 0, dim = 2, lwork = 10, itype = 1;
    double work[10];
    char jobz = 'V';
    char uplo = 'L';
 
    /* Solve eigenvalue problem */
    dsyev(&jobz, &uplo, &dim, *B_init, &dim, lambda_init, work, &lwork, &info);
    lambda_min = lambda_init[0];
    y = B_init[0];
    #else
    double y[3];
    double b, c, di;
    b = -(lambda_min + wAw);
    c = lambda_min*wAw - xAw*xAw;
    di = b*b - 4*c;
    lambda_min = (-b-c_sqrt(di))/2;
    B_init[0][0] -= lambda_min;
    B_init[1][1] -= lambda_min;
    if (c_absval(B_init[0][0]) < RREF_TOL)
    {
        y[0] = 1; y[1] = 0;
    }
    else
    {
        B_init[0][1] /= B_init[0][0];
        b = 1/c_sqrt(1 + B_init[0][1]*B_init[0][1]);
        y[0] = -B_init[0][1]*b;
        y[1] = b;
    }
    
    
    #endif
    // mexPrintf("Iter: %d, Lam = %e, y = [%e, %e]\n", 0, lambda_min, y[0], y[1]);

    if (nlhs > 3)
    {
        lam_vec[0] = lambda_min;
    }

    /* Compute first p */
    vec_mult_scalar(w, y[1], p, n);
    vec_mult_scalar(Aw, y[1], Ap, n);
    vec_add_scaled(p, x, x, y[0], n);
    vec_add_scaled(Ap, Ax, Ax, y[0], n);
    
    #ifdef DLAPACK
    dim = 3; /* From now on, the dimension of the eigenproblem to solve will be 3 */
    #endif
    size_t max_iter = 10000;
    for (i = 0; i < max_iter; i++) {

        /* Update w */
        vec_add_scaled(Ax, x, w, -lambda_min, n);
        if (vec_norm_inf(w, n) < TOL) {
            *lambda_min_result = lambda_min - c_sqrt(2.0)*norm(w, n);
            if (nlhs > 1)
            {
                *iter = (double) ++i;
            }
            if (nlhs > 2)
            {
                prea_vec_copy(x, eig_vec, n);
            }
            lobcpg_cleanup(x, Ax, w, Aw, p, Ap);

            // mexPrintf("Iter = %d\n", i);
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

        #ifdef DLAPACK
        /* Solve eigenproblem B*x = lambda*C*x */
        dsygv(&itype, &jobz, &uplo, &dim, *B, &dim, *C, &dim, lambda_B, work, &lwork, &info);
        lambda_min = lambda_B[0];
        y = B[0];
        #else 
        lambda_min = custom_eig(B, C, y);
        
        #endif
        
        if (nlhs > 3)
        {
            lam_vec[i+1] = lambda_min;
        }
        
        /* Update p and x */
        vec_mult_add_scaled(p, w, y[2], y[1], n);
        vec_mult_add_scaled(Ap, Aw, y[2], y[1], n);
        vec_mult_add_scaled(x, p, y[0], 1, n);
        vec_mult_add_scaled(Ax, Ap, y[0], 1, n);

        #ifndef DLAPACK
        if (mod(i, 50) == 0)
        {
            lambda_min = vec_prod(x, Ax, n);
        }
        #endif

        // mexPrintf("Iter: %d, Lam = %e, y = [%e, %e, %e]\n", i+1, lambda_min, y[0], y[1], y[2]);


    }

    *lambda_min_result = lambda_min;
    if (nlhs > 1)
    {
        *iter = (double) max_iter;
    }
    mexPrintf("Not converged!\n");
    lobcpg_cleanup(x, Ax, w, Aw, p, Ap);
    return;
}