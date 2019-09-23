#include "mex.h"
#include <stdlib.h>
#include <math.h>

#define TOL 1e-5
#  define c_absval(x) (((x) < 0) ? -(x) : (x)) /**< absolute value */
// #define MIN(a,b,c) (a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c))

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

void mat_vec(const double *Ax, const long *Ap, const long *Ai, size_t n, const double *x, double *y) {
    /* NB: Assume A is square */
    size_t Anzmax = Ap[n];
    size_t k, i;
    for (k = 0; k < n; k++) {
        y[k] = 0.0;
    }
    for (k = 0; k < Anzmax; k++) {
        i = Ai[k];
        y[i] += Ax[k]*x[i];
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

void orthonormal_basis(double *V1, double *V2, double *V3, double *x1, double *x2, double *x3, size_t n) {
    vec_mult_scalar(x1, 1/norm(x1, n), V1, n);
    vec_add_scaled(x2, V1, V2, -vec_prod(x2, V1, n), n);
    vec_self_mult_scalar(V2, 1/norm(V2, n), n);
    vec_double_add_scaled(x3, V1, V2, V3, -vec_prod(x3, V1, n), -vec_prod(x3, V2, n), n);
    vec_self_mult_scalar(V3, 1/norm(V3, n), n);
}

void vec_t_mat(const double *x, const double *Ax, const long *Ap, const long *Ai, size_t n, double *y) {
    /* y = x^T A*/
    size_t k, i, p;
    for (k = 0; k < n; k++) {
        y[k] = 0.0;
        for (p = Ap[k]; p < Ap[k+1]; p++) {
            i = Ai[p];
            y[k] += x[i]*Ax[p];
        }
    }
}

void vec_mat_compression(const double* V1, const double* V2, const double* V3, 
                            const double* Ax, const long *Ap, const long* Ai, size_t n,
                            double *V1A, double *V2A, double *V3A, double B[3][3]){

    vec_t_mat(V1, Ax, Ap, Ai, n, V1A);
    vec_t_mat(V2, Ax, Ap, Ai, n, V2A);
    vec_t_mat(V3, Ax, Ap, Ai, n, V3A);
    B[0][0] = vec_prod(V1A, V1, n);
    B[0][1] = vec_prod(V1A, V2, n);
    B[0][2] = vec_prod(V1A, V3, n);
    B[1][0] = vec_prod(V2A, V1, n);
    B[1][1] = vec_prod(V2A, V2, n);
    B[1][2] = vec_prod(V2A, V3, n);
    B[2][0] = vec_prod(V3A, V1, n);
    B[2][1] = vec_prod(V3A, V2, n);
    B[2][2] = vec_prod(V3A, V3, n);
}

void QR_decompostion(const double B[3][3], double Q[3][3], double R[3][3]) {
    double inv_norm, vec_product, vec_product2;
    inv_norm = 1.0/sqrt(B[0][0]*B[0][0] + B[1][0]*B[1][0] + B[2][0]*B[2][0]);
    Q[0][0] = B[0][0]*inv_norm; Q[1][0] = B[1][0]*inv_norm; Q[2][0] = B[2][0]*inv_norm;

    vec_product = Q[0][0]*B[0][1] + Q[1][0]*B[1][1] + Q[2][0]*B[2][1];
    Q[0][1] = B[0][1] - vec_product*Q[0][0];
    Q[1][1] = B[1][1] - vec_product*Q[1][0];
    Q[2][1] = B[2][1] - vec_product*Q[2][0];
    inv_norm = 1.0/sqrt(Q[0][1]*Q[0][1] + Q[1][1]*Q[1][1] + Q[2][1]*Q[2][1]);
    Q[0][1] = Q[0][1]*inv_norm; Q[1][1] = Q[1][1]*inv_norm; Q[2][1] = Q[2][1]*inv_norm;


    vec_product = Q[0][0]*B[0][2] + Q[1][0]*B[1][2] + Q[2][0]*B[2][2];
    vec_product2 = Q[0][1]*B[0][2] + Q[1][1]*B[1][2] + Q[2][1]*B[2][2];
    Q[0][2] = B[0][2] - vec_product*Q[0][0] - vec_product2*Q[0][1];
    Q[1][2] = B[1][2] - vec_product*Q[1][0] - vec_product2*Q[1][1];
    Q[2][2] = B[2][2] - vec_product*Q[2][0] - vec_product2*Q[2][1];
    inv_norm = 1.0/sqrt(Q[0][2]*Q[0][2] + Q[1][2]*Q[1][2] + Q[2][2]*Q[2][2]);
    Q[0][2] = Q[0][2]*inv_norm; Q[1][2] = Q[1][2]*inv_norm; Q[2][2] = Q[2][2]*inv_norm;

    /* R = Q^T B */
    size_t i, j, k;
    for(i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            R[i][j] = 0;
            for (k = 0; k < 3; k++) {
                R[i][j] += Q[k][i]*B[k][j];
            }
        }
    }

    // vec_mult_scalar(x, 1/norm(x, n), V1, n);
    // vec_add_scaled(x2, V1, V2, -vec_prod(x2, V1, n), n);
    // vec_self_mult_scalar(V2, 1/norm(V2, n), n);
    // vec_double_add_scaled(x3, V1, V2, V3, -vec_prod(x3, V1, n), -vec_prod(x3, V2, n), n);
    // vec_self_mult_scalar(V3, 1/norm(V3, n), n);
}

void QR_algorithm(double B[3][3], double eig_vec_B[3][3], double *lambda_B) {
    double temp[3][3];
    double Q[3][3];
    double R[3][3];
    eig_vec_B[0][0] = 1;
    eig_vec_B[0][1] = 0;
    eig_vec_B[0][2] = 0;
    eig_vec_B[1][0] = 0;
    eig_vec_B[1][1] = 1;
    eig_vec_B[1][2] = 0;
    eig_vec_B[2][0] = 0;
    eig_vec_B[2][1] = 0;
    eig_vec_B[2][2] = 1;
    
    size_t iter, i, j, k;

    // for (int row = 0; row < 3; row++) {
    //         mexPrintf("B(%d,:) = [%.4f %.4f %.4f] \n", row+1, B[row][0], B[row][1], B[row][2]);
    // }

    for (iter = 0; iter < 10; iter++) {
        // for (int row = 0; row < 3; row++) {
        //     for (int col = 0; col < 3; col++) {
        //         mexPrintf("B[%d][%d] = %.4f\n", row+1, col+1, B[row][col]);
        //     }
        // }
        // for (int row = 0; row < 3; row++) {
        //     mexPrintf("B(%d,:) = [%.4f %.4f %.4f] \n", row+1, B[row][0], B[row][1], B[row][2]);
        // }
        QR_decompostion(B, Q, R);
        
        // for (int row = 0; row < 3; row++) {
        //     mexPrintf("Q(%d,:) = [%.4f %.4f %.4f] \n", row+1, Q[row][0], Q[row][1], Q[row][2]);
        // }
        // for (int row = 0; row < 3; row++) {
        //     mexPrintf("R(%d,:) = [%.4f %.4f %.4f] \n", row+1, R[row][0], R[row][1], R[row][2]);
        // }
        for(i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                B[i][j] = 0;
                temp[i][j] = 0; 
                for (k = 0; k < 3; k++) {
                    temp[i][j] += eig_vec_B[i][k] * Q[k][j];
                    B[i][j] += R[i][k] * Q[k][j];
                }
            }
        }
        for(i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                eig_vec_B[i][j] = temp[i][j];
            }
        }
        
        // for (int row = 0; row < 3; row++) {
        //     mexPrintf("B(%d,:) = %.4f %.4f %.4f \n", row+1, B[row][0], B[row][1], B[row][2]);
        // }
    }
    lambda_B[0] = B[0][0];
    lambda_B[1] = B[1][1];
    lambda_B[2] = B[2][2];

    // for (int col = 0; col < 3; col++) {
    //         mexPrintf("eig_vec_B(:,%d) = [%.4f; %.4f; %.4f] \n", col+1, eig_vec_B[0][col], eig_vec_B[1][col], eig_vec_B[2][col]);
    // }
    // mexPrintf("Lambda: %.4f, %.4f, %.4f\n", lambda_B[0], lambda_B[1], lambda_B[2]);


}

void lobcpg_cleanup(double *x, double *x_prev, double *Ax,double *w, double *V1, double *V2, double *V3, double *V1A, double *V2A, double *V3A){
    mxFree(x);
    mxFree(x_prev);
    mxFree(Ax);
    mxFree(w);
    mxFree(V1);
    mxFree(V2);
    mxFree(V3);
    mxFree(V1A);
    mxFree(V2A);
    mxFree(V3A);
}

void mexFunction(int nlhs, mxArray * plhs [], int nrhs, const mxArray * prhs []) {

    if (nrhs != 1) {
        mexErrMsgTxt("Wrong number of input arguments.");
    }

    const mxArray* A = prhs[0];
    size_t n, m;
    m = mxGetM(A);
    n = mxGetN(A);
    if (n != m) {
        mexErrMsgTxt("Input needs to be a square matrix.");
    }
    long *Ap = (long *) mxGetJc (A);
    long *Ai = (long *) mxGetIr (A);
    double *Aval = mxGetPr(A);

    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *lambda_min_result; 
    lambda_min_result = mxGetPr(plhs[0]);

    double *x = mxCalloc(n, sizeof(double));
    double *x_prev = mxCalloc(n, sizeof(double));

    size_t i;
    for (i = 0; i < n; i++) {
        x[i] = (double) rand()/RAND_MAX;
        x_prev[i] = (double) rand()/RAND_MAX;
    }
    double x_norm, x_prev_norm;
    x_norm = norm(x, n);
    x_prev_norm = norm(x_prev, n);
    for (i = 0; i < n; i++) {
        x[i] = x[i]/x_norm;
        x_prev[i] = x_prev[i]/x_prev_norm;
    }

    double *Ax = mxCalloc(n, sizeof(double));
    double *w = mxCalloc(n, sizeof(double));
    double *V1 = mxCalloc(n, sizeof(double));
    double *V2 = mxCalloc(n, sizeof(double));
    double *V3 = mxCalloc(n, sizeof(double));
    double *V1A = mxCalloc(n, sizeof(double));
    double *V2A = mxCalloc(n, sizeof(double));
    double *V3A = mxCalloc(n, sizeof(double));
    double lambda_min;    
    double B[3][3];
    double eig_vec_B[3][3];
    double lambda_B[3];
    double y[3];
    double norm_w, abs_w_k;

    size_t max_iter = 10000;
    size_t k;
    mat_vec(Aval, Ap, Ai, n, x, Ax);
    lambda_min = vec_prod(x, Ax, n)/vec_prod(x, x, n);
    vec_add_scaled(Ax, x, w, -lambda_min, n);
    for (i = 0; i < max_iter; i++) {
        
        orthonormal_basis(V1, V2, V3, x, w, x_prev, n);
        vec_mat_compression(V1, V2, V3, Aval, Ap, Ai, n, V1A, V2A, V3A, B);

        QR_algorithm(B, eig_vec_B, lambda_B);

        if (lambda_B[0] < lambda_B[1]) {
            if (lambda_B[0] < lambda_B[2]) {
                lambda_min = lambda_B[0];
                y[0] = eig_vec_B[0][0]; y[1] = eig_vec_B[1][0]; y[2] = eig_vec_B[2][0];
            } else {
                lambda_min = lambda_B[2];
                y[0] = eig_vec_B[0][2]; y[1] = eig_vec_B[1][2]; y[2] = eig_vec_B[2][2];
            }
        } else if (lambda_B[1] < lambda_B[2]) {
            lambda_min = lambda_B[1];
            y[0] = eig_vec_B[0][1]; y[1] = eig_vec_B[1][1]; y[2] = eig_vec_B[2][1];
        } else {
            lambda_min = lambda_B[2];
            y[0] = eig_vec_B[0][2]; y[1] = eig_vec_B[1][2]; y[2] = eig_vec_B[2][2];
        }

        norm_w = 0.0;
        for (k = 0; k < n; k++) {
            x_prev[k] = x[k];
            x[k] = y[0]*V1[k] + y[1]*V2[k] + y[2]*V3[k];
            Ax[k] = y[0]*V1A[k] + y[1]*V2A[k] + y[2]*V3A[k];
            w[k] = Ax[k] - lambda_min*x[k];
            // mexPrintf("x[k] = %.4f, Ax[k] = %.4f, w[k] = %.4f \n", x[k], Ax[k], w[k]);
            abs_w_k = c_absval(w[k]);
            norm_w = abs_w_k > norm_w ? abs_w_k : norm_w;
        }
        // mat_vec(Aval, Ap, Ai, n, x, Ax);
        // mexPrintf("norm_w: %.4f \n", norm_w);
        if (norm_w < TOL) {
            *lambda_min_result = lambda_min;
            mexPrintf("Iter: %lu\n", i);
            lobcpg_cleanup(x, x_prev, Ax, w, V1, V2, V3, V1A, V2A, V3A );
            return;
        }
    }
    *lambda_min_result = lambda_min;
    mexPrintf("Not converged!\n", i);
    lobcpg_cleanup(x, x_prev, Ax, w, V1, V2, V3, V1A, V2A, V3A );


}