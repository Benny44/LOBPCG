# LOBPCG
This code uses LOBPCG to find the minimum eigenvalue of a matrix. Knowledge of the minimum eigenvalue can be used to regularize the matrix to be positive (semi)definite.

## Install
Add LOBPCG to your matlab folder, and then run 
```
lobpcg_setup.m
```
## Usage
Using LOBPCG is very simple. For a sparse symmetrical matrix A, do 
```
lambda = lobpcg_mex(A);
```
If A is dense, do 
```
lambda = lobpcg_mex(sparse(A));
```

## References
The algorithm implemented here is described in detail in Algorithm 4.1 in 
``` 
Knyazev, A. V. (2001). Toward the optimal preconditioned eigensolver: Locally optimal block preconditioned conjugate gradient method. SIAM journal on scientific computing, 23(2), 517-541.
```
The method also makes use of the eigenvalue problem solving implemented in LAPACK, particularly of [`dsyev`](http://www.netlib.org/lapack/explore-html/d2/d8a/group__double_s_yeigen_ga442c43fca5493590f8f26cf42fed4044.html) and [`dsygv`](http://www.netlib.org/lapack/explore-html/d2/d8a/group__double_s_yeigen_ga007d33bcdcc697e17c6d15432f159b73.html).