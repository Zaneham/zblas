/*
 * zblas.h — z/OS Native BLAS Library
 *
 * C prototypes for calling HLASM BLAS kernels from XL C/C++.
 * All routines follow standard BLAS signatures (Fortran-style:
 * everything by reference).
 *
 * Link with: -lzblas
 *
 * Apache 2.0 — see LICENSE.
 */

#ifndef ZBLAS_H
#define ZBLAS_H

#ifdef __cplusplus
extern "C" {
#endif

/* ---- BLAS Level 1 (vector operations) ---- */

/* DAXPY: y = alpha*x + y */
void daxpy_(const int *n, const double *alpha,
            const double *x, const int *incx,
            double *y, const int *incy);

/* DDOT: dot product */
double ddot_(const int *n,
             const double *x, const int *incx,
             const double *y, const int *incy);

/* DSCAL: x = alpha*x */
void dscal_(const int *n, const double *alpha,
            double *x, const int *incx);

/* DNRM2: euclidean norm */
double dnrm2_(const int *n,
              const double *x, const int *incx);

/* DCOPY: y = x */
void dcopy_(const int *n,
            const double *x, const int *incx,
            double *y, const int *incy);

/* DSWAP: swap x and y */
void dswap_(const int *n,
            double *x, const int *incx,
            double *y, const int *incy);

/* IDAMAX: index of element with max abs value (1-based) */
int idamax_(const int *n,
            const double *x, const int *incx);

/* DASUM: sum of absolute values */
double dasum_(const int *n,
              const double *x, const int *incx);

/* ---- BLAS Level 2 (matrix-vector) ---- */

/* DGEMV: y = alpha*A*x + beta*y  (or A' variant) */
void dgemv_(const char *trans, const int *m, const int *n,
            const double *alpha, const double *a, const int *lda,
            const double *x, const int *incx,
            const double *beta, double *y, const int *incy);

/* DTRSV: solve A*x = b (triangular) */
void dtrsv_(const char *uplo, const char *trans, const char *diag,
            const int *n, const double *a, const int *lda,
            double *x, const int *incx);

/* DGER: A = alpha*x*y' + A (rank-1 update) */
void dger_(const int *m, const int *n,
           const double *alpha,
           const double *x, const int *incx,
           const double *y, const int *incy,
           double *a, const int *lda);

/* ---- BLAS Level 3 (matrix-matrix) ---- */

/* DGEMM: C = alpha*A*B + beta*C */
void dgemm_(const char *transa, const char *transb,
            const int *m, const int *n, const int *k,
            const double *alpha, const double *a, const int *lda,
            const double *b, const int *ldb,
            const double *beta, double *c, const int *ldc);

/* DTRSM: solve A*X = alpha*B (triangular, matrix RHS) */
void dtrsm_(const char *side, const char *uplo,
            const char *transa, const char *diag,
            const int *m, const int *n,
            const double *alpha, const double *a, const int *lda,
            double *b, const int *ldb);

/* DSYRK: C = alpha*A*A' + beta*C (symmetric rank-k update) */
void dsyrk_(const char *uplo, const char *trans,
            const int *n, const int *k,
            const double *alpha, const double *a, const int *lda,
            const double *beta, double *c, const int *ldc);

#ifdef __cplusplus
}
#endif

#endif /* ZBLAS_H */
