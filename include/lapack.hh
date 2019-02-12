#ifndef _LAPACK_HH_
#define _LAPACK_HH_

#include <cstdio>
#include <cstdlib>

extern "C" {

extern int dgesv_(
    const int *n, const int *nrhs, double *A, const int *lda, 
    int *ipiv, double *B, const int *ldb, int *info);

// Function Argument(s)
//
//
// [INPUT]
//
// `n`: The number of linear equations, 
//      i.e., the order of the matrix A. `n` >= 0.
//
// `nrhs`: The number of right hand sides, 
//         i.e., the number of columns of the matrix `B`. `nrhs` >= 0.
//
// `lda`: The leading dimension of the array A.  `lda` >= max(1,`n`).
//
// `ldb`: The leading dimension of the array `B`.  `ldb` >= max(1,`n`).
// 
//
// [INPUT/OUTPUT]
//
// `A`: array of shape (`lda`,`n`)
//      n entry, the N-by-N coefficient matrix `A`.
//      On exit, the factors L and U from the factorization
//      A = P*L*U; the unit diagonal elements of L are not stored.
//
// `B`: 2D array of shape (`ldb`,`nrhs`)
//      On entry, the `n`-by-`nrhs` matrix of right hand side matrix B.
//      On exit, if `info` = 0, the N-by-NRHS solution matrix X.
//
//
// [OUTPUT]
//
// `ipiv`: array of shape (`n`,)
//         The pivot indices that define the permutation matrix P;
//         row i of the matrix was interchanged with row `ipiv[i]`.
//
// `info` = 0:  successful exit
//        < 0:  if INFO = -i, the i-th argument had an illegal value
//        > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
//              has been completed, but the factor U is exactly
//              singular, so the solution could not be computed.
//

}

int handle_gesv_info(int info);

#endif // _LAPACK_HH_
