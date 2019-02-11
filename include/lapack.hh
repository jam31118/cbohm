#ifndef _LAPACK_HH_
#define _LAPACK_HH_

#include <cstdio>
#include <cstdlib>

extern "C" {

extern int dgesv_(
    const int *n, const int *nrhs, double *A, const int *lda, 
    const int *ipiv, double *B, const int *ldb, int *info);

}

int handle_gesv_info(int info);

#endif // _LAPACK_HH_
