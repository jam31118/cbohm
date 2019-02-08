#ifndef _LAPACK_HH_
#define _LAPACK_HH_

#include <cstdio>
#include <cstdlib>

extern "C" {

extern int dgesv_(
    int *n, int *nrhs, double *A, int *lda, int *ipiv, double *B, 
    int *ldb, int *info);

}

int handle_gesv_info(int info);

#endif // _LAPACK_HH_
