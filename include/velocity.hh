#ifndef _VELOCITY_HH_
#define _VELOCITY_HH_

#include <cstdio>
#include <cstdlib>
#include <complex>

#ifndef DEBUG
#define NDEBUG
#endif // DEBUG
#include <cassert>

extern "C" {

extern int dgesv_(
    int *n, int *nrhs, double *A, int *lda, int *ipiv, double *B, 
    int *ldb, int *info);

}

int handle_gesv_info(int info);

int return_with_mesg(const char *mesg, int return_code=EXIT_FAILURE);

template <typename T>
int eval_psi_and_dpsidx_arr(
    double *x_p_arr, T *psi_arr, double *x_arr, 
    int N_s, int N_p, int N_x, const double *x_p_lim, 
    T *psi_p_arr, T *dpsidx_p_arr);

#endif // _VELOCITY_HH_
