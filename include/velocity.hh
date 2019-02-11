#ifndef _VELOCITY_HH_
#define _VELOCITY_HH_

#include <cstdio>
#include <cstdlib>
#include <complex>
#include <cstring>

#ifndef DEBUG
#define NDEBUG
#endif // DEBUG
#include <cassert>

#include <cmath>

#define CODE_OUT_OF_RANGE 256

template <typename T>
int eval_psi_and_dpsidx_p(
    const double x_p, T *psi_arr, const double *x_arr, 
    const int N_s, const int N_x, const double *x_p_lim, 
    T *psi_p_arr_p, T *dpsidx_p_arr_p);

template <typename T>
int eval_psi_and_dpsidx_arr(
    const double *x_p_arr, T *psi_arr, const double *x_arr, 
    const int N_s, const int N_p, const int N_x, const double *x_p_lim, 
    T *psi_p_arr, T *dpsidx_p_arr);

int eval_v_p_arr_for_sph_harm_basis(
    const int N_s, const int N_p, const int N_r_dim, 
    const int N_rho, const int N_lm,
    double **r_p_arr, std::complex<double> **psi_in_sph_harm_basis_arr,
    double *rho_arr, int *l_arr, int *m_arr, const double *rho_p_lim, 
    double **v_p_arr);

#endif // _VELOCITY_HH_
