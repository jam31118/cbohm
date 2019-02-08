#ifndef _VELOCITY_HH_
#define _VELOCITY_HH_

#include <cstdio>
#include <cstdlib>
#include <complex>

#ifndef DEBUG
#define NDEBUG
#endif // DEBUG
#include <cassert>

#include <cmath>

template <typename T>
int eval_psi_and_dpsidx_arr(
    double *x_p_arr, T *psi_arr, double *x_arr, 
    int N_s, int N_p, int N_x, const double *x_p_lim, 
    T *psi_p_arr, T *dpsidx_p_arr);

int eval_v_p_arr_for_sph_harm_basis(
    const int N_s, const int N_p, const int N_r_dim, 
    const int N_rho, const int N_lm,
    double **r_p_arr, std::complex<double> **psi_in_sph_harm_basis_arr,
    double *rho_arr, int *l_arr, int *m_arr, const double *rho_p_lim, 
    double **v_p_arr);

#endif // _VELOCITY_HH_
