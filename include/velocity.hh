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

#include "velocity-def.hh"


double vec_norm(vec_t vec);

int move_to_canonical_range_of_sph_coord(vec_t r_vec);

template <typename T>
int eval_psi_deriv_p(
    const double x_p, const T *psi_arr, const double *x_arr, 
    const int N_s, const int N_x, const double *x_p_lim, const int N_o,
    const int *deriv_order_arr, T *psi_deriv_p_arr);


template <typename T>
int eval_psi_and_dpsidx_p(
    const double x_p, const T *psi_arr, const double *x_arr, 
    const int N_s, const int N_x, const double *x_p_lim, 
    T *psi_p_arr_p, T *dpsidx_p_arr_p);


template <typename T>
int eval_psi_and_dpsidx_arr(
    const double *x_p_arr, T *psi_arr, const double *x_arr, 
    const int N_s, const int N_p, const int N_x, const double *x_p_lim, 
    T *psi_p_arr, T *dpsidx_p_arr);

int eval_psi_and_di_psi_and_dj_di_psi_for_sph_harm_basis(
    const int N_s, const int N_rho, const int N_lm, 
    const double r_p_vec[DIM_R], 
    const std::complex<double> **psi_in_sph_harm_basis_arr,
    const double *rho_arr, const int *l_arr, const int *m_arr, 
    const double *rho_p_lim, 
    z_t *psi_p_p, z_t di_psi[DIM_R], z_t dj_di_psi[DIM_R][DIM_R]=NULL);

int eval_jac_v_3D(
    const jac_t dj_hi, const vec_t hi, const z_t dj_di_psi[DIM_R][DIM_R],
    const z_t di_psi[DIM_R], const z_t psi, jac_t dj_vi);

int eval_v_p_for_sph_harm_basis(
    const int N_s, const int N_rho, const int N_lm, 
    double r_p_vec[DIM_R], 
    const std::complex<double> **psi_in_sph_harm_basis_arr,
    const double *rho_arr, const int *l_arr, const int *m_arr, 
    const double *rho_p_lim, 
    double v_p_vec[DIM_R], double jac[DIM_R][DIM_R]=NULL);


int eval_v_p_vec_arr_for_sph_harm_basis(
    const int N_s, const int N_p, const int N_rho, const int N_lm,
    double **r_p_vec_arr, 
    const std::complex<double> **psi_in_sph_harm_basis_arr,
    const double *rho_arr, const int *l_arr, const int *m_arr, 
    const double *rho_p_lim, 
    double **v_p_vec_arr, jac_t *jac_p_arr=NULL);


int eval_v_p_arr_for_sph_harm_basis(
    const int N_s, const int N_p, const int N_rho, const int N_lm,
    double **r_p_arr, const std::complex<double> **psi_in_sph_harm_basis_arr,
    const double *rho_arr, const int *l_arr, const int *m_arr, 
    const double *rho_p_lim, 
    double **v_p_arr, jac_t *jac_p_arr=NULL);


#endif // _VELOCITY_HH_
