#ifndef _ODE_HH_
#define _ODE_HH_

#include <cstdlib>
#include <complex>

#include "velocity-def.hh"

#define NEWTON_ITER_MAX 1000
#define NEWTON_LINE_SEARCH_ITER_MAX 100

int prop_implicit_euler_in_sph_harm_basis(
    const int N_s, const int N_rho, const int N_lm,
    const std::complex<double> **psi_t_next_arr_arr,
    const double *rho_arr, const int *l_arr, const int *m_arr,
    const double *rho_p_lim, double delta_t, double thres, 
    double r_p_t_vec[DIM_R], bool verbose=false);

#endif // _ODE_HH_

