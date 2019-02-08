#ifndef _HYDROGEN_HH_
#define _HYDROGEN_HH_

#include <complex>

int propa_psi_arr(
    std::complex<double> **psi_arr, double delta_t, 
    const int N_rho, const int N_lm, const int qprop_dim, const int initial_m);

#endif // _HYDROGEN_HH_
