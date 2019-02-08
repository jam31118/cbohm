#include "../include/hydrogen.hh"
#include "base/common.hh"
#include "../include/log.hh"

int propa_psi_arr(
    std::complex<double> **psi_arr, double delta_t, 
    const int N_rho, const int N_lm, const int qprop_dim, const int initial_m)
{

  long n, l, m;
  int return_code;
  double E_n, phase;
  std::complex<double> U;

  for (int i_lm = 0; i_lm < N_lm; i_lm++) {

    return_code = get_ell_and_m_from_lm_index(
        i_lm, &l, &m, initial_m, qprop_dim);
    if (return_code != 0) { 
      return return_with_mesg("Failed to evaluate l, m"); }

    n = l + 1; 
    E_n = -0.5 / (n*n);
    phase = -E_n * delta_t;
    U = std::complex<double>(cos(phase),sin(phase));

    for (int i_rho=0; i_rho<N_rho; i_rho++) {
      psi_arr[i_lm][i_rho] *= U;
    }

  }
  return EXIT_SUCCESS;
}

