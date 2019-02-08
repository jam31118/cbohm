#include "../include/lm.hh"

int eval_N_lm(int qprop_dim, int N_l, int *N_lm) {
  if (qprop_dim == 34) { *N_lm = N_l; }
  else if (qprop_dim == 44) { *N_lm = N_l * N_l; }
  else { return EXIT_FAILURE; }
  return EXIT_SUCCESS;
}


int eval_l_m_arr(
    int qprop_dim, int N_l, int *l_arr, int *m_arr, int initial_m) 
{

  if (qprop_dim == 34) { 
    if (initial_m >= N_l) { return EXIT_FAILURE; }
    for (int i_l = 0; i_l < N_l; i_l++) {
      l_arr[i_l] = i_l; m_arr[i_l] = initial_m;
    } 
  }
  else if (qprop_dim == 44) {
    int _i_lm = 0;
    for (int i_l = 0; i_l < N_l; i_l++) {
      for (int i_m = -i_l; i_m < i_l+1; i_m++) {
        l_arr[_i_lm] = i_l; m_arr[_i_lm] = i_m;
        _i_lm++;
      }
    }
  }
  else { return EXIT_FAILURE; }

  return EXIT_SUCCESS;

} 
