#ifndef _LM_HH_
#define _LM_HH_

#include <cstdlib>

int eval_N_lm(int qprop_dim, int N_l, int *N_lm);

int eval_l_m_arr(
    int qprop_dim, int N_l, int *l_arr, int *m_arr, int initial_m);

#endif // _LM_HH_
