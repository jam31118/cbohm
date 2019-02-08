#ifndef _LOG_HH_
#define _LOG_HH_

#include <stdio.h>
#include <stdlib.h>

int return_with_mesg(const char *mesg, int return_code=EXIT_FAILURE);

void print_bar();

void print_3vec_p_arr(double *vec_p_arr, const int N_p, const int N_r_dim);

#endif // _LOG_HH_

