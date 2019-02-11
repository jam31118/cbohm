#ifndef _LOG_HH_
#define _LOG_HH_

#include <stdio.h>
#include <stdlib.h>

#define MAX_STR_LEN 1024

int return_with_mesg(const char *mesg, int return_code=EXIT_FAILURE);

int return_with_debug_mesg_func(
    const char *mesg, const char *file_name, const int line_num, 
    const char *func_name, int return_code=EXIT_FAILURE);

#define debug_mesg(mesg) return_with_debug_mesg_func(\
    mesg, __FILE__, __LINE__, __func__)
#define debug_mesg_with_code(mesg, return_code) return_with_debug_mesg_func(\
    mesg, __FILE__, __LINE__, __func__, return_code)

void print_bar();

void print_3vec_p_arr(double *vec_p_arr, const int N_p, const int N_r_dim);

#endif // _LOG_HH_

