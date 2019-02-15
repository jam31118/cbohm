#include "../include/log.hh"


int return_with_mesg(const char *mesg, int return_code) {
  fprintf(stderr,"[ERROR] %s\n", mesg);
  return return_code;
}


int return_with_debug_mesg_func(
    const char *mesg, const char *file_name, const int line_num, 
    const char *func_name, int return_code) 
{
  char debug_mesg[MAX_STR_LEN];
  snprintf(
      debug_mesg, MAX_STR_LEN, 
      "@ '%s:%d' in '%s()' : %s", file_name, line_num, func_name, mesg);
  return return_with_mesg(debug_mesg, return_code);
}


void print_bar() {
  printf("\n----------------------\n\n");
}


void print_3vec_p_arr(double *vec_p_arr, const int N_p, const int N_r_dim) {
  for (int i_p = 0; i_p < N_p; i_p++) {
    for (int i_r_dim = 0; i_r_dim < N_r_dim; i_r_dim++) {
      printf("%15.5f",vec_p_arr[i_r_dim*N_p+i_p]);
    } printf("\n");
  }
}


void print_jac_t(jac_t jac) {

    for (int i_row = 0; i_row < DIM_R; i_row++) {
      for (int i_col = 0; i_col < DIM_R; i_col++)
      { printf("%15.5f", jac[i_row][i_col]); }
      printf("\n");
    } 
  
}

void print_jac_t_complex(z_t z_jac[DIM_R][DIM_R]) {

  for (int i=0; i<DIM_R; i++) {
    for (int j=0; j<DIM_R; j++) {
      printf("(%15.10f,%15.10f)", z_jac[i][j].real(), z_jac[i][j].imag()) ;
    } printf("\n");
  } 

}

void print_vec_t_complex(vec_z_t vec_z) {
  for (int i=0; i<DIM_R; i++) { 
    printf("(%15.10f,%15.10f)", vec_z[i].real(), vec_z[i].imag()); 
  } printf("\n");
}

void print_vec_t(vec_t vec) {
  for (int i=0; i<DIM_R; i++) { 
    printf("%15.10f", vec[i]); 
  } printf("\n");
}
