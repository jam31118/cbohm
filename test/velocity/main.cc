#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "../../include/velocity.hh"

void print_bar() {
  printf("\n----------------------\n\n");
}

int main(int argc, char *argv[]) {

  const int N_s = 4, N_p = 3, N_x = 20;
  double x_p_arr[N_p] = { -0.1, 1, 2.1 };
  double psi_arr[N_x], x_arr[N_x];
  double psi_p_arr[N_p], dpsidx_p_arr[N_p];
  double x_lim[2] = { -3, 3 };
  const double delta_x = (x_lim[1] - x_lim[0]) / (N_x - 1);

  print_bar();

  printf("%-5s= %-5d%s\n", "N_s", N_s, "the number of stencils");
  printf("%-5s= %-5d%s\n", "N_p", N_p, "the number of particles");
  printf("%-5s= %-5d%s\n", "N_x", N_x, "the number of grid points");

  print_bar();

  printf("%10s","x_arr");
  for (int i_x = 0; i_x < N_x; i_x++) {
    x_arr[i_x] = x_lim[0] + i_x * delta_x;
    psi_arr[i_x] = exp(-x_arr[i_x] * x_arr[i_x]);
    printf("%7.3f",x_arr[i_x]);
  } printf("\n");
  printf("%10s","psi_arr");
  for (int i_x = 0; i_x < N_x; i_x++) {
    printf("%7.3f",psi_arr[i_x]);
  } printf("\n");

  print_bar();

  x_lim[1] = x_arr[N_x-1];
  const double x_p_lim[2] = { x_lim[0], x_lim[1] };


  int return_code = EXIT_FAILURE;

  return_code = eval_psi_and_dpsidx_arr<double>(
      x_p_arr, psi_arr, x_arr, 
      N_s, N_p, N_x, x_p_lim, 
      psi_p_arr, dpsidx_p_arr);

  if (return_code != EXIT_SUCCESS) {
    fprintf(stderr, "[ERROR] Failed to run 'eval_psi_and_dpsidx_arr()'");
    return return_code;
  }


  //// Evaluate analytical results
  double psi_p_ana_arr[N_p];
  double dpsidx_p_ana_arr[N_p];
  double x_p;
  for (int i_p = 0; i_p < N_p; i_p++) {
    x_p = x_p_arr[i_p];
    psi_p_ana_arr[i_p] = exp(-x_p*x_p);
    dpsidx_p_ana_arr[i_p] = -2.0 * x_p * exp(-x_p*x_p);
  }


  //// Print results
  
  // print column names
  printf(
      "%5s%15s%15s%15s%15s\n", 
      "i_p", "psi-estim", "psi", "dpsidx-estim", "dpsidx");

  // print columns
  for (int i_p = 0; i_p < N_p; i_p++) {
    printf("%5d%15.5f%15.5f%15.5f%15.5f\n", 
        i_p, psi_p_arr[i_p], psi_p_ana_arr[i_p], 
        dpsidx_p_arr[i_p], dpsidx_p_ana_arr[i_p]);
  }

  print_bar();


  //// Return from this program
  return EXIT_SUCCESS;
}

