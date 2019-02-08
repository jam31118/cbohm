#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "../../include/velocity.hh"

typedef std::complex<double> m_t;

void print_bar() {
  printf("\n----------------------\n\n");
}

int main(int argc, char *argv[]) {

  const int N_s = 4, N_p = 4, N_rho = 20;
  const int N_r_dim = 3, N_l = 2, qprop_dim=44;
  int initial_m = N_l + 1; // initialized with forbidden m value
  int N_lm;

  if (qprop_dim == 34) { N_lm = N_l; }
  else if (qprop_dim == 44) { N_lm = N_l * N_l; }
  else { return EXIT_FAILURE; }
  
  int l_arr[N_lm], m_arr[N_lm];

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


  double **r_p_arr = new double*[N_r_dim];
  r_p_arr[0] = new double[N_r_dim*N_p];
  r_p_arr[1] = r_p_arr[0] + N_p;
  r_p_arr[2] = r_p_arr[0] + 2*N_p;
  r_p_arr[0][0] = 0.5; r_p_arr[0][1] = 1.0; r_p_arr[0][2] = 2.1; r_p_arr[0][3] = 2.4;
  r_p_arr[1][0] = 0.5 * M_PI; r_p_arr[1][1] = 0.5 * M_PI; r_p_arr[1][2] = 0.5 * M_PI; r_p_arr[1][3] = 0.5 * M_PI;
  r_p_arr[2][0] = 0.0; r_p_arr[2][0] = 0.0; r_p_arr[2][2] = 0.0; r_p_arr[2][3] = 0.0;
//  { 
//    { 0.5, 1, 2.1 }, 
//    { 0.5 * M_PI, 0.5 * M_PI, 0.5 * M_PI },
//    { 0.0, 0.0, 0.0 }
//  };
  //double v_p_arr[N_r_dim][N_p];
  double **v_p_arr = new double*[N_r_dim];
  v_p_arr[0] = new double[N_r_dim*N_p];
  v_p_arr[1] = v_p_arr[0] + N_p;
  v_p_arr[2] = v_p_arr[0] + 2*N_p;

  double v_p_ana_arr[N_r_dim][N_p];
//  m_t psi_arr[N_lm][N_rho];
  m_t **psi_arr = new m_t*[N_lm];
  psi_arr[0] = new m_t[N_lm*N_rho];
  for (int i_lm=0; i_lm<N_lm; i_lm++) {
    psi_arr[i_lm] = psi_arr[0] + i_lm * N_rho;
  }

  double rho_arr[N_rho];
//  m_t psi_p_arr[N_p], dpsidx_p_arr[N_p];
  const double rho_max = 5.0;
  const double delta_rho = (rho_max - 0) / (N_rho);
  double rho_lim[2] = { delta_rho, rho_max };

  print_bar();

  printf("%-5s= %-5d%s\n", "N_s", N_s, "the number of stencils");
  printf("%-5s= %-5d%s\n", "N_p", N_p, "the number of particles");
  printf("%-5s= %-5d%s\n", "N_rho", N_rho, "the number of grid points");
  printf("%-5s= %-5d%s\n", "N_lm", N_lm, "the number of spherical harmonics");
  printf("%-5s= %-5d%s\n", "N_r_dim", N_r_dim, 
      "the dimension of position vector space");

  print_bar();

  printf("%10s","rho_arr: \n");
  for (int i_rho = 0; i_rho < N_rho; i_rho++) {
    rho_arr[i_rho] = rho_lim[0] + i_rho * delta_rho;
    printf("%7.3f",rho_arr[i_rho]);
  } printf("\n");

  double rho;
  for (int i_rho = 0; i_rho < N_rho; i_rho++) {
    rho = rho_arr[i_rho];
    psi_arr[0][i_rho] = 0.0;
//    psi_arr[0][i_rho] = rho * 2.0 * exp(-rho);  // s-orbital
 //   printf("rho: %.4f / psi_arr[0][i_rho]: %.4f\n", rho, psi_arr[0][i_rho].real());
    psi_arr[1][i_rho] = 0.0;
    psi_arr[2][i_rho] = 0.0;
    psi_arr[3][i_rho] = rho * 1.0 / (2.0*sqrt(6)) * rho * exp(-rho/2.0);
//    psi_arr[3][i_rho] = 0.0;
  }

  printf("%10s","psi_arr: \n");
  for (int i_lm = 0; i_lm < N_lm; i_lm++) {
    for (int i_rho = 0; i_rho < N_rho; i_rho++) {
      printf("%7.3f",psi_arr[i_lm][i_rho].real());
//    printf("(%7.3f,%7.3f)",psi_arr[i_x].real(),psi_arr[i_x].imag());
    } printf("\n");
  } printf("\n");

  print_bar();

  rho_lim[1] = rho_arr[N_rho-1];
  const double rho_p_lim[2] = { 0, rho_lim[1] };


  int return_code = EXIT_FAILURE;

  printf("before main routine\n");

  return_code = eval_v_p_arr_for_sph_harm_basis(
      N_s, N_p, N_r_dim, 
      N_rho, N_lm,
      (double **) r_p_arr, (m_t **) psi_arr, 
      rho_arr, l_arr, m_arr, rho_p_lim, 
      (double **) v_p_arr);

  if (return_code != EXIT_SUCCESS) {
    fprintf(stderr, "[ERROR] Failed to run 'eval_psi_and_dpsidx_arr()'");
    return return_code;
  }

  printf("after main routine\n");


  //// Evaluate analytical results
//  m_t psi_p_ana_arr[N_p];
//  m_t dpsidx_p_ana_arr[N_p];
  for (int i_p = 0; i_p < N_p; i_p++) {
    v_p_ana_arr[0][i_p] = 0.0;
    v_p_ana_arr[1][i_p] = 0.0;
    v_p_ana_arr[2][i_p] = 1.0 / r_p_arr[0][i_p];
  }


  //// Print results
  
  // print column names
  printf("real part: \n");
  printf(
      "%5s%15s%45s\n", 
      "i_p", "v_p", "v_p_ana");

  // print columns
  for (int i_p = 0; i_p < N_p; i_p++) {
    printf("%5d%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f\n", 
        i_p, 
        v_p_arr[0][i_p], v_p_arr[1][i_p], v_p_arr[2][i_p],
        v_p_ana_arr[0][i_p], v_p_ana_arr[1][i_p], v_p_ana_arr[2][i_p]
        );
//    printf("%5d%15.5f%15.5f%15.5f%15.5f\n", 
//        i_p, psi_p_arr[i_p].real(), psi_p_ana_arr[i_p].real(), 
//        dpsidx_p_arr[i_p].real(), dpsidx_p_ana_arr[i_p].real());
  }

  print_bar();
//
//  printf("imaginary part: \n");
//  printf(
//      "%5s%15s%15s%15s%15s\n", 
//      "i_p", "psi-estim", "psi", "dpsidx-estim", "dpsidx");
//  for (int i_p = 0; i_p < N_p; i_p++) {
//    printf("%5d%15.5f%15.5f%15.5f%15.5f\n", 
//        i_p, psi_p_arr[i_p].imag(), psi_p_ana_arr[i_p].imag(), 
//        dpsidx_p_arr[i_p].imag(), dpsidx_p_ana_arr[i_p].imag());
//  }
//
//  print_bar();
//

  delete [] r_p_arr[0];
  delete [] r_p_arr;
  delete [] v_p_arr[0];
  delete [] v_p_arr;


  //// Return from this program
  return EXIT_SUCCESS;
}

