#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "../../include/velocity.hh"
#include "../../include/log.hh"
#include "../../include/lm.hh"

typedef std::complex<double> m_t;


void print_bar() {
  printf("\n----------------------\n\n");
}

int main(int argc, char *argv[]) {

  //// Configuration
  const int N_s = 4, N_p = 4, N_rho = 20;
  const int N_r_dim = 3, N_l = 2, qprop_dim=44;
  int initial_m = N_l + 1; // initialized with forbidden m value
  const double rho_max = 5.0;


  int N_lm;
  if ( eval_N_lm(qprop_dim, N_l, &N_lm) != EXIT_SUCCESS) 
  { return return_with_mesg("Failed to evaluate 'N_lm'"); }

  int l_arr[N_lm], m_arr[N_lm];
  if ( eval_l_m_arr(qprop_dim, N_l, l_arr, m_arr, initial_m) != EXIT_SUCCESS )
  { return return_with_mesg("Failed to evaluate 'l_arr', 'm_arr'"); }



  //// Allocate memory for data

  // `r_p_arr`
  double **r_p_vec_arr = new double*[N_p];
  double *r_p_vec_arr_1d = new double[N_p*DIM_R];
//  for (int i_r_dim = 0; i_r_dim < N_r_dim; i_r_dim++) {
  for (int i_p = 0; i_p < N_p; i_p++) {
    r_p_vec_arr[i_p] = r_p_vec_arr_1d + i_p * DIM_R;
  }

  // `v_p_arr`
  double **v_p_vec_arr = new double*[N_p];
  double *v_p_vec_arr_1d = new double[N_p*DIM_R];
//  for (int i_r_dim = 0; i_r_dim < N_r_dim; i_r_dim++) {
  for (int i_p = 0; i_p < N_p; i_p++) {
    v_p_vec_arr[i_p] = v_p_vec_arr_1d + i_p * DIM_R;
  }

  // `v_p_ana_arr`
  double v_p_ana_arr[N_p][DIM_R];

  // `psi_arr`
  m_t **psi_arr = new m_t*[N_lm];
  psi_arr[0] = new m_t[N_lm*N_rho];
  for (int i_lm=0; i_lm<N_lm; i_lm++) {
    psi_arr[i_lm] = psi_arr[0] + i_lm * N_rho;
  }

  // array of jacobian
  jac_t *jac_p_arr = new jac_t[N_p];
  jac_t *jac_p_ana_arr = new jac_t[N_p];


  //// Initialization

  // `rho_p_arr`
  double rho_p_arr[N_p] = {0.5, 1.0, 2.1, 2.4};
  for (int i_p = 0; i_p < N_p; i_p++) {
    r_p_vec_arr[i_p][0] = rho_p_arr[i_p];
    r_p_vec_arr[i_p][1] = 0.3 * M_PI;
    r_p_vec_arr[i_p][2] = 0.0;
  }

  // Construct radial coordinate grid points `rho_arr`
  const double delta_rho = (rho_max - 0) / (N_rho);
  double rho_lim[2] = { delta_rho, rho_max };
  double rho_arr[N_rho];
  for (int i_rho = 0; i_rho < N_rho; i_rho++) {
    rho_arr[i_rho] = rho_lim[0] + i_rho * delta_rho;
  } 
  rho_lim[1] = rho_arr[N_rho-1];
  const double rho_p_lim[2] = { 0.0, rho_lim[1] };

  // Initialize the state funciton array `psi_arr`
  double rho;
  for (int i_rho = 0; i_rho < N_rho; i_rho++) {
    rho = rho_arr[i_rho];
    psi_arr[0][i_rho] = 0.0;
//    psi_arr[0][i_rho] = rho * 2.0 * exp(-rho);  // s-orbital
    psi_arr[1][i_rho] = 0.0;
    psi_arr[2][i_rho] = 0.0;
    psi_arr[3][i_rho] = 
      rho * 1.0 / (2.0*sqrt(6)) * rho * exp(-rho/2.0); // p(m=1) orbital
//    psi_arr[3][i_rho] = 0.0;
  }

  // Evaluate analytical results
  for (int i_p = 0; i_p < N_p; i_p++) {
    v_p_ana_arr[i_p][0] = 0.0;
    v_p_ana_arr[i_p][1] = 0.0;
    v_p_ana_arr[i_p][2] = 1.0 / (r_p_vec_arr[i_p][0] * sin(r_p_vec_arr[i_p][1]));
  }

  // Evaluate analytical jacobian
  for (int i_p = 0; i_p < N_p; i_p++) {
    for (int i_row = 0; i_row < DIM_R-1; i_row++) {
      for (int i_col = 0; i_col < DIM_R-1; i_col++) {
        jac_p_ana_arr[i_p][i_row][i_col] = 0.0;
      }
    }
    jac_p_ana_arr[i_p][2][0] = - 1.0 / r_p_vec_arr[i_p][0] * v_p_ana_arr[i_p][2];
    jac_p_ana_arr[i_p][2][1] = - cos(r_p_vec_arr[i_p][1]) / sin(r_p_vec_arr[i_p][1]) * v_p_ana_arr[i_p][2];
  }

  

  //// Evaluate velocity vector for each particle
  int return_code = EXIT_FAILURE;

//  return_code = eval_v_p_arr_for_sph_harm_basis(
//      N_s, N_p, N_rho, N_lm,
//      (double **) r_p_arr, (const m_t **) psi_arr, 
//      rho_arr, l_arr, m_arr, rho_p_lim, 
//      (double **) v_p_arr);

//  for (int i_p = 0; i_p < N_p; i_p++) {
  return_code = eval_v_p_vec_arr_for_sph_harm_basis(
      N_s, N_p, N_rho, N_lm,
      (double **) r_p_vec_arr, (const m_t **) psi_arr, 
      rho_arr, l_arr, m_arr, rho_p_lim, 
      (double **) v_p_vec_arr, jac_p_arr);
  if (return_code != EXIT_SUCCESS) {
    return debug_mesg_with_code(
        "Failed to run 'eval_psi_and_dpsidx_arr()'", return_code);
//    fprintf(stderr, "[ERROR] Failed to run 'eval_psi_and_dpsidx_arr()'");
//    return return_code;
  }

//  }



  //// Print results
  print_bar();

  printf("%-5s= %-5d%s\n", "N_s", N_s, "the number of stencils");
  printf("%-5s= %-5d%s\n", "N_p", N_p, "the number of particles");
  printf("%-5s= %-5d%s\n", "N_rho", N_rho, "the number of grid points");
  printf("%-5s= %-5d%s\n", "N_lm", N_lm, "the number of spherical harmonics");
  printf("%-5s= %-5d%s\n", "N_r_dim", N_r_dim, 
      "the dimension of position vector space");

  print_bar();

  printf("%7s","rho_arr: \n");
  for (int i_rho = 0; i_rho < N_rho; i_rho++) {
    printf("%7.3f",rho_arr[i_rho]);
  } printf("\n");

  printf("\n");

  printf("%10s","psi_arr (real): \n");
  for (int i_lm = 0; i_lm < N_lm; i_lm++) {
    for (int i_rho = 0; i_rho < N_rho; i_rho++) {
      printf("%7.3f",std::real(psi_arr[i_lm][i_rho]));
    } printf("\n");
  }

  printf("\n");

  printf("%10s","psi_arr (imag): \n");
  for (int i_lm = 0; i_lm < N_lm; i_lm++) {
    for (int i_rho = 0; i_rho < N_rho; i_rho++) {
      printf("%7.3f",std::imag(psi_arr[i_lm][i_rho]));
    } printf("\n");
  }

  print_bar();


  // Print 'v_p' and its analytical counterpart
  printf("real part: \n");
  printf(
      "%5s%15s%45s\n", 
      "i_p", "v_p", "v_p_ana");

  for (int i_p = 0; i_p < N_p; i_p++) {
    printf("%5d%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f\n", 
        i_p, 
        v_p_vec_arr[i_p][0], v_p_vec_arr[i_p][1], v_p_vec_arr[i_p][2],
        v_p_ana_arr[i_p][0], v_p_ana_arr[i_p][1], v_p_ana_arr[i_p][2]
        );
  }

  print_bar();
  
  // Print 'v_p' and its analytical counterpart
  printf("jacobian[N_p]: \n");
  printf(
      "%5s\t%-45s%-45s\n", 
      "i_p", "jac[i_p]", "jac_analytic[i_p]");

  for (int i_p = 0; i_p < N_p; i_p++) {
    printf("%5d\n",i_p);
    printf("     \t\n");
    for (int i_row = 0; i_row < DIM_R; i_row++) {
      for (int i_col = 0; i_col < DIM_R; i_col++)
      { printf("%15.5f", jac_p_arr[i_p][i_row][i_col]); }
      for (int i_col = 0; i_col < DIM_R; i_col++)
      { printf("%15.5f", jac_p_ana_arr[i_p][i_row][i_col]); }
      printf("\n");
    } printf("\n");
  }
  

  // Deallocate memory for arrays
  delete [] r_p_vec_arr[0];
  delete [] r_p_vec_arr;
  delete [] v_p_vec_arr[0];
  delete [] v_p_vec_arr;
  delete [] psi_arr[0];
  delete [] psi_arr;

  delete [] jac_p_arr;
  delete [] jac_p_ana_arr;


  //// Return from this program
  return EXIT_SUCCESS;
}

