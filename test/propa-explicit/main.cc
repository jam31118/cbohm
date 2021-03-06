#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "../../include/velocity.hh"
#include "../../include/log.hh"
#include "../../include/lm.hh"
#include "../../include/hydrogen.hh"
#include "../../include/ode.hh"

typedef std::complex<double> m_t;

int main(int argc, char *argv[]) {

  //// Configuration
  const int N_s = 4, N_p = 4, N_rho = 20;
  const int N_r_dim = 3, N_l = 2, qprop_dim=44;
  int initial_m = N_l + 1; // initialized with forbidden m value
  const double rho_max = 5.0;

  const int N_t = 4;
  const double delta_t = 0.2;
  const double t_0 = 0.2;

  double t_arr[N_t];
  for (int i_t = 0; i_t < N_t; i_t++) {
    t_arr[i_t] = t_0 + i_t * delta_t;
  }

  int N_lm;
  if ( eval_N_lm(qprop_dim, N_l, &N_lm) != EXIT_SUCCESS) 
  { return return_with_mesg("Failed to evaluate 'N_lm'"); }

  int l_arr[N_lm], m_arr[N_lm];
  if ( eval_l_m_arr(qprop_dim, N_l, l_arr, m_arr, initial_m) != EXIT_SUCCESS )
  { return return_with_mesg("Failed to evaluate 'l_arr', 'm_arr'"); }


  //// Define useful variables
  int return_code = EXIT_FAILURE;
  


  //// Allocate memory for data

  // Some variable
  const int r_p_arr_size = N_r_dim * N_p;

  // `r_p_arr`
  double **r_p_arr = new double*[N_r_dim];
  double *r_p_arr_1d = new double[r_p_arr_size];
  for (int i_r_dim = 0; i_r_dim < N_r_dim; i_r_dim++) {
    r_p_arr[i_r_dim] = r_p_arr_1d + i_r_dim * N_p;
  }

  // `v_p_arr`
  double **v_p_arr = new double*[N_r_dim];
  double *v_p_arr_1d = new double[r_p_arr_size];
  for (int i_r_dim = 0; i_r_dim < N_r_dim; i_r_dim++) {
    v_p_arr[i_r_dim] = v_p_arr_1d + i_r_dim * N_p;
  }

  // `v_p_ana_arr`
  double v_p_ana_arr[N_r_dim][N_p];

  // `psi_arr`
  m_t **psi_arr = new m_t*[N_lm];
  psi_arr[0] = new m_t[N_lm*N_rho];
  for (int i_lm=0; i_lm<N_lm; i_lm++) {
    psi_arr[i_lm] = psi_arr[0] + i_lm * N_rho;
  }

  // `r_p_t_arr`
  double **r_p_t_arr = new double*[N_t];
  double *r_p_t_arr_1d = new double[N_t*r_p_arr_size];
  for (int i_t = 0; i_t < N_t; i_t++) {
    r_p_t_arr[i_t] = r_p_t_arr_1d + i_t * r_p_arr_size;
  }
  
  double *r_p_t0_arr_1d = new double[r_p_arr_size];
  
  // `v_p_t_arr`
  double **v_p_t_arr = new double*[N_t];
  double *v_p_t_arr_1d = new double[N_t*r_p_arr_size];
  for (int i_t = 0; i_t < N_t; i_t++) {
    v_p_t_arr[i_t] = v_p_t_arr_1d + i_t * r_p_arr_size;
  }


  // Define variables for looping
  double *r_p_p;
  double *r_p_p_max = r_p_p_max = r_p_arr_1d + r_p_arr_size;


  //// Initialization

  // `rho_p_arr`
  double rho_p_arr[N_p] = {0.5, 1.0, 2.1, 2.4};
  for (int i_p = 0; i_p < N_p; i_p++) {
    r_p_arr[0][i_p] = rho_p_arr[i_p];
    r_p_arr[1][i_p] = 0.3 * M_PI;
    r_p_arr[2][i_p] = 0.0;
  }

  // `r_p_t0_arr`
  std::copy(r_p_arr_1d,r_p_arr_1d+r_p_arr_size,r_p_t0_arr_1d);

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
    v_p_ana_arr[0][i_p] = 0.0;
    v_p_ana_arr[1][i_p] = 0.0;
    v_p_ana_arr[2][i_p] = 1.0 / (r_p_arr[0][i_p] * sin(r_p_arr[1][i_p]));
  }

  // `v_p_arr`
 
  return_code = eval_v_p_arr_for_sph_harm_basis(
      N_s, N_p, N_rho, N_lm,
      r_p_arr, (const m_t **) psi_arr, 
      rho_arr, l_arr, m_arr, rho_p_lim, 
      v_p_arr);

  if (return_code != EXIT_SUCCESS) {
    fprintf(stderr, "[ERROR] Failed to run 'eval_psi_and_dpsidx_arr()'");
    return return_code;
  }

  // `r_p_t_arr_1d` and `v_p_t_arr_1d` at t=t_0
  std::copy(r_p_arr_1d,r_p_arr_1d+r_p_arr_size,r_p_t_arr_1d+0*r_p_arr_size);
  std::copy(v_p_arr_1d,v_p_arr_1d+r_p_arr_size,v_p_t_arr_1d+0*r_p_arr_size);
  


  const int N_st = 4;

  double A[N_st][N_st] = { {0,0,0,0},{0.5,0,0,0},{0,0.5,0,0,},{0,0,1,0} };
  const double c[N_st] = {0.0, 0.5, 0.5, 1.0};
  const double b[N_st] = {1.0/6.0,1.0/3.0,1.0/3.0,1.0/6.0};
  
  //// Allocate temporary buffer
  double ***k_st_arr = new double**[N_st];
  double **k_st_arr_2d = new double*[N_st*N_r_dim];
  double *k_st_arr_1d = new double[N_st*N_r_dim*N_p];
  for (int i_st = 0; i_st < N_st; i_st++) {
    k_st_arr[i_st] = k_st_arr_2d + i_st * N_r_dim;
    for (int i_st_r_dim = 0; i_st_r_dim < N_st*N_r_dim; i_st_r_dim++) {
      k_st_arr_2d[i_st_r_dim] = k_st_arr_1d + i_st_r_dim * N_p;
    }
  }

  double **r_p_i_st_arr = new double*[N_r_dim];
  double *r_p_i_st_arr_1d = new double[r_p_arr_size];
  for (int i_r_dim = 0; i_r_dim < N_r_dim; i_r_dim++) {
    r_p_i_st_arr[i_r_dim] = r_p_i_st_arr_1d + i_r_dim * N_p;
  }
  
  double *k_p; // *k_p_max;
  double a_i_k;
  int i_st, i_k;
 
  double _delta_t, _delta_inner_t; // _inner_t

  for (int i_t = 1; i_t < N_t; i_t++) {

    _delta_t = t_arr[i_t] - t_arr[i_t-1];

    //// Eval k1 (assuming v_p_arr is already k1, provided c[0] == 0)
    i_st = 0;
    std::copy(v_p_arr_1d, v_p_arr_1d + r_p_arr_size, k_st_arr[i_st][0]);


    for (i_st = 1; i_st < N_st; i_st++) {
    
      _delta_inner_t = c[i_st-1] - c[i_st];

      // Eval y_i_st
      std::copy(r_p_arr_1d, r_p_arr_1d + r_p_arr_size, r_p_i_st_arr_1d);
      for (i_k=0; i_k < i_st; i_k++) 
      {
        a_i_k = A[i_st][i_k];

        for(r_p_p = r_p_i_st_arr_1d, 
            r_p_p_max = r_p_i_st_arr_1d + r_p_arr_size, 
            k_p = k_st_arr[i_k][0]; 
            r_p_p < r_p_p_max; r_p_p++, k_p++) 
        {
          *r_p_p += a_i_k * *k_p;
        }

      } // for-loop : i_k
      

      // Eval psi at t_i_st = t_n + c_i_st * _delta_t
      if (_delta_inner_t != 0) {
        return_code = propa_psi_arr(
            psi_arr, _delta_inner_t*_delta_t, N_rho,N_lm,qprop_dim,initial_m);
        if (return_code != EXIT_SUCCESS) {
          return return_with_mesg("Failed during propagation of 'psi_arr'");
        }
      }


      // Eval k
      return_code = eval_v_p_arr_for_sph_harm_basis(
          N_s, N_p, N_rho, N_lm, r_p_i_st_arr, (const m_t **) psi_arr, 
          rho_arr, l_arr, m_arr, rho_p_lim, k_st_arr[i_st]);
      if (return_code != EXIT_SUCCESS) {
        return return_with_mesg(
            "Failed to run 'Failed to run 'eval_psi_and_dpsidx_arr()");
      }

    } // for-loop : i_st


    for (i_st = 0; i_st < N_st; i_st++) {
      for(r_p_p = r_p_arr_1d, 
          r_p_p_max = r_p_arr_1d + r_p_arr_size, 
          k_p = k_st_arr[i_st][0]; 
          r_p_p < r_p_p_max; r_p_p++, k_p++) 
      {
        *r_p_p += _delta_t * b[i_st] * *k_p;
      }
    } // for-loop : i_st

    //// Evaluate velocity vector for each particle
    return_code = eval_v_p_arr_for_sph_harm_basis(
        N_s, N_p, N_rho, N_lm, r_p_arr, (const m_t **) psi_arr, 
        rho_arr, l_arr, m_arr, rho_p_lim, v_p_arr);
    if (return_code != EXIT_SUCCESS) {
      return return_with_mesg("Failed to run 'eval_psi_and_dpsidx_arr()");
    }

    //// Store `r_p_arr` and `v_p_arr`
    std::copy(
        r_p_arr_1d, r_p_arr_1d + r_p_arr_size,
        r_p_t_arr_1d + (i_t) * r_p_arr_size);
    std::copy(
        v_p_arr_1d, v_p_arr_1d + r_p_arr_size,
        v_p_t_arr_1d + (i_t) * r_p_arr_size);

  } // for-loop : i_t

  delete [] k_st_arr_1d;
  delete [] k_st_arr_2d;
  delete [] k_st_arr;
  delete [] r_p_i_st_arr_1d;
  delete [] r_p_i_st_arr;



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

  printf("r_p_t0_arr: \n");
  print_3vec_p_arr(r_p_t0_arr_1d,N_p,N_r_dim);

  print_bar();

  // Print 'v_p' and its analytical counterpart
  printf(
      "%5s%15s%45s\n", 
      "i_p", "v_p", "v_p_ana");

  for (int i_p = 0; i_p < N_p; i_p++) {
    printf("%5d%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f\n", 
        i_p, 
        v_p_arr[0][i_p], v_p_arr[1][i_p], v_p_arr[2][i_p],
        v_p_ana_arr[0][i_p], v_p_ana_arr[1][i_p], v_p_ana_arr[2][i_p]
        );
  }

  print_bar();

  printf("r_p_t_arr: \n");
  for (int i_t = 0; i_t < N_t; i_t++) {
    print_3vec_p_arr(r_p_t_arr[i_t],N_p,N_r_dim);
    printf("\n");
  } 

  printf("v_p_t_arr: \n");
  for (int i_t = 0; i_t < N_t; i_t++) {
    print_3vec_p_arr(v_p_t_arr[i_t],N_p,N_r_dim);
    printf("\n");
  } 
  

  // Deallocate memory for arrays
  delete [] r_p_arr[0];
  delete [] r_p_arr;
  delete [] v_p_arr[0];
  delete [] v_p_arr;
  delete [] psi_arr[0];
  delete [] psi_arr;
  delete [] r_p_t_arr_1d;
  delete [] r_p_t_arr;
  delete [] r_p_t0_arr_1d;


  //// Return from this program
  return EXIT_SUCCESS;
}

