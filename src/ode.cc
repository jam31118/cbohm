#include "../include/ode.hh"
#include "../include/velocity.hh"
#include "../include/lapack.hh"
#include "../include/log.hh"

int prop_implicit_euler_in_sph_harm_basis(
    const int N_s, const int N_rho, const int N_lm,
    const std::complex<double> **psi_t_next_arr_arr,
    const double *rho_arr, const int *l_arr, const int *m_arr,
    const double *rho_p_lim, double delta_t, double thres, 
    double r_p_t_vec[DIM_R]) 
{

  //// Function Argument(s)
  // 
  // [INPUT]
  // `N_s`
  // `N_rho`: the number of radial grid points
  // `N_lm`: the number of spherical harmonics
  // `psi_t_next_arr_arr`: 2D array of shape (`N_lm`,`N_rho`)
  // `rho_arr`
  // `l_arr`
  // `m_arr`
  // `rho_p_lim`
  // `delta_t`: stepsize of propagation
  // `thres`
  //
  // [INPUT/OUTPUT]
  // `r_p_t_vec`: position vector of a particle of shape (`DIM_R`,)
  //
  // [OUTPUT]
  

  //== Description ==
  //
  // Target function g(r), of root-finding with respect to 'r'
  //
  // g(r) = r - r_n - h*v(t_{n+1},r)
  //
  // ( I - h*J_v ) delta_r = - g(r)
  //
  //=================
  

  if (thres <= 0.0)
  { return debug_mesg("invalid argument(s)"); }

  int return_code = EXIT_FAILURE;
  
  vec_t v_p_vec;
  jac_t jac_v;
  jac_t jac_g;
  vec_t negative_g;
  double *delta_r_p_vec;

  double r_p_vec[DIM_R];
  const int n = DIM_R;
  const int nrhs = 1;
  int ipiv[n];
  int gesv_info;

  //// Initialize `r_p_vec`
  for (int i_dim=0; i_dim < DIM_R; i_dim++)
  { r_p_vec[i_dim] = r_p_t_vec[i_dim]; }


  int i_iter;
  for (i_iter = 0; i_iter < NEWTON_ITER_MAX; i_iter++) {

    eval_v_p_for_sph_harm_basis(
        N_s, N_rho, N_lm, r_p_vec, psi_t_next_arr_arr, rho_arr, 
        l_arr, m_arr, rho_p_lim, v_p_vec, jac_v
        );
  
    for (int i_dim = 0; i_dim < DIM_R; i_dim++) {
      negative_g[i_dim] = - (
          r_p_vec[i_dim] - r_p_t_vec[i_dim] - delta_t*v_p_vec[i_dim]);
    }
    
//    printf("[i_iter=%d] norm of `negative_g`: %.20f / r_p_vec=[%.20f,%.20f,%.20f]\n", 
//        i_iter, vec_norm(negative_g), r_p_vec[0], r_p_vec[1], r_p_vec[2]);
    if (vec_norm(negative_g) < thres) { break; }
  
    for (int i_row = 0; i_row < DIM_R; i_row++) {
      for (int i_col = 0; i_col < DIM_R; i_col++) {
        jac_g[i_col][i_row] = (i_row==i_col) + delta_t * jac_v[i_row][i_col];
      }
    }
  
    dgesv_(&n, &nrhs, jac_g[0], &n, ipiv, negative_g, &n, &gesv_info);
    return_code = handle_gesv_info(gesv_info);
    if (return_code != EXIT_SUCCESS)
    { return debug_mesg("Failed during checking `Xgesv_` info"); }
  
    delta_r_p_vec = negative_g;
  
    for (int i_dim = 0; i_dim < DIM_R; i_dim++) {
      r_p_vec[i_dim] += delta_r_p_vec[i_dim];
    }

  }

  for (int i_dim = 0; i_dim < DIM_R; i_dim++) {
    r_p_t_vec[i_dim] = r_p_vec[i_dim];
  }

  //// Check whether the number of iteration has been exceeded a given value
  if (i_iter >= NEWTON_ITER_MAX) 
  { return debug_mesg("NEWTON_ITER_MAX exceeded"); }

  return EXIT_SUCCESS;
}

