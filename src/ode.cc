#include "../include/ode.hh"
#include "../include/velocity.hh"
#include "../include/lapack.hh"
#include "../include/log.hh"

int prop_implicit_euler_in_sph_harm_basis(
    const int N_s, const int N_rho, const int N_lm,
    const std::complex<double> **psi_t_next_arr_arr,
    const double *rho_arr, const int *l_arr, const int *m_arr,
    const double *rho_p_lim, double delta_t, double thres, 
    double r_p_t_vec[DIM_R], bool verbose)
{


//  printf("[ LOG @ %s:%d:%s() ] r_p_t_vec %7.3f%7.3f%7.3f: \n", __FILE__,__LINE__,__func__, r_p_t_vec[0], r_p_t_vec[1], r_p_t_vec[2]);
//  printf("%7.3f%7.3f%7.3f\n", r_p_t_vec[0], r_p_t_vec[1], r_p_t_vec[2]);

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
  

  //== Description ==
  //
  // Target function g(r), of root-finding with respect to 'r'
  //
  // g(r) = r - r_n - h*v(t_{n+1},r)
  // 
  // J_g = ( I - h*J_v )
  //
  // J_v * delta_r = - g(r)
  //
  //=================
  
 
  //// Check argument
  if (thres <= 0.0)
  { return debug_mesg("invalid argument(s)"); }


  //// Declare variables
  int return_code = EXIT_FAILURE;
  
  vec_t v_p_vec, negative_g, r_p_vec, r_p_vec_trial, v_p_vec_trial, negative_g_trial;
  double *delta_r_p_vec;
  jac_t jac_v, jac_g;

  const int n = DIM_R, nrhs = 1;
  int ipiv[n], gesv_info;
  float attenuation = 1.0;
  const float attenuation_ratio = 0.8;


  //// Initialize `r_p_vec`
  for (int i_dim=0; i_dim < DIM_R; i_dim++)
  { r_p_vec[i_dim] = r_p_t_vec[i_dim]; }


//  printf("%60s\n", "r_p_vec: ");

  //// Iteration starts here
  int i_iter;
  for (i_iter = 0; i_iter < NEWTON_ITER_MAX; i_iter++) {

//    printf("[ LOG @ %s:%d ] right before eval_v_p_for_sph_harm_basis()\n", __FILE__,__LINE__);
//    if (jac_v == NULL) {printf("[ LOG @ %s:%d ] jac_v is NULL\n", __FILE__,__LINE__);}
//    else {printf("[ LOG @ %s:%d ] jac_v is not NULL\n", __FILE__,__LINE__);}
    
    //// Evaluate velocity vector and Jacobian of vector
    eval_v_p_for_sph_harm_basis(
        N_s, N_rho, N_lm, r_p_vec, psi_t_next_arr_arr, rho_arr, 
        l_arr, m_arr, rho_p_lim, v_p_vec, jac_v);
  
    //// Evaluate negative_g array (the right hand side of linear system eq.)
    for (int i_dim = 0; i_dim < DIM_R; i_dim++) {
      negative_g[i_dim] = - (
          r_p_vec[i_dim] - r_p_t_vec[i_dim] - delta_t * v_p_vec[i_dim]);
    }


    if (verbose) {
      printf("[ LOG ][i_iter=%05d] r_p_vec: \n", i_iter);
      for (int i_dim = 0; i_dim < DIM_R; i_dim++) {
        printf("%20.15f", r_p_vec[i_dim]);
      } // end-for-loop : `i_dim`
      printf("%20.15f\n", vec_norm(negative_g));
//      printf("\n");
//      printf("[ LOG ][i_iter=%05d] jac_v: \n", i_iter);
//      print_jac_t(jac_v);
    }
    

    //// Terminate iteration
    if (vec_norm(negative_g) < thres) { break; }
  
    //// Evaluate Jacobian for function 'g'
    for (int i_row = 0; i_row < DIM_R; i_row++) {
      for (int i_col = 0; i_col < DIM_R; i_col++) {
        jac_g[i_col][i_row] = (i_row==i_col) + delta_t - jac_v[i_row][i_col];
      }
    }
  
    //// Solve linear system (Jacobian of 'g') to get delta_r_p_vec
    dgesv_(&n, &nrhs, jac_g[0], &n, ipiv, negative_g, &n, &gesv_info);
    return_code = handle_gesv_info(gesv_info);
    if (return_code != EXIT_SUCCESS)
    { return debug_mesg("Failed during checking `Xgesv_` info"); }
  
    //// aliasing
    delta_r_p_vec = negative_g;


    //// Determine step size and commit newton iteration
    int i_trial;
    for (i_trial = 0; i_trial < NEWTON_LINE_SEARCH_ITER_MAX; i_trial++) {
  
      for (int i_dim = 0; i_dim < DIM_R; i_dim++) {
        r_p_vec_trial[i_dim] = r_p_vec[i_dim] + attenuation * delta_r_p_vec[i_dim];
      } // end-for-loop : `i_dim`


      //// Test vec_norm
      eval_v_p_for_sph_harm_basis(
          N_s, N_rho, N_lm, r_p_vec_trial, psi_t_next_arr_arr, rho_arr, 
          l_arr, m_arr, rho_p_lim, v_p_vec_trial, NULL);

      for (int i_dim = 0; i_dim < DIM_R; i_dim++) {
        negative_g_trial[i_dim] = - (
            r_p_vec_trial[i_dim] - r_p_t_vec[i_dim] - delta_t*v_p_vec_trial[i_dim]);
      }
//      if (vec_norm(r_p_vec_trial) < vec_norm(r_p_vec)) { break; } 
      if (vec_norm(negative_g_trial) < vec_norm(negative_g)) { break; } 
      else { attenuation *= attenuation_ratio; }

    }
    if (i_trial >= NEWTON_LINE_SEARCH_ITER_MAX) 
    { return debug_mesg("NEWTON_LINE_SEARCH_ITER_MAX exceeded"); }


    //// Update `r_p_vec` (with determined attenuation)
    for (int i_dim = 0; i_dim < DIM_R; i_dim++) {
      r_p_vec[i_dim] = r_p_vec_trial[i_dim];
    } // end-for-loop : `i_dim`


  } // end-for-loop : `i_iter`



//  for (int i_trial = 0; i_trial < NEWTON_LINE_SEARCH_ITER_MAX; i_trial++) {
//    for (int i_dim = 0; i_dim < DIM_R; i_dim++) {
//      r_p_vec_trial[i_dim] = attenuation * r_p_vec[i_dim
//    } // end-for-loop : `i_dim`
//  }


  for (int i_dim = 0; i_dim < DIM_R; i_dim++) {
    r_p_t_vec[i_dim] = r_p_vec[i_dim];
  } // end-for-loop : `i_dim`


  //// Check whether the number of iteration has been exceeded a given value
  if (i_iter >= NEWTON_ITER_MAX) 
  { 
    fprintf(stderr,
        "[i_iter=%03d] r_p_vec: %20.15f%20.15f%20.15f / norm: %20.15f\n",
        i_iter, r_p_vec[0],r_p_vec[1],r_p_vec[2],vec_norm(negative_g));
    fprintf(stderr,
        "[i_iter=%03d] v_p_vec: %20.15f%20.15f%20.15f / norm: %20.15f\n",
        i_iter, v_p_vec[0],v_p_vec[1],v_p_vec[2],vec_norm(negative_g));
    return debug_mesg("NEWTON_ITER_MAX exceeded"); 
  }

  //// Return if everything is alright
  return EXIT_SUCCESS;

}

