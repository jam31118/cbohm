#include "../include/ode.hh"
#include "../include/velocity.hh"
#include "../include/lapack.hh"
#include "../include/log.hh"
#include "gsl/gsl_multiroots.h"

struct g_param {
  const int N_s;
  const int N_rho;
  const int N_lm;
  const std::complex<double> **psi_t_next_arr_arr;
  const double *rho_arr;
  const int *l_arr;
  const int *m_arr;
  const double *rho_p_lim;
  double delta_t;
  double r_p_t_vec[DIM_R];
};


int implicit_euler_fzero_func(
    const gsl_vector *r_p_gsl_vec, void *g_param_p, gsl_vector *g)
{
  //// Process function arguments 
  struct g_param *p = (struct g_param *) g_param_p;
  vec_t r_p_vec;
  for (int i=0; i<DIM_R; i++) 
  { r_p_vec[i] = gsl_vector_get(r_p_gsl_vec, i); }

  //// Evaluate v_p_vec
  int return_code; vec_t v_p_vec;
  return_code = eval_v_p_for_sph_harm_basis(
      p->N_s, p->N_rho, p->N_lm, r_p_vec, p->psi_t_next_arr_arr, p->rho_arr, 
      p->l_arr, p->m_arr, p->rho_p_lim, v_p_vec, NULL);
  if (return_code != EXIT_SUCCESS) 
  { return debug_mesg("Failed to eval_v_p_for_sph_harm_basis()"); } 

  //// Evaluate target function value(s)
  double g_i;
  for (int i=0; i<DIM_R; i++) {
    g_i = r_p_vec[i] - p->r_p_t_vec[i] - p->delta_t*v_p_vec[i];
    gsl_vector_set(g, i, g_i);
  }

  //// Return if everything goes well
  return GSL_SUCCESS;

}


void print_state(size_t i_iter, gsl_multiroot_fsolver *s) {                     
  printf("[ i_iter = %03lu ] x = %7.3f %7.3f %7.3f / f(x) = %10.3e %10.3e %10.3e\n", 
      i_iter,
      gsl_vector_get(s->x,0), gsl_vector_get(s->x,1), gsl_vector_get(s->x,2),
      gsl_vector_get(s->f,0), gsl_vector_get(s->f,1), gsl_vector_get(s->f,2));    
}         



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
  

  //// Set up solver
  const gsl_multiroot_fsolver_type *T; 
  gsl_multiroot_fsolver *s;
  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc(T, DIM_R);
  struct g_param g_param = {
    N_s, N_rho, N_lm, psi_t_next_arr_arr, rho_arr, l_arr, m_arr, rho_p_lim, 
    delta_t, {r_p_t_vec[0], r_p_t_vec[1], r_p_t_vec[2]}
  };
  gsl_multiroot_function f = {&implicit_euler_fzero_func, DIM_R, &g_param};
//  const double x_init[2] = {-10.0, -5.0};
  gsl_vector *r_p_initial = gsl_vector_alloc(DIM_R);
  for (int i=0; i<DIM_R; i++) 
  { gsl_vector_set(r_p_initial, i, r_p_t_vec[i]); }
  gsl_multiroot_fsolver_set(s, &f, r_p_initial);



  //// Iterate                                                                  
  int iter_status, resi_status;                                                 
  int i_iter = 0;                                                            
  print_state(i_iter, s);                                                       
  do                                                                            
  {                                                                             
    i_iter++;                                                                   
    iter_status = gsl_multiroot_fsolver_iterate(s);                             
    print_state(i_iter, s);                                                     
    if (iter_status != 0) { break; } // check if the solver is stuck            
    resi_status = gsl_multiroot_test_residual(s->f, 1e-7);                      
  }                                                                             
  while (resi_status == GSL_CONTINUE && i_iter < 1000);                         
  printf("residual status = %s\n", gsl_strerror(resi_status));

  //// Store result
  for (int i=0; i<DIM_R; i++) 
  { r_p_t_vec[i] = gsl_vector_get(s->x,i); }
  

 //// Deallocate
  gsl_multiroot_fsolver_free(s);
  gsl_vector_free(r_p_initial);
       
  return EXIT_SUCCESS;

  
// 
//  //// Check argument
//  if (thres <= 0.0)
//  { return debug_mesg("invalid argument(s)"); }
//
//
//  //// Declare variables
//  int return_code = EXIT_FAILURE;
//  
//  vec_t v_p_vec, negative_g, r_p_vec;
//  double *delta_r_p_vec;
//  jac_t jac_v, jac_g;
//
//  const int n = DIM_R, nrhs = 1;
//  int ipiv[n], gesv_info;
//
//
//  //// Initialize `r_p_vec`
//  for (int i_dim=0; i_dim < DIM_R; i_dim++)
//  { r_p_vec[i_dim] = r_p_t_vec[i_dim]; }
//
//
////  printf("%60s\n", "r_p_vec: ");
//
//  //// Iteration starts here
//  for (i_iter = 0; i_iter < NEWTON_ITER_MAX; i_iter++) {
//
////    printf("[ LOG @ %s:%d ] right before eval_v_p_for_sph_harm_basis()\n", __FILE__,__LINE__);
////    if (jac_v == NULL) {printf("[ LOG @ %s:%d ] jac_v is NULL\n", __FILE__,__LINE__);}
////    else {printf("[ LOG @ %s:%d ] jac_v is not NULL\n", __FILE__,__LINE__);}
//    
//    //// Evaluate velocity vector and Jacobian of vector
//    eval_v_p_for_sph_harm_basis(
//        N_s, N_rho, N_lm, r_p_vec, psi_t_next_arr_arr, rho_arr, 
//        l_arr, m_arr, rho_p_lim, v_p_vec, jac_v);
//  
//    //// Evaluate negative_g array (the right hand side of linear system eq.)
//    for (int i_dim = 0; i_dim < DIM_R; i_dim++) {
//      negative_g[i_dim] = - (
//          r_p_vec[i_dim] - r_p_t_vec[i_dim] - delta_t*v_p_vec[i_dim]);
//    }
//
//if (verbose) {
//
//    printf("[ LOG ][i_iter=%05d] r_p_vec: \n", i_iter);
//    for (int i_dim = 0; i_dim < DIM_R; i_dim++) {
//      printf("%20.15f", r_p_vec[i_dim]);
//    } // end-for-loop : `i_dim`
//    printf("%20.15f\n", vec_norm(negative_g));
////    printf("\n");
//    
////    printf("[ LOG ][i_iter=%05d] jac_v: \n", i_iter);
////    print_jac_t(jac_v);
//
//}
//    
//
//    //// Terminate iteration
//    if (vec_norm(negative_g) < thres) { break; }
//  
//    //// Evaluate Jacobian for function 'g'
//    for (int i_row = 0; i_row < DIM_R; i_row++) {
//      for (int i_col = 0; i_col < DIM_R; i_col++) {
//        jac_g[i_col][i_row] = (i_row==i_col) - delta_t * jac_v[i_row][i_col];
//      }
//    }
//  
//    //// Solve linear system (Jacobian of 'g') to get delta_r_p_vec
//    dgesv_(&n, &nrhs, jac_g[0], &n, ipiv, negative_g, &n, &gesv_info);
//    return_code = handle_gesv_info(gesv_info);
//    if (return_code != EXIT_SUCCESS)
//    { return debug_mesg("Failed during checking `Xgesv_` info"); }
//  
//    //// aliasing
//    delta_r_p_vec = negative_g;
//  
//    for (int i_dim = 0; i_dim < DIM_R; i_dim++) {
//      r_p_vec[i_dim] += delta_r_p_vec[i_dim];
//    } // end-for-loop : `i_dim`
//
//  } // end-for-loop : `i_iter`
//
//
//  for (int i_dim = 0; i_dim < DIM_R; i_dim++) {
//    r_p_t_vec[i_dim] = r_p_vec[i_dim];
//  } // end-for-loop : `i_dim`
//
//
//  //// Check whether the number of iteration has been exceeded a given value
//  if (i_iter >= NEWTON_ITER_MAX) 
//  { 
//    fprintf(stderr,
//        "[i_iter=%03d] r_p_vec: %20.15f%20.15f%20.15f / norm: %20.15f\n",
//        i_iter, r_p_vec[0],r_p_vec[1],r_p_vec[2],vec_norm(negative_g));
//    fprintf(stderr,
//        "[i_iter=%03d] v_p_vec: %20.15f%20.15f%20.15f / norm: %20.15f\n",
//        i_iter, v_p_vec[0],v_p_vec[1],v_p_vec[2],vec_norm(negative_g));
//    return debug_mesg("NEWTON_ITER_MAX exceeded"); 
//  }
//
//  //// Return if everything is alright
//  return EXIT_SUCCESS;
//
}

