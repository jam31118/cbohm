#include "../include/velocity.hh"
#include "../include/log.hh"
#include "../include/lapack.hh"
#include "gsl/gsl_sf_legendre.h"


int eval_shift_offset(
    const int num_of_stencils, const int i_x_s_at_nls, 
    const int N_x, int *const offset) {

  //// Declaration
  int shift_offset_to_right, shift_offset_to_left;
  bool do_shift_to_right, do_shift_to_left;

  //// Define required values
  const int i_x_max = N_x - 1;
  const int N_s_on_left = (num_of_stencils / 2) - 1;
  const int N_s_on_right = num_of_stencils - N_s_on_left - 1;

  shift_offset_to_right = N_s_on_left - i_x_s_at_nls;
  shift_offset_to_left = (i_x_max - i_x_s_at_nls) - N_s_on_right;

  do_shift_to_right = shift_offset_to_right > 0;
  do_shift_to_left = shift_offset_to_left < 0;

  //// Check if shift occurs at both direction, which is abnormal
  if (do_shift_to_right and do_shift_to_left) 
  { return debug_mesg("Cannot shift to both direction"); }

  //// Evaluate the offset
  *offset 
    = do_shift_to_right * shift_offset_to_right 
    + do_shift_to_left * shift_offset_to_left;

  //// Return if everything is fine
  return EXIT_SUCCESS;

}


template <typename T>
int eval_psi_deriv_p(
    const double x_p, const T *psi_arr, const double *x_arr, 
    const int N_s, const int N_x, const double *x_p_lim, const int N_o,
    const int *deriv_order_arr, T *psi_deriv_p_arr) {
  
  //// Check arguments
  if ( !(N_s > 0 and N_x > N_s and N_o > 0 and N_o <= N_s) )
  { return debug_mesg("Invalid argument"); }
  for (int i_o = 0; i_o < N_o; i_o++) {
    if (deriv_order_arr[i_o] >= N_s) { 
      fprintf(stderr,"[ERROR] invalid deriv_order_arr[%d]: '%d' (>=N_S=%d)\n",
          i_o, deriv_order_arr[i_o], N_s);
      return debug_mesg("Invalid derivative order found"); 
    }
  }

  //// Define some useful variables
  const double delta_x = x_arr[1] - x_arr[0];
  
  //// Aliasing
  const int num_of_stencils = N_s;


  //// Explanation on stencil system
  //
  // Consider a case where four stencils are used for each finite
  // difference approximation, denoting each stencil with 'O'
  // in the following figure:
  //
  //  [ ]  [i_x_0][i_x_1][i_x_2][i_x_3]  [ ]    [ ]   (spatial grid indice)
  //   .      .      .      .      .      .      .    (spatial grid points)
  //         [0]    [1]    [2]    [3]                 (stencil indice)
  //          O      O      O      O                  (stencils)
  //                   X                              (particle position)
  //
  // where the numbers in parentheses `[]` on the third line 
  // refer to the index of each stencil, 
  // `X` denotes the position of the particle,
  // and several '.' represent the spatial grid points,
  // identified by index `i_x`. The `i_x_s` with `s=0,1,...,N_s-1`
  // denotes the index of spatial grid points where the stencils are placed.
  // 
  // The `index_of_nearest_left_stencil` is the index of 
  // the nearest-left stencil ('nls'), which is `1` in this case.
  // Using an alias, the following is true: `i_nls == 1`,
  // and thus, the spatial index at the nearest-left-stencil,
  // denoted by `i_x_at_nls` is `i_x_1` in this case.
  // 
  // The `N_s_on_left` is the number of stencils on the left
  // of the nearest-left-stencil, and there is one stencil, `[0]`
  // in this case.
  // 
  // The `N_s_on_right` is the number of stencils on the right
  // of the nearest-left-stencil, and there are two: `[2]` and `[3]`
  // in this case.


  const int index_of_nearest_left_stencil = (num_of_stencils / 2) - 1;
  const int i_nlp = index_of_nearest_left_stencil; // aliasing

  int i_x_s_arr[num_of_stencils];
  int i_s, *i_x_s_arr_p, *i_x_s_arr_p_max = i_x_s_arr + N_s;
  int i_x_s_at_nls;

  int shift_offset;

  double power_matrix[num_of_stencils][num_of_stencils];
  double *power_matrix_1d = power_matrix[0];

  const int b_vec_matrix_col_num = N_o; // for psi and dpsidx respectively.
  double *b_vec_matrix[b_vec_matrix_col_num];
  double b_vec_matrix_1d[b_vec_matrix_col_num * num_of_stencils];
  for (int i_col = 0; i_col < b_vec_matrix_col_num; i_col++)
  { b_vec_matrix[i_col] = b_vec_matrix_1d + i_col * num_of_stencils; }
  std::memset(b_vec_matrix_1d, 0, sizeof(b_vec_matrix_1d));

  // An alias for `b_vec_matrix`
  // Note that the `b_vec_matrix` holds the coefficients after gesv() routine
  // This `coef_vec_matrix` will be used after the gesv() routine.
  double **coef_vec_matrix = b_vec_matrix;
  
  double delta_x_s;

  // The pivot indices that define the permutation matrix P;
  // row i of the matrix was interchanged with row IPIV(i).
  int pivot_indices[num_of_stencils];
  int gesv_info;

  T psi_x_s;
  
  int deriv_order;
  
  //// Check in-range
  if (x_p_lim[0] > x_p or x_p >= x_p_lim[1]) { return CODE_OUT_OF_RANGE; }

  //// Evaluate the index of grid point of nearest left stencil
  i_x_s_at_nls = (x_p - x_arr[0]) / delta_x;

  //// Evaluate shift offset for stencils
  if (eval_shift_offset(N_s,i_x_s_at_nls,N_x,&shift_offset) != EXIT_SUCCESS)
  { return debug_mesg("Failed to evaluate 'shift_offset'"); }


  for (
      i_s = 0, i_x_s_arr_p = i_x_s_arr;
      i_x_s_arr_p < i_x_s_arr_p_max;
      ++i_x_s_arr_p, ++i_s
      ) 
  {
    *i_x_s_arr_p = i_x_s_at_nls + (i_s - i_nlp) + shift_offset; 
    
    // Construct the power matrix - in column major for FORTRAN LAPACK
    power_matrix[i_s][0] = 1.0; 
    delta_x_s = x_arr[*i_x_s_arr_p] - x_p;
    for (int i_row = 1; i_row < N_s; ++i_row) {
      power_matrix[i_s][i_row] = power_matrix[i_s][i_row-1] * delta_x_s;
    } // end for loop : `i_row`

  } // for-loop `i_s`

  for (int i_order = 0; i_order < b_vec_matrix_col_num; i_order++) {
    deriv_order = deriv_order_arr[i_order];
    b_vec_matrix[i_order][deriv_order] = 1.0; 
  } // end for loop : `i_order`

 
#ifdef DEBUG
  printf("x_p = %7.3f\n",x_p);
  for (int i_row = 0; i_row < N_s; i_row++) {
    printf("%7.3f",x_arr[i_x_s_arr[i_row]]);
    for (int i_col = 0; i_col < N_s; i_col++)
    { printf("%7.3f", power_matrix[i_col][i_row]); }
    for (int i_o = 0; i_o < N_o; i_o++) 
    { printf("%7.3f", b_vec_matrix[i_o][i_row]); }
    printf("\n");
  }

  printf("power_matrix_1d: ");
  for (int i=0; i < N_s*N_s; i++) {
    printf("%7.3f", power_matrix_1d[i]);
  } printf("\n");
  printf("b_vec_matrix_1d: ");
  for (int i=0; i < N_s*b_vec_matrix_col_num; i++) {
    printf("%7.3f", b_vec_matrix_1d[i]);
  } printf("\n");

#endif // DEBUG

  //// Solve linear system for obtaining finite-difference coefficients
  dgesv_(
      &N_s, &b_vec_matrix_col_num, power_matrix_1d, &N_s, 
      pivot_indices, b_vec_matrix_1d, &N_s, &gesv_info
  );
  if ( handle_gesv_info(gesv_info) != EXIT_SUCCESS)
  { return return_with_mesg("Failed solving for coeffcients"); }
  
  //// Evaluate the finite-difference-approximated values: psi, dpsidx
  for (int i_o = 0; i_o < N_o; i_o++) {
    psi_deriv_p_arr[i_o] = 0.0;
    for (int i_s = 0; i_s < N_s; i_s++) {
      psi_x_s = psi_arr[i_x_s_arr[i_s]];
      psi_deriv_p_arr[i_o] += coef_vec_matrix[i_o][i_s] * psi_x_s;
    } // end for loop : i_o
  } // end for loop : i_s

  return EXIT_SUCCESS;
}



template <typename T>
int eval_psi_and_dpsidx_p(
    const double x_p, const T *psi_arr, const double *x_arr, 
    const int N_s, const int N_x, const double *x_p_lim, 
    T *psi_p_arr_p, T *dpsidx_p_arr_p) {
  
  const int N_o = 2;
  const int deriv_order_arr[N_o] = { 0, 1 };
  T psi_deriv_p_arr[N_o];

  int return_code = EXIT_FAILURE;
  return_code = eval_psi_deriv_p<T>(
      x_p, psi_arr, x_arr, N_s, N_x, x_p_lim, N_o,
      deriv_order_arr, psi_deriv_p_arr);
  if (return_code != EXIT_SUCCESS)
  { return debug_mesg("Failed to `eval_psi_deriv_p()`"); }
  
  *psi_p_arr_p = psi_deriv_p_arr[0];
  *dpsidx_p_arr_p = psi_deriv_p_arr[1];

  return EXIT_SUCCESS;
}



template <typename T>
int eval_psi_and_dpsidx_arr(
    const double *x_p_arr, T *psi_arr, const double *x_arr, 
    const int N_s, const int N_p, const int N_x, const double *x_p_lim, 
    T *psi_p_arr, T *dpsidx_p_arr) {
  
  //// Check arguments
  if ( !(N_s > 0 and N_x > N_s and N_p > 0) )
  { return debug_mesg("Invalid argument"); }

  //// Function Arguments
  //
  // `N_s`: the number of stencils used for finite difference expression
  // `N_p`: the number of particles
  // `N_x`: the number of spatial grid points 
  //        on which a state function `psi_arr` is defined
  // `x_p_arr`: an array of positions of particles of length `N_p`
  // `psi_arr`: an array of state function of length `N_x`
  // `x_arr`: an array of spatial grid points of length `N_x`


  //// Check function arguments
  assert(N_s > 0 and N_p > 0 and N_x > N_s);

  const double *x_p_arr_p, *x_p_arr_p_max = x_p_arr + N_p;
  T *psi_p_arr_p, *dpsidx_p_arr_p;
  
  int return_code = EXIT_FAILURE;


  for (
      x_p_arr_p = x_p_arr, 
      psi_p_arr_p = psi_p_arr, dpsidx_p_arr_p = dpsidx_p_arr;
      x_p_arr_p < x_p_arr_p_max; 
      ++x_p_arr_p, ++psi_p_arr_p, ++dpsidx_p_arr_p
      )
  {

    return_code = eval_psi_and_dpsidx_p(
          *x_p_arr_p, psi_arr, x_arr, 
          N_s, N_x, x_p_lim, 
          psi_p_arr_p, dpsidx_p_arr_p
        );
    if (return_code != EXIT_SUCCESS) 
    { return debug_mesg("Failed during 'eval_psi_and_dpsidx_p'"); }

  } // for-loop `x_p_arr_p`

  // Returns if everything is fine.
  return EXIT_SUCCESS;

}



//// Instantiation

// Instantiation of `eval_psi_and_dpsidx_p()` 
template int eval_psi_deriv_p<double>(
    const double x_p, const double *psi_arr, const double *x_arr, 
    const int N_s, const int N_x, const double *x_p_lim, const int N_o,
    const int *deriv_order_arr, double *psi_deriv_p_arr);

template int eval_psi_deriv_p<std::complex<double>>(
    const double x_p, const std::complex<double> *psi_arr, const double *x_arr, 
    const int N_s, const int N_x, const double *x_p_lim, const int N_o,
    const int *deriv_order_arr, std::complex<double> *psi_deriv_p_arr);

// Instantiation of `eval_psi_and_dpsidx_p()` 
template int eval_psi_and_dpsidx_p<double>(
    const double x_p, const double *psi_arr, const double *x_arr, 
    const int N_s, const int N_x, const double *x_p_lim, 
    double *psi_p, double *dpsidx_p);

template int eval_psi_and_dpsidx_p< std::complex<double> >(
    const double x_p, const std::complex<double> *psi_arr, const double *x_arr, 
    const int N_s, const int N_x, const double *x_p_lim, 
    std::complex<double> *psi_p, std::complex<double> *dpsidx_p);

// Instantiation of `eval_psi_and_dpsidx_arr()`
template int eval_psi_and_dpsidx_arr<double>(
    const double *x_p_arr, double *psi_arr, const double *x_arr, 
    const int N_s, const int N_p, const int N_x, const double *x_p_lim, 
    double *psi_p_arr, double *dpsidx_p_arr);

template int eval_psi_and_dpsidx_arr< std::complex<double> >(
    const double *x_p_arr, std::complex<double> *psi_arr, const double *x_arr, 
    const int N_s, const int N_p, const int N_x, const double *x_p_lim, 
    std::complex<double> *psi_p_arr, std::complex<double> *dpsidx_p_arr);






//// Implementation of `eval_v_p_for_sph_harm_basis`
int eval_v_p_for_sph_harm_basis(
    const int N_s, const int N_rho, const int N_lm, 
    const double r_p_vec[DIM_R], 
    const std::complex<double> **psi_in_sph_harm_basis_arr,
    const double *rho_arr, int *l_arr, int *m_arr, const double *rho_p_lim, 
    double v_p_vec[DIM_R], double jac[DIM_R][DIM_R])
{
  
  //// Function Argument
  //
  // [NOTE]
  // `DIM_R` : the dimension of paritlce's position vector space
  //
  // [INPUT]
  // `N_s` : the number of stencils used in finite-difference approximation
  //         for estimating psi and dpsidrho at particle's radial coodinate
  // `N_rho` : the number of radial grid points
  // `N_lm` : the number of spherical harmonics basis
  //
  // `r_p_vec` : 1D array of shape (`DIM_R`,)
  // `psi_in_sph_harm_basis_arr` : 2D array of shape (`N_lm`,`N_rho`)
  // `rho_arr` : 1D array of shape (`N_rho`,)
  // `l_arr` : 1D array of shape (`l_arr`,)
  // `m_arr` : 1D array of shape (`m_arr`,)
  // `rho_p_lim` : 1D array of shape (2,)
  //
  // [OUTPUT]
  // `v_p_vec` : 1D array of shape (`DIM_R`,)
  //

  //// Check arguments
  if ( !(N_s > 0 and N_rho > N_s and N_lm > 0) )
  { return debug_mesg("Invalid argument"); }

  
  //// Define useful variables
  const double rho_p_min = rho_p_lim[0], rho_p_max = rho_p_lim[1];
  int return_code = EXIT_FAILURE;

  //// Declare variables
  std::complex<double> 
    psi_p_Y = 0, dpsi_p_Y = 0, psi_p_dtheta_Y = 0, psi_p_dphi_Y = 0,
    d2psi_p_Y = 0, dpsi_p_dtheta_Y = 0, dpsi_p_dphi_Y = 0, psi_p_d2theta_Y = 0,
    psi_p_dphi_dtheta_Y = 0, psi_p_d2phi_Y = 0;
//    psi_p_dtheta_dphi_Y = 0, 

  double rho_p = r_p_vec[0], theta_p = r_p_vec[1], phi_p = r_p_vec[2];
  bool out_of_range;

  //// Check singularity
  if (rho_p == 0 or sin(theta_p) == 0) 
  { return debug_mesg("Zero 'rho_p' or 'sin(theta_p)' occurred"); }

  //// If out-of-range, the velocity is returned as zero
  out_of_range = rho_p_min > rho_p or rho_p >= rho_p_max; 
  if (out_of_range) { 
    for (int i_p_dim = 0; i_p_dim < DIM_R; i_p_dim++) 
    { v_p_vec[i_p_dim] = 0.0; }
    return EXIT_SUCCESS;
  }
  
  //// Define useful variables
  const double cos_theta_p = cos(theta_p);
  const double sin_theta_p = sin(theta_p);
  const double cot_theta_p = cos_theta_p / sin_theta_p;

  //// Aliasing
  const std::complex<double> **psi_arr_arr = psi_in_sph_harm_basis_arr;

  //// Determine the number of derivative orders to be evaluated
  const bool eval_jac = (jac != NULL);
  const int N_order_max = 3;
  int N_order;
  if (eval_jac) { N_order = N_order_max; }
  else { N_order = 2; }
  int deriv_order_arr[N_order];
  deriv_order_arr[0] = 0, deriv_order_arr[1] = 1;
  if (eval_jac) { deriv_order_arr[2] = 2; }

  //// Allocate memory
  std::complex<double> *Ylm_arr, *Yl1m_arr, *dtheta_Ylm_arr;
  std::complex<double> **psi_deriv_p_lm_arr;
  std::complex<double> *psi_deriv_p_lm_arr_1d;
  std::complex<double> **psi_deriv_p_lm_arr_p, **psi_deriv_p_lm_arr_p_max;
  std::complex<double> *psi_deriv_p_lm_arr_1d_p;
  psi_deriv_p_lm_arr = new std::complex<double>*[N_lm];
  psi_deriv_p_lm_arr_1d = new std::complex<double>[N_lm*N_order_max];

  for(psi_deriv_p_lm_arr_p = psi_deriv_p_lm_arr, 
      psi_deriv_p_lm_arr_p_max = psi_deriv_p_lm_arr + N_lm,
      psi_deriv_p_lm_arr_1d_p = psi_deriv_p_lm_arr_1d;
      psi_deriv_p_lm_arr_p < psi_deriv_p_lm_arr_p_max; 
      psi_deriv_p_lm_arr_p++, psi_deriv_p_lm_arr_1d_p += N_order_max) 
  { *psi_deriv_p_lm_arr_p = psi_deriv_p_lm_arr_1d_p; }

  Ylm_arr = new std::complex<double>[N_lm];
  Yl1m_arr = new std::complex<double>[N_lm];
  dtheta_Ylm_arr = new std::complex<double>[N_lm];

  //// Variables to be used in loop
  int i_lm;
  const std::complex<double> **psi_arr_p;
  int *l_p, *m_p;
  std::complex<double> 
    exp_phi, Ylm, Yl1m, dtheta_Ylm, 
    psi_p, dpsi_p, d2psi_p, *p_Ylm, *p_Yl1m, *p_dtheta_Ylm;
  double l, m, m_power_of_minus_1;
  int li, mi, m_sign, m_abs;

  
  //// Evaluate Ylm etc. and store them into arrays
  for (
      i_lm = 0, psi_arr_p = psi_arr_arr, l_p=l_arr, m_p=m_arr, 
      p_Ylm=Ylm_arr, p_Yl1m = Yl1m_arr, p_dtheta_Ylm = dtheta_Ylm_arr;
      i_lm < N_lm;
      i_lm++, psi_arr_p++, l_p++, m_p++, p_Ylm++, p_Yl1m++, p_dtheta_Ylm++
      ) 
  {

    li = *l_p, mi = *m_p;
    l = (double) li; m = (double) mi;
    m_sign = 1 - 2 * (1 - (mi>=0));
    m_abs = m_sign * mi;
    m_power_of_minus_1 = 1.0 - 2.0 * (mi%2);

    //// Estimate `psi_p` and `dpsidrho_p` (and more) using finite difference
    return_code = eval_psi_deriv_p< std::complex<double> >(
        rho_p, *psi_arr_p, rho_arr, N_s, N_rho, rho_p_lim, N_order,
        deriv_order_arr, psi_deriv_p_lm_arr[i_lm]);
    if (return_code != EXIT_SUCCESS)
    { return debug_mesg("Failed to 'eval_psi_and_dpsidx_p()'"); }

    //// Evaluate ingredients
    exp_phi = std::complex<double>(cos(m_abs*phi_p), sin(m_abs*phi_p));
    Ylm = gsl_sf_legendre_sphPlm(li, m_abs, cos(theta_p)) * exp_phi;
    Yl1m = gsl_sf_legendre_sphPlm(li+1, m_abs, cos(theta_p)) * exp_phi;
    if (mi < 0) {
      Ylm = m_power_of_minus_1 * std::conj(Ylm);
      Yl1m = m_power_of_minus_1 * std::conj(Yl1m);
    }
    
    //// Store results to arrays
    *p_Ylm = Ylm;
    *p_Yl1m = Yl1m;
    *p_dtheta_Ylm = (
        - cos_theta_p * (l+1.0) * Ylm
        + sqrt( (2.0*l+1.0)/(2.0*l+3.0) * (l+1.0+m) * (l+1.0-m) ) * Yl1m
        ) / sin_theta_p;

  } // end for loop : `i_lm`


  //// Evaluate summations for velocity 
  for (
      i_lm = 0, l_p=l_arr, m_p=m_arr, 
      p_Ylm = Ylm_arr, p_Yl1m = Yl1m_arr, p_dtheta_Ylm = dtheta_Ylm_arr,
      psi_deriv_p_lm_arr_1d_p = psi_deriv_p_lm_arr_1d;
      i_lm < N_lm;
      i_lm++, l_p++, m_p++, psi_deriv_p_lm_arr_1d_p += N_order_max,
      p_Ylm++, p_Yl1m++, p_dtheta_Ylm++
      )
  {
    l = (double) *l_p; m = (double) *m_p;
    Ylm = *p_Ylm; Yl1m = *p_Yl1m; dtheta_Ylm = *p_dtheta_Ylm;
    psi_p = *(psi_deriv_p_lm_arr_1d_p);
    dpsi_p = *(psi_deriv_p_lm_arr_1d_p+1);

    //// Add to target summations
    psi_p_Y += psi_p * Ylm;
    dpsi_p_Y += dpsi_p * Ylm;

    psi_p_dtheta_Y += psi_p * dtheta_Ylm;
//      += - cos_theta_p * (l+1.0) * psi_p * Ylm
//      + sqrt( (2.0*l+1.0)/(2.0*l+3.0) * (l+1.0+m) * (l+1.0-m) )
//        * psi_p * Yl1m;
//    psi_p_dtheta_Y /= sin_theta_p;

    psi_p_dphi_Y += m * psi_p * Ylm;
    
  } // i_lm
  


  std::complex<double> d2theta_Y;
  //// Evaluate summations for evaluating Jacobian
  if (eval_jac) {
    
    for (
        i_lm = 0, l_p=l_arr, m_p=m_arr, 
        p_Ylm = Ylm_arr, p_Yl1m = Yl1m_arr, p_dtheta_Ylm = dtheta_Ylm_arr,
        psi_deriv_p_lm_arr_1d_p = psi_deriv_p_lm_arr_1d;
        i_lm < N_lm;
        i_lm++, l_p++, m_p++, psi_deriv_p_lm_arr_1d_p += N_order_max,
        p_Ylm++, p_Yl1m++, p_dtheta_Ylm++
        )
    {
      l = (double) *l_p; m = (double) *m_p;
      Ylm = *p_Ylm; Yl1m = *p_Yl1m; dtheta_Ylm = *p_dtheta_Ylm;
      d2theta_Y = 
        - cot_theta_p * dtheta_Ylm 
        - ( l*(l+1) - m*m/(sin_theta_p*sin_theta_p) ) * Ylm;
      psi_p = *(psi_deriv_p_lm_arr_1d_p);
      dpsi_p = *(psi_deriv_p_lm_arr_1d_p+1);
      d2psi_p = *(psi_deriv_p_lm_arr_1d_p+2);
  
      d2psi_p_Y += d2psi_p * Ylm;
      dpsi_p_dtheta_Y += dpsi_p * dtheta_Ylm;
      dpsi_p_dphi_Y += dpsi_p * m * Ylm;
      psi_p_d2theta_Y += psi_p * d2theta_Y;
      psi_p_dphi_dtheta_Y += psi_p * m * dtheta_Ylm;
//      psi_p_dtheta_dphi_Y += psi_p * m * dtheta_Ylm;
      psi_p_d2phi_Y += psi_p * m * m * Ylm;
  
    } // i_lm

  } // end if : `eval_jac`



  //// Deallocate arrays after use
  delete [] psi_deriv_p_lm_arr_1d;
  delete [] psi_deriv_p_lm_arr;
  delete [] Ylm_arr;
  delete [] Yl1m_arr;
  delete [] dtheta_Ylm_arr;


#ifdef DEBUG

  //// Print real part
  printf("real part: \n");
  printf(
      "%15s%15s%15s%15s\n", 
      "denum","numer_p_rho","numer_p_theta","numer_p_phi");

  printf("%15.7f%15.7f%15.7f%15.7f\n",
      std::real(psi_p_Y), std::real(dpsi_p_Y),
      std::real(psi_p_dtheta_Y), std::real(psi_p_dphi_Y));
  printf("\n");


  //// Print imaginary part
  printf("imag part: \n");
  printf(
      "%15s%15s%15s%15s\n", 
      "denum","numer_p_rho","numer_p_theta","numer_p_phi");

  printf("%15.7f%15.7f%15.7f%15.7f\n",
      std::imag(psi_p_Y), std::imag(dpsi_p_Y),
      std::imag(psi_p_dtheta_Y), std::imag(psi_p_dphi_Y));
  printf("\n");

#endif // DEBUG



#ifdef DEBUG
  printf("v_p_arr: \n");
#endif // DEBUG
      

  // Check singularity
  if (psi_p_Y == 0.0) { 

    fprintf(stderr,"[ LOG ] %-25s\n","psi_p_Y_arr: ");
    fprintf(stderr,"(%20s,%20s)\n","real","imag");
    fprintf(stderr,"(%20.15f,%20.15f)\n",psi_p_Y.real(),psi_p_Y.imag());

    fprintf(stderr,"\n");

    fprintf(stderr,
        "[ LOG ] r_p_vec\n"
        "        = (%19s,%19s,%19s)\n"
        "        = (%19.15f,%19.15f,%19.15f)\n",
        "rho_p","theta_p","phi_p",
        r_p_vec[0],r_p_vec[1],r_p_vec[2]);

    fprintf(stderr,"\n");

    fprintf(stderr,
        "[ LOG ] v_p_vec\n"
        "        = (%19s,%19s,%19s)\n"
        "        = (%19.15f,%19.15f,%19.15f)\n",
        "v_rho_p","v_theta_p","v_phi_p",
        v_p_vec[0],v_p_vec[1],v_p_vec[2]);

    fprintf(stderr,"\n");
    
    const double delta_rho = rho_arr[1] - rho_arr[0];
    const int i_rho_nls = (int) (r_p_vec[0] - rho_arr[0]) / delta_rho;

    fprintf(stderr,
        "[ LOG ] radial grid point index of nearest left stencil: "
        "%d (= %.15f a.u.)\n",
        i_rho_nls, delta_rho * i_rho_nls);
    fprintf(stderr,"\n");

    fprintf(stderr,"[ LOG ] psi_arr_arr[i_lm][i_rho_nls]: \n");
    fprintf(stderr,"\n");
    fprintf(stderr,"%5s (%19s,%19s)\n","i_lm","real","imag");
    std::complex<double> _psi;
    for (int i_lm = 0; i_lm < N_lm; i_lm++) {
      _psi = psi_arr_arr[i_lm][i_rho_nls];
      fprintf(stderr,"%5d (%19.15f,%19.15f)\n",
          i_lm,std::real(_psi),std::imag(_psi));
    } fprintf(stderr,"\n");

    return debug_mesg("Zero 'psi_p' occurred"); 
  }  // endif denum == 0


  //// Evaluate velocity if not singular
  v_p_vec[0]
    = (dpsi_p_Y / psi_p_Y).imag();

  v_p_vec[1]
    = (psi_p_dtheta_Y / psi_p_Y).imag()
    / (rho_p);

  v_p_vec[2] 
    = (psi_p_dphi_Y / psi_p_Y).real()
    / (rho_p * sin_theta_p);


#ifdef DEBUG
  printf("%15.5f%15.5f%15.5f\n", 
      v_p_vec[0], v_p_vec[1], v_p_vec[2]);
#endif // DEBUG


  //// Evaluate Jacobian
  const std::complex<double> psi_p_Y_sq = psi_p_Y * psi_p_Y;
  if (eval_jac) {
    jac[0][0] = ( d2psi_p_Y / psi_p_Y - (dpsi_p_Y*dpsi_p_Y)/(psi_p_Y_sq) ).imag();
    jac[0][1] = ( dpsi_p_dtheta_Y / psi_p_Y - dpsi_p_Y*psi_p_dtheta_Y/(psi_p_Y_sq) ).imag();
    jac[0][2] = ( dpsi_p_dphi_Y / psi_p_Y - dpsi_p_Y * psi_p_dphi_Y / (psi_p_Y_sq) ).imag();
//    jac[1][0] = (-v_p_vec[0] + dpsi_p_dtheta_Y/psi_p_Y - (psi_p_dtheta_Y*dpsi_p_Y)/psi_p_Y_sq).imag() / rho_p;
    jac[1][0] = (-v_p_vec[1] + jac[0][1]) / rho_p;
    jac[1][1] = ( psi_p_d2theta_Y / psi_p_Y - psi_p_dtheta_Y*psi_p_dtheta_Y/psi_p_Y_sq ).imag() / rho_p;
    jac[1][2] = ( psi_p_dphi_dtheta_Y / psi_p_Y - psi_p_dtheta_Y*psi_p_dphi_Y / psi_p_Y_sq ).imag() / rho_p;
    jac[2][0] = - 1.0 / rho_p * (v_p_vec[2] + 1.0 / sin_theta_p * jac[0][2]);
    jac[2][1] = - cot_theta_p * v_p_vec[2] + 1.0 / sin_theta_p * jac[1][2];
    jac[2][2] = ( psi_p_d2phi_Y / psi_p_Y - psi_p_dphi_Y*psi_p_dphi_Y / psi_p_Y_sq ).imag() / (rho_p * sin_theta_p);
  }

  //// Return from thie program
  return EXIT_SUCCESS;

}




//// `eval_v_p_vec_arr_for_sph_harm_basis`
int eval_v_p_vec_arr_for_sph_harm_basis(
    const int N_s, const int N_p, const int N_rho, const int N_lm,
    double **r_p_vec_arr, const std::complex<double> **psi_in_sph_harm_basis_arr,
    double *rho_arr, int *l_arr, int *m_arr, const double *rho_p_lim, 
    double **v_p_vec_arr, jac_t *jac_p_arr)
{

  //// Function Argument
  //
  // [NOTE]
  // `DIM_R` : the dimension of paritlce's position vector space
  //
  // [INPUT]
  // `N_s` : the number of stencils used in finite-difference approximation
  //         for estimating psi and dpsidrho at particle's radial coodinate
  // `N_p` : the number of particles
  // `N_rho` : the number of radial grid points
  // `N_lm` : the number of spherical harmonics basis
  //
  // `r_p_vec_arr` : 2D array of shape (`N_p`,`DIM_R`)
  // `psi_in_sph_harm_basis_arr` : 2D array of shape (`N_lm`,`N_rho`)
  // `rho_arr` : 1D array of shape (`N_rho`,)
  // `l_arr` : 1D array of shape (`l_arr`,)
  // `m_arr` : 1D array of shape (`m_arr`,)
  // `rho_p_lim` : 1D array of shape (2,)
  //
  // [OUTPUT]
  // `v_p_vec_arr` : 2D array of shape (`N_p`,`DIM_R`)
  //

  int return_code = EXIT_FAILURE;

//  if (jac_p_arr == NULL) {
//    jac_p_arr = new jac_t[N_p];
//    for (int i_p = 0; i_p < N_p; i_p++) {
//      jac_p_arr[i_p] = NULL;
//    }
//  }
//
//  jac_t jac_p;
//  double r_p_vec[DIM_R], v_p_vec[DIM_R];

  for (int i_p = 0; i_p < N_p; i_p++) {

//    for (int i_dim = 0; i_dim < DIM_R; i_dim++)
//    { r_p_vec[i_dim] = r_p_arr[i_dim][i_p]; }
    if (jac_p_arr == NULL) { // { jac_p = jac_p_arr[i_p]; }
      return_code = eval_v_p_for_sph_harm_basis(
          N_s, N_rho, N_lm, r_p_vec_arr[i_p], psi_in_sph_harm_basis_arr,
          rho_arr, l_arr, m_arr, rho_p_lim, v_p_vec_arr[i_p]);

    } else {
      return_code = eval_v_p_for_sph_harm_basis(
          N_s, N_rho, N_lm, r_p_vec_arr[i_p], psi_in_sph_harm_basis_arr,
          rho_arr, l_arr, m_arr, rho_p_lim, v_p_vec_arr[i_p], jac_p_arr[i_p]);
    
    }

    if (return_code != EXIT_SUCCESS) 
    { return debug_mesg("Failed during 'eval_v_p_for_sph_harm_basis()'"); }

//    for (int i_dim = 0; i_dim < DIM_R; i_dim++)
//    { v_p_arr[i_dim][i_p] = v_p_vec[i_dim]; }

  } // for-loop : i_p

  return EXIT_SUCCESS;

}





//// `eval_v_p_arr_for_sph_harm_basis`
int eval_v_p_arr_for_sph_harm_basis(
    const int N_s, const int N_p, const int N_rho, const int N_lm,
    double **r_p_arr, const std::complex<double> **psi_in_sph_harm_basis_arr,
    double *rho_arr, int *l_arr, int *m_arr, const double *rho_p_lim, 
    double **v_p_arr)
{

  //// Function Argument
  //
  // [NOTE]
  // `DIM_R` : the dimension of paritlce's position vector space
  //
  // [INPUT]
  // `N_s` : the number of stencils used in finite-difference approximation
  //         for estimating psi and dpsidrho at particle's radial coodinate
  // `N_p` : the number of particles
  // `N_rho` : the number of radial grid points
  // `N_lm` : the number of spherical harmonics basis
  //
  // `r_p_arr` : 2D array of shape (`DIM_R`,`N_p`)
  // `psi_in_sph_harm_basis_arr` : 2D array of shape (`N_lm`,`N_rho`)
  // `rho_arr` : 1D array of shape (`N_rho`,)
  // `l_arr` : 1D array of shape (`l_arr`,)
  // `m_arr` : 1D array of shape (`m_arr`,)
  // `rho_p_lim` : 1D array of shape (2,)
  //
  // [OUTPUT]
  // `v_p_arr` : 2D array of shape (`DIM_R`,`N_p`)
  //

  int return_code = EXIT_FAILURE;

  double r_p_vec[DIM_R], v_p_vec[DIM_R];

  for (int i_p = 0; i_p < N_p; i_p++) {

    for (int i_dim = 0; i_dim < DIM_R; i_dim++)
    { r_p_vec[i_dim] = r_p_arr[i_dim][i_p]; }

    return_code = eval_v_p_for_sph_harm_basis(
        N_s, N_rho, N_lm, r_p_vec, psi_in_sph_harm_basis_arr,
        rho_arr, l_arr, m_arr, rho_p_lim, v_p_vec);
    if (return_code != EXIT_SUCCESS) 
    { return debug_mesg("Failed during 'eval_v_p_for_sph_harm_basis()'"); }

    for (int i_dim = 0; i_dim < DIM_R; i_dim++)
    { v_p_arr[i_dim][i_p] = v_p_vec[i_dim]; }

  } // for-loop : i_p

  return EXIT_SUCCESS;

}

