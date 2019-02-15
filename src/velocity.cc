#include "../include/velocity.hh"
#include "../include/log.hh"
#include "../include/lapack.hh"
#include "gsl/gsl_sf_legendre.h"
#include "gsl/gsl_errno.h"


double vec_norm(vec_t vec) {
  double _norm = 0.0;
  for (int i=0; i<DIM_R; i++) {
    _norm += vec[i] * vec[i];
  }
  return sqrt(_norm);
}


int factorial(int n, int *n_fac) {
  if (n < 0) { return EXIT_FAILURE; }
  int _n_fac = 1;
  for (int i=1; i<=n; i++) { _n_fac *= i; }
  *n_fac = _n_fac;
  return EXIT_SUCCESS;
}


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


//  printf("[ LOG @ %s:%d:%s() ] N_order: %d / deriv_order_arr: (%d,%d,%d)\n", __FILE__,__LINE__,__func__, N_o, deriv_order_arr[0], deriv_order_arr[1], deriv_order_arr[2]);


//    printf("[ LOG @ %s:%d:%s() ] psi_deriv_p_arr : \n", __FILE__,__LINE__,__func__);
//    for (int i_o=0; i_o < N_o; i_o++) {
//      printf("(%7.3f,%7.3f)", std::real(psi_deriv_p_arr[i_o]),std::imag(psi_deriv_p_arr[i_o]));
//    } printf("\n");
  

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

  const int b_vec_matrix_col_num = N_o; // for psi and dpsidx (and more) respectively.
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
//    printf("[ LOG @ %s:%d:%s() ] delta_x_s: %f\n", __FILE__,__LINE__,__func__,delta_x_s);
    for (int i_row = 1; i_row < N_s; ++i_row) {
      power_matrix[i_s][i_row] = power_matrix[i_s][i_row-1] * delta_x_s;
    } // end for loop : `i_row`

  } // for-loop `i_s`

  int n_fac = -1;
  for (int i_order = 0; i_order < b_vec_matrix_col_num; i_order++) {
    deriv_order = deriv_order_arr[i_order];
    if ( factorial(deriv_order, &n_fac) != EXIT_SUCCESS or n_fac <= 0) 
    { return debug_mesg("Failed to evaluate factorial"); }
    b_vec_matrix[i_order][deriv_order] = 1.0 * n_fac; 
//    printf("n_fac: %d / deriv_order: %d\n",n_fac, deriv_order);
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
//      printf("psi_arr[i_x_s_arr[i_s]]: %7.5f\n", std::real(psi_arr[i_x_s_arr[i_s]]));
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
int eval_psi_and_di_psi_and_dj_di_psi_for_sph_harm_basis(
    const int N_s, const int N_rho, const int N_lm, 
    const double r_p_vec[DIM_R], 
    const std::complex<double> **psi_in_sph_harm_basis_arr,
    const double *rho_arr, const int *l_arr, const int *m_arr, 
    const double *rho_p_lim, 
    z_t *psi_p_p, z_t di_psi[DIM_R], z_t dj_di_psi[DIM_R][DIM_R])
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
//  std::complex<double> 
//    psi_p_Y = 0, dpsi_p_Y = 0, psi_p_dtheta_Y = 0, psi_p_dphi_Y = 0,
//    d2psi_p_Y = 0, dpsi_p_dtheta_Y = 0, dpsi_p_dphi_Y = 0, psi_p_d2theta_Y = 0,
//    psi_p_dphi_dtheta_Y = 0, psi_p_d2phi_Y = 0;
//    psi_p_dtheta_dphi_Y = 0, 

  double rho_p = r_p_vec[0], theta_p = r_p_vec[1], phi_p = r_p_vec[2];
  bool out_of_range;

//    printf("[ LOG @ %s:%d:%s() ] rho_p : %f\n", __FILE__,__LINE__,__func__, rho_p);

  //// Check singularity
  if (rho_p == 0 or sin(theta_p) == 0) 
  { return debug_mesg("Zero 'rho_p' or 'sin(theta_p)' occurred"); }

  //// If out-of-range, the velocity is returned as zero
  out_of_range = rho_p_min > rho_p or rho_p >= rho_p_max; 
  if (out_of_range) { 
    for (int i_p_dim = 0; i_p_dim < DIM_R; i_p_dim++) 
    { di_psi[i_p_dim] = 0.0; }
    return EXIT_SUCCESS;
  }
  
  //// Define useful variables
  const double cos_theta_p = cos(theta_p);
//  const double sin_theta_p = sin(theta_p);
//  const double cot_theta_p = cos_theta_p / sin_theta_p;

  //// Aliasing
  const std::complex<double> **psi_arr_arr = psi_in_sph_harm_basis_arr;

  //// Determine the number of derivative orders to be evaluated
  const bool eval_jac = (dj_di_psi != NULL);
  const int N_order_max = 3;
  int N_order;
  if (eval_jac) { N_order = N_order_max; }
  else { N_order = 2; }
  int deriv_order_arr[N_order];
  deriv_order_arr[0] = 0, deriv_order_arr[1] = 1;
  if (eval_jac) { deriv_order_arr[2] = 2; }

  

//  printf("[ LOG @ %s:%d:%s() ] N_order: %d / N_order_max: %d / deriv_order_arr: (%d,%d,%d)\n", __FILE__,__LINE__,__func__, N_order, N_order_max, deriv_order_arr[0], deriv_order_arr[1], deriv_order_arr[2]);



  //// make a container
  z_t F[DIM_R][DIM_R];


  //// Allocate memory
//  std::complex<double> *Ylm_arr, *Yl1m_arr, *dtheta_Ylm_arr;
//  std::complex<double> d2theta_Y;
//  std::complex<double> **psi_deriv_p_lm_arr;
//  std::complex<double> *psi_deriv_p_lm_arr_1d;
//  std::complex<double> **psi_deriv_p_lm_arr_p, **psi_deriv_p_lm_arr_p_max;
//  std::complex<double> *psi_deriv_p_lm_arr_1d_p;
//  psi_deriv_p_lm_arr = new std::complex<double>*[N_lm];
//  psi_deriv_p_lm_arr_1d = new std::complex<double>[N_lm*N_order_max];
    z_t psi_deriv_p_lm_arr[N_order_max];

//  for(psi_deriv_p_lm_arr_p = psi_deriv_p_lm_arr, 
//      psi_deriv_p_lm_arr_p_max = psi_deriv_p_lm_arr + N_lm,
//      psi_deriv_p_lm_arr_1d_p = psi_deriv_p_lm_arr_1d;
//      psi_deriv_p_lm_arr_p < psi_deriv_p_lm_arr_p_max; 
//      psi_deriv_p_lm_arr_p++, psi_deriv_p_lm_arr_1d_p += N_order_max) 
//  { *psi_deriv_p_lm_arr_p = psi_deriv_p_lm_arr_1d_p; }

//  Ylm_arr = new std::complex<double>[N_lm];
//  Yl1m_arr = new std::complex<double>[N_lm];
//  dtheta_Ylm_arr = new std::complex<double>[N_lm];

  
  const size_t gsl_l_max = l_arr[N_lm-1];
  const size_t N_gsl_sphPlm = gsl_sf_legendre_array_n(gsl_l_max); 

  double *sphPlm_arr, *dtheta_sphPlm_arr, *d2theta_sphPlm_arr;
  sphPlm_arr = new double[N_gsl_sphPlm];
  dtheta_sphPlm_arr = new double[N_gsl_sphPlm];
  d2theta_sphPlm_arr = new double[N_gsl_sphPlm];


  return_code = gsl_sf_legendre_deriv2_alt_array(
      GSL_SF_LEGENDRE_SPHARM, gsl_l_max, cos_theta_p, 
      sphPlm_arr, dtheta_sphPlm_arr, d2theta_sphPlm_arr);
  if (return_code != GSL_SUCCESS) 
  { return debug_mesg("Failed during evaluationg of Plm and its derivatives");}
  

  //// Variables to be used in loop
  int i_lm;
  const std::complex<double> **psi_arr_p;
  const int *l_p, *l_p_max, *m_p;
//  std::complex<double> 
//    Ylm, Yl1m, dtheta_Ylm, d2theta_Y;
  std::complex<double> 
    di_psi_lm, dj_di_psi_lm, exp_phi, 
    psi_p, dpsi_p, d2psi_p;
//  std::complex<double> 
//    *p_Ylm, *p_Yl1m, *p_dtheta_Ylm;
  double Plm, dtheta_Plm, d2theta_Plm; // Pl1m;
  double l, m, m_power_of_minus_1;
  int li, mi, m_sign, m_abs;
  bool inconsistent_l_m = true;
  size_t gsl_lm_index; // gsl_l1m_index;


  //// Check whether there is any negative l value
  for (l_p=l_arr, l_p_max=l_arr+N_lm; l_p<l_p_max; l_p++) 
  { if (*l_p < 0) { return debug_mesg("negative l value found"); } }


  //// Initialize
  *psi_p_p = 0.0;
  for (int i=0; i<DIM_R; i++) {
    di_psi[i] = 0.0;
  }
  if (eval_jac) {
    for (int i=0; i<DIM_R; i++) {
      for (int j=0; j<DIM_R; j++)
      { dj_di_psi[i][j] = 0.0; }
    }
  }

  const std::complex<double> imag_unit = std::complex<double>(0.0, 1.0);

  //// Evaluate Ylm etc. and store them into arrays
  for (
      i_lm = 0, psi_arr_p = psi_arr_arr, l_p=l_arr, m_p=m_arr; 
//      p_Ylm=Ylm_arr, p_Yl1m = Yl1m_arr, p_dtheta_Ylm = dtheta_Ylm_arr,
//      psi_deriv_p_lm_arr_1d_p = psi_deriv_p_lm_arr_1d;
      i_lm < N_lm;
      i_lm++, psi_arr_p++, l_p++, m_p++ 
//      p_Ylm++, p_Yl1m++, p_dtheta_Ylm++,
//      psi_deriv_p_lm_arr_1d_p++
      ) 
  {

    li = *l_p, mi = *m_p;
    inconsistent_l_m = (mi > li or mi < -li);
    l = (double) li; m = (double) mi;
    m_sign = 1 - 2 * (1 - (mi>=0));
    m_abs = m_sign * mi;
    m_power_of_minus_1 = 1.0 - 2.0 * (mi%2);


//    for (int i_o=0; i_o < N_order_max; i_o++) {
//      psi_deriv_p_lm_arr[i_o] = 0.0; 
//    }

//    printf("[ LOG @ %s:%d:%s() ] N_order: %d / N_order_max: %d\n", __FILE__,__LINE__,__func__, N_order, N_order_max);
//    printf("[ LOG @ %s:%d:%s() ] psi_deriv_p_lm_arr (before) : \n", __FILE__,__LINE__,__func__);
//    for (int i_o=0; i_o < N_order_max; i_o++) {
//      printf("(%7.3f,%7.3f)", psi_deriv_p_lm_arr[i_o].real(),psi_deriv_p_lm_arr[i_o].imag());
//    } printf("\n");


    //// Estimate `psi_p` and `dpsidrho_p` (and more) using finite difference
    return_code = eval_psi_deriv_p< std::complex<double> >(
        rho_p, *psi_arr_p, rho_arr, N_s, N_rho, rho_p_lim, N_order,
        deriv_order_arr, psi_deriv_p_lm_arr);
    if (return_code != EXIT_SUCCESS)
    { return debug_mesg("Failed to 'eval_psi_and_dpsidx_p()'"); }

 //   printf("[ LOG @ %s:%d:%s() ] psi_deriv_p_lm_arr (after) : \n", __FILE__,__LINE__,__func__);
//    for (int i_o=0; i_o < N_order_max; i_o++) {
//      printf("(%7.3f,%7.3f)", psi_deriv_p_lm_arr[i_o].real(),psi_deriv_p_lm_arr[i_o].imag());
//    } printf("\n");


    //// Deal with the case where the |m|>l, possibly by qprop-dim == 34
    if (inconsistent_l_m) {
      for (int i_order=0; i_order < N_order; i_order++) {
        if (psi_deriv_p_lm_arr[i_order] != 0.0) {
          return debug_mesg(
              "non zero psi_deriv values when l and m are inconsistent");
        }
      }
      continue;
    }

    //// Evaluate ingredients
    exp_phi = std::complex<double>(cos(m*phi_p), sin(m*phi_p));
//    exp_phi = std::complex<double>(cos(m_abs*phi_p), sin(m_abs*phi_p));

//    Plm = gsl_sf_legendre_sphPlm(li, m_abs, cos_theta_p);
    gsl_lm_index = gsl_sf_legendre_array_index(li, m_abs);
//    gsl_l1m_index = gsl_sf_legendre_array_index(li+1, m_abs);
    Plm = sphPlm_arr[gsl_lm_index];
//    Pl1m = dtheta_sphPlm_arr[gsl_l1m_index];

    dtheta_Plm = dtheta_sphPlm_arr[gsl_lm_index];
    d2theta_Plm = d2theta_sphPlm_arr[gsl_lm_index];
//    Pl1m = gsl_sf_legendre_sphPlm(li+1, m_abs, cos_theta_p);

//    Ylm = gsl_sf_legendre_sphPlm(li, m_abs, cos(theta_p)) * exp_phi;
//    Yl1m = gsl_sf_legendre_sphPlm(li+1, m_abs, cos(theta_p)) * exp_phi;
    if (mi < 0) {
      Plm *= m_power_of_minus_1;
//      Pl1m *= m_power_of_minus_1;

      dtheta_Plm *= m_power_of_minus_1;
      d2theta_Plm *= m_power_of_minus_1;
//      Ylm = m_power_of_minus_1 * std::conj(Ylm);
//      Yl1m = m_power_of_minus_1 * std::conj(Yl1m);
    }
    
    //// Store results to arrays
//    *p_Ylm = Ylm;
//    *p_Yl1m = Yl1m;
//    *p_dtheta_Ylm = (
//        - cos_theta_p * (l+1.0) * Ylm
//        + sqrt( (2.0*l+1.0)/(2.0*l+3.0) * (l+1.0+m) * (l+1.0-m) ) * Yl1m
//        ) / sin_theta_p;

//    dtheta_Plm = (
//        - cos_theta_p * (l+1.0) * Plm
//        + sqrt( (2.0*l+1.0)/(2.0*l+3.0) * (l+1.0+m) * (l+1.0-m) ) * Pl1m
//        ) / sin_theta_p;


//    psi_p = psi_deriv_p_lm_arr[0];
//    dpsi_p = psi_deriv_p_lm_arr[1];
//    dtheta_Ylm = *p_dtheta_Ylm;
//    d2theta_Y = 
//      - cot_theta_p * dtheta_Ylm 
//      - ( l*(l+1) - m*m/(sin_theta_p*sin_theta_p) ) * Ylm;
//    d2theta_Plm = 
//      - cot_theta_p * dtheta_Plm 
//      - ( l*(l+1) - m*m/(sin_theta_p*sin_theta_p) ) * Plm;

//    printf("psi_deriv_p_lm_arr @ i_lm=%d: (%5.3f,%5.3f)\n",psi_deriv_p_lm_arr[0].real(), psi_deriv_p_lm_arr[0].imag());
//    printf("Plm: %5.3f / dtheta_Plm: %5.3f\n",Plm, dtheta_Plm);

//    printf("[ LOG @ %s:%d:%s ][i_lm=%03d] psi_deriv_p_lm_arr:\n", __FILE__, __LINE__, __func__, i_lm);
//    print_vec_t_complex(psi_deriv_p_lm_arr);
//  printf("\n");

    F[0][0] = psi_deriv_p_lm_arr[0] / rho_p; 
    F[0][1] = 1.0 / rho_p * ( - F[0][0] + psi_deriv_p_lm_arr[1] );
    F[1][0] = Plm; 
    F[1][1] = dtheta_Plm; 
    F[2][0] = exp_phi; // std::complex<double>(cos(m*phi_p), sin(m*phi_p));
    F[2][1] = F[2][0] * imag_unit * m;
    if (eval_jac) { 
      F[0][2] = 1.0 / rho_p * ( - 2.0 * F[0][1] + psi_deriv_p_lm_arr[2] ); 
      F[1][2] = d2theta_Plm;
      F[2][2] = F[2][1] * imag_unit * m;
    }

//    printf("[ LOG @ %s:%d ][i_lm=%03d] F: \n", __FILE__, __LINE__, i_lm);
//    print_jac_t_complex(F);

    z_t psi_p_lm;

    //// Eval psi
    psi_p_lm = 1.0;
    for (int k=0; k<DIM_R; k++) 
    { psi_p_lm *= F[k][0]; }
    *psi_p_p += psi_p_lm;

//    printf("[ LOG @ %s:%d:%s ][i_lm=%03d] *psi_p_p: %7.5f\n", __FILE__, __LINE__, __func__, i_lm, std::real(*psi_p_p));

    //// Eval di_psi
    for (int i=0; i<DIM_R; i++) {
      di_psi_lm = 1.0;
      for (int k=0; k<DIM_R; k++) 
      { di_psi_lm *= F[k][i==k]; }
      di_psi[i] += di_psi_lm;
    }
//    printf("[ LOG @ %s:%d:%s ][i_lm=%03d] di_psi[i]: %7.5f, %7.5f, %7.5f\n", __FILE__, __LINE__, __func__, i_lm, di_psi[0].real(), di_psi[1].real(), di_psi[2].real());
    
    //// Eval dj_di_psi
    if (eval_jac) {
      for (int i=0; i<DIM_R; i++) {
        for (int j=0; j<DIM_R; j++) {
          dj_di_psi_lm = 1.0;
          for (int k=0; k<DIM_R; k++) 
          { dj_di_psi_lm *= F[k][(i==k)+(j==k)]; }
          dj_di_psi[i][j] += dj_di_psi_lm;
        } // `j`
      } // `i`
    }

  } // end for loop : `i_lm`



  delete [] sphPlm_arr;
  delete [] dtheta_sphPlm_arr;
  delete [] d2theta_sphPlm_arr;



  //// Evaluate summations for velocity 
//  for (
//      i_lm = 0, l_p=l_arr, m_p=m_arr, 
//      p_Ylm = Ylm_arr, p_Yl1m = Yl1m_arr, p_dtheta_Ylm = dtheta_Ylm_arr,
//      psi_deriv_p_lm_arr_1d_p = psi_deriv_p_lm_arr_1d;
//      i_lm < N_lm;
//      i_lm++, l_p++, m_p++, psi_deriv_p_lm_arr_1d_p += N_order_max,
//      p_Ylm++, p_Yl1m++, p_dtheta_Ylm++
//      )
//  {
//    l = (double) *l_p; m = (double) *m_p;
//    Ylm = *p_Ylm; Yl1m = *p_Yl1m; dtheta_Ylm = *p_dtheta_Ylm;
//    psi_p = *(psi_deriv_p_lm_arr_1d_p);
//    dpsi_p = *(psi_deriv_p_lm_arr_1d_p+1);
//
//    //// Add to target summations
//    psi_p_Y += psi_p * Ylm;
//    dpsi_p_Y += dpsi_p * Ylm;
//    psi_p_dtheta_Y += psi_p * dtheta_Ylm;
//    psi_p_dphi_Y += m * psi_p * Ylm;
//    
//  } // i_lm
// 
//  psi_p_dphi_Y *= imag_unit;
//
//
//  //// Evaluate summations for evaluating Jacobian
//  if (eval_jac) {
//    
//    for (
//        i_lm = 0, l_p=l_arr, m_p=m_arr, 
//        p_Ylm = Ylm_arr, p_Yl1m = Yl1m_arr, p_dtheta_Ylm = dtheta_Ylm_arr,
//        psi_deriv_p_lm_arr_1d_p = psi_deriv_p_lm_arr_1d;
//        i_lm < N_lm;
//        i_lm++, l_p++, m_p++, psi_deriv_p_lm_arr_1d_p += N_order_max,
//        p_Ylm++, p_Yl1m++, p_dtheta_Ylm++
//        )
//    {
//      l = (double) *l_p; m = (double) *m_p;
//      Ylm = *p_Ylm; Yl1m = *p_Yl1m; dtheta_Ylm = *p_dtheta_Ylm;
//      d2theta_Y = 
//        - cot_theta_p * dtheta_Ylm 
//        - ( l*(l+1) - m*m/(sin_theta_p*sin_theta_p) ) * Ylm;
//      psi_p = *(psi_deriv_p_lm_arr_1d_p);
//      dpsi_p = *(psi_deriv_p_lm_arr_1d_p+1);
//      d2psi_p = *(psi_deriv_p_lm_arr_1d_p+2);
//  
//      d2psi_p_Y += d2psi_p * Ylm;
//      dpsi_p_dtheta_Y += dpsi_p * dtheta_Ylm;
//      dpsi_p_dphi_Y += dpsi_p * m * Ylm;
//      psi_p_d2theta_Y += psi_p * d2theta_Y;
//      psi_p_dphi_dtheta_Y += psi_p * m * dtheta_Ylm;
//      psi_p_d2phi_Y += psi_p * m * m * Ylm;
//  
//    } // i_lm
//
//  } // end if : `eval_jac`
//
//  dpsi_p_dphi_Y *= imag_unit;
//  psi_p_dphi_dtheta_Y *= imag_unit;
//  psi_p_d2phi_Y *= imag_unit * imag_unit;
//


  //// Deallocate arrays after use
//  delete [] psi_deriv_p_lm_arr_1d;
//  delete [] psi_deriv_p_lm_arr;
//  delete [] Ylm_arr;
//  delete [] Yl1m_arr;
//  delete [] dtheta_Ylm_arr;

//
//#ifdef DEBUG
//
//  //// Print real part
//  printf("real part: \n");
//  printf(
//      "%15s%15s%15s%15s\n", 
//      "denum","numer_p_rho","numer_p_theta","numer_p_phi");
//
//  printf("%15.7f%15.7f%15.7f%15.7f\n",
//      std::real(psi_p_Y), std::real(dpsi_p_Y),
//      std::real(psi_p_dtheta_Y), std::real(psi_p_dphi_Y));
//  printf("\n");
//
//
//  //// Print imaginary part
//  printf("imag part: \n");
//  printf(
//      "%15s%15s%15s%15s\n", 
//      "denum","numer_p_rho","numer_p_theta","numer_p_phi");
//
//  printf("%15.7f%15.7f%15.7f%15.7f\n",
//      std::imag(psi_p_Y), std::imag(dpsi_p_Y),
//      std::imag(psi_p_dtheta_Y), std::imag(psi_p_dphi_Y));
//  printf("\n");
//
//#endif // DEBUG
//
//
//
//#ifdef DEBUG
//  printf("v_p_arr: \n");
//#endif // DEBUG
//      
//
//  // Check singularity
//  if (psi_p_Y == 0.0) { 
//
//    fprintf(stderr,"[ LOG ] %-25s\n","psi_p_Y_arr: ");
//    fprintf(stderr,"(%20s,%20s)\n","real","imag");
//    fprintf(stderr,"(%20.15f,%20.15f)\n",psi_p_Y.real(),psi_p_Y.imag());
//
//    fprintf(stderr,"\n");
//
//    fprintf(stderr,
//        "[ LOG ] r_p_vec\n"
//        "        = (%19s,%19s,%19s)\n"
//        "        = (%19.15f,%19.15f,%19.15f)\n",
//        "rho_p","theta_p","phi_p",
//        r_p_vec[0],r_p_vec[1],r_p_vec[2]);
//
//    fprintf(stderr,"\n");
//
//    fprintf(stderr,
//        "[ LOG ] v_p_vec\n"
//        "        = (%19s,%19s,%19s)\n"
//        "        = (%19.15f,%19.15f,%19.15f)\n",
//        "v_rho_p","v_theta_p","v_phi_p",
//        v_p_vec[0],v_p_vec[1],v_p_vec[2]);
//
//    fprintf(stderr,"\n");
//    
//    const double delta_rho = rho_arr[1] - rho_arr[0];
//    const int i_rho_nls = (int) (r_p_vec[0] - rho_arr[0]) / delta_rho;
//
//    fprintf(stderr,
//        "[ LOG ] radial grid point index of nearest left stencil: "
//        "%d (= %.15f a.u.)\n",
//        i_rho_nls, delta_rho * i_rho_nls);
//    fprintf(stderr,"\n");
//
//    fprintf(stderr,"[ LOG ] psi_arr_arr[i_lm][i_rho_nls]: \n");
//    fprintf(stderr,"\n");
//    fprintf(stderr,"%5s (%19s,%19s)\n","i_lm","real","imag");
//    std::complex<double> _psi;
//    for (int i_lm = 0; i_lm < N_lm; i_lm++) {
//      _psi = psi_arr_arr[i_lm][i_rho_nls];
//      fprintf(stderr,"%5d (%19.15f,%19.15f)\n",
//          i_lm,std::real(_psi),std::imag(_psi));
//    } fprintf(stderr,"\n");
//
//    return debug_mesg("Zero 'psi_p' occurred"); 
//  }  // endif denum == 0
//
//
//  //// Evaluate velocity if not singular
//  v_p_vec[0]
//    = (dpsi_p_Y / psi_p_Y).imag();
//
//  v_p_vec[1]
//    = (psi_p_dtheta_Y / psi_p_Y).imag()
//    / (rho_p);
//
//  v_p_vec[2] 
//    = (psi_p_dphi_Y / psi_p_Y).imag()
//    / (rho_p * sin_theta_p);
//
//
//#ifdef DEBUG
//  printf("%15.5f%15.5f%15.5f\n", 
//      v_p_vec[0], v_p_vec[1], v_p_vec[2]);
//#endif // DEBUG
//
//
//  //// Evaluate Jacobian
//  const std::complex<double> psi_p_Y_sq = psi_p_Y * psi_p_Y;
//  if (eval_jac) {
//
//    jac[0][0] = ( d2psi_p_Y / psi_p_Y - (dpsi_p_Y*dpsi_p_Y)/(psi_p_Y_sq) ).imag();
//    jac[0][1] = ( dpsi_p_dtheta_Y / psi_p_Y - dpsi_p_Y*psi_p_dtheta_Y/(psi_p_Y_sq) ).imag();
//    jac[0][2] = ( dpsi_p_dphi_Y / psi_p_Y - dpsi_p_Y * psi_p_dphi_Y / (psi_p_Y_sq) ).imag();
//
//    jac[1][0] = (-v_p_vec[1] + jac[0][1]) / rho_p;
//    jac[1][1] = ( psi_p_d2theta_Y / psi_p_Y - psi_p_dtheta_Y*psi_p_dtheta_Y/psi_p_Y_sq ).imag() / rho_p;
//    jac[1][2] = ( psi_p_dphi_dtheta_Y / psi_p_Y - psi_p_dtheta_Y*psi_p_dphi_Y / psi_p_Y_sq ).imag() / rho_p;
//
//    jac[2][0] = - 1.0 / rho_p * (v_p_vec[2] - 1.0 / sin_theta_p * jac[0][2]);
//    jac[2][1] = - cot_theta_p * v_p_vec[2] + 1.0 / sin_theta_p * jac[1][2];
//    jac[2][2] = ( psi_p_d2phi_Y / psi_p_Y - psi_p_dphi_Y*psi_p_dphi_Y / psi_p_Y_sq ).imag() / (rho_p * sin_theta_p);
//
//  }

  //// Return from thie program
  return EXIT_SUCCESS;

}



int eval_jac_v_3D(
    const jac_t dj_hi, const vec_t hi, const z_t dj_di_psi[DIM_R][DIM_R],
    const z_t di_psi[DIM_R], const z_t psi, jac_t dj_vi)
{
  for (int i=0; i < DIM_R; i++) {
    for (int j=0; j < DIM_R; j++) {
      dj_vi[i][j] = 1.0 / hi[i] * (
          - dj_hi[i][j] / hi[i] * di_psi[i] / psi 
          + dj_di_psi[i][j] / psi 
          - di_psi[i] / psi * di_psi[j] / psi
          ).imag();
    }
  }
  return EXIT_SUCCESS;
}




int eval_hi_and_dj_hi_for_sph_harm_basis(
    const vec_t r_p_vec, vec_t hi, jac_t dj_hi)
{
  const double rho_p = r_p_vec[0], theta_p = r_p_vec[1];

  //// Evaluate hi
  hi[0] = 1.0; hi[1] = rho_p; hi[2] = rho_p * sin(theta_p);

  //// Evaluate dj_hi
  dj_hi[0][0] = 0.0; dj_hi[0][1] = 0.0; dj_hi[0][2] = 0.0;
  dj_hi[1][0] = 1.0; dj_hi[1][1] = 0.0, dj_hi[1][2] = 0.0;
  dj_hi[2][0] = sin(theta_p); 
  dj_hi[2][1] = rho_p * cos(theta_p); 
  dj_hi[2][2] = 0.0;

  return EXIT_SUCCESS;
}



int eval_v_3D(const vec_t hi, const z_t di_psi[DIM_R], const z_t psi, vec_t vi)
{
  for (int i=0; i<DIM_R; i++) 
  { vi[i] = (di_psi[i] / psi).imag() / hi[i]; }
  return EXIT_SUCCESS;
}



//// Implementation of `eval_v_p_for_sph_harm_basis`
int eval_v_p_for_sph_harm_basis(
    const int N_s, const int N_rho, const int N_lm, 
    const double r_p_vec[DIM_R], 
    const std::complex<double> **psi_in_sph_harm_basis_arr,
    const double *rho_arr, const int *l_arr, const int *m_arr, 
    const double *rho_p_lim, 
    double v_p_vec[DIM_R], double jac[DIM_R][DIM_R])
{

  int return_code = EXIT_FAILURE;
  z_t psi, di_psi[DIM_R], dj_di_psi[DIM_R][DIM_R]; 
  vec_t hi; jac_t dj_hi;

//  printf("[ LOG @ %s:%d:%s() ][!] jac is NULL: '%d'\n", __FILE__, __LINE__, __func__, jac==NULL);

  if (jac == NULL) {
    return_code = eval_psi_and_di_psi_and_dj_di_psi_for_sph_harm_basis(
        N_s, N_rho, N_lm, r_p_vec, psi_in_sph_harm_basis_arr, rho_arr,
        l_arr, m_arr, rho_p_lim, &psi, di_psi, NULL);
  } else {
    return_code = eval_psi_and_di_psi_and_dj_di_psi_for_sph_harm_basis(
        N_s, N_rho, N_lm, r_p_vec, psi_in_sph_harm_basis_arr, rho_arr,
        l_arr, m_arr, rho_p_lim, &psi, di_psi, dj_di_psi);
  }

//  printf("[ LOG @ %s:%d:%s() ][!] psi: (%10.5f,%10.5f)\n", __FILE__, __LINE__, __func__, psi.real(), psi.imag());

//  printf("[ LOG @ %s:%d:%s() ][!] di_psi: \n", __FILE__, __LINE__, __func__);
//  print_vec_t_complex(di_psi);
//  for (int i=0; i<DIM_R; i++) { printf("(%7.5f,%7.5f)", di_psi[i].real(), di_psi[i].imag()); }
//  printf("\n");

//  printf("[ LOG @ %s:%d:%s() ][!] dj_di_psi: \n", __FILE__, __LINE__, __func__);
//  print_jac_t_complex(dj_di_psi);
//  for (int i=0; i<DIM_R; i++) {
//    for (int j=0; j<DIM_R; j++) {
//      printf("(%5.3f,%5.3f)", dj_di_psi[i][j].real(), dj_di_psi[i][j].imag()) ;
//    } printf("\n");
//  } printf("\n");

  if (return_code != EXIT_SUCCESS) 
  { return debug_mesg("Failed evaluating deriv psi and its Jacobian"); }

  return_code = eval_hi_and_dj_hi_for_sph_harm_basis(r_p_vec, hi, dj_hi);
  if (return_code != EXIT_SUCCESS) 
  { return debug_mesg("Failed evaluating hi and dj_hi"); }

//  printf("[ LOG @ %s:%d:%s() ][!] hi: \n", __FILE__, __LINE__, __func__);
//  for (int i=0; i<DIM_R; i++) { printf("%7.5f", hi[i]); }
//  printf("\n");

//  printf("[ LOG @ %s:%d:%s() ][!] dj_hi: \n", __FILE__, __LINE__, __func__);
//  print_jac_t(dj_hi);

  return_code = eval_v_3D(hi, di_psi, psi, v_p_vec);
  if (return_code != EXIT_SUCCESS) 
  { return debug_mesg("Failed evaluating velocity vector"); }
//  printf("[ LOG @ %s:%d:%s() ][!] v_p_vec: \n", __FILE__, __LINE__, __func__);
//  for (int i=0; i<DIM_R; i++) { printf("%7.5f", v_p_vec[i]); }
//  printf("\n");

  if (jac != NULL) {
    return_code = eval_jac_v_3D(dj_hi, hi, dj_di_psi, di_psi, psi, jac);
    if (return_code != EXIT_SUCCESS) 
    { return debug_mesg("Failed evaluating velocity Jacobian"); }
  }
  
//  if (jac == NULL) { printf("jac is NULL\n"); }
//  else { printf("jac is not NULL\n"); }
//  if (jac != NULL) {
//    printf("[ LOG @ %s:%d ] jac: \n", __FILE__, __LINE__);
//    print_jac_t(jac);
//  }

  return EXIT_SUCCESS;

}


//// `eval_v_p_vec_arr_for_sph_harm_basis`
int eval_v_p_vec_arr_for_sph_harm_basis(
    const int N_s, const int N_p, const int N_rho, const int N_lm,
    double **r_p_vec_arr, 
    const std::complex<double> **psi_in_sph_harm_basis_arr,
    const double *rho_arr, const int *l_arr, const int *m_arr, 
    const double *rho_p_lim, 
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
  // `jac_p_arr` : 1D array of `jac_t` type of shape (`N_p`,)
  //             : It is effectivly a 3D array of shape (`N_p`,`DIM_R`,`DIM_R`)
  //

  int return_code = EXIT_FAILURE;

  for (int i_p = 0; i_p < N_p; i_p++) {

    if (jac_p_arr == NULL) { 
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

  } // for-loop : i_p

  return EXIT_SUCCESS;

}





//// `eval_v_p_arr_for_sph_harm_basis`
int eval_v_p_arr_for_sph_harm_basis(
    const int N_s, const int N_p, const int N_rho, const int N_lm,
    double **r_p_arr, const std::complex<double> **psi_in_sph_harm_basis_arr,
    const double *rho_arr, const int *l_arr, const int *m_arr, 
    const double *rho_p_lim, 
    double **v_p_arr, jac_t *jac_p_arr)
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
  // `jac_p_arr` : 1D array of `jac_t` type of shape (`N_p`,)
  //             : It is effectivly a 3D array of shape (`N_p`,`DIM_R`,`DIM_R`)
  //

  int return_code = EXIT_FAILURE;

  double r_p_vec[DIM_R], v_p_vec[DIM_R];

  for (int i_p = 0; i_p < N_p; i_p++) {

    for (int i_dim = 0; i_dim < DIM_R; i_dim++)
    { r_p_vec[i_dim] = r_p_arr[i_dim][i_p]; }

    if (jac_p_arr == NULL) { 
      return_code = eval_v_p_for_sph_harm_basis(
          N_s, N_rho, N_lm, r_p_vec, psi_in_sph_harm_basis_arr,
          rho_arr, l_arr, m_arr, rho_p_lim, v_p_vec);
    } else {
      return_code = eval_v_p_for_sph_harm_basis(
          N_s, N_rho, N_lm, r_p_vec, psi_in_sph_harm_basis_arr,
          rho_arr, l_arr, m_arr, rho_p_lim, v_p_vec, jac_p_arr[i_p]);
    }
    if (return_code != EXIT_SUCCESS) 
    { return debug_mesg("Failed during 'eval_v_p_for_sph_harm_basis()'"); }

    for (int i_dim = 0; i_dim < DIM_R; i_dim++)
    { v_p_arr[i_dim][i_p] = v_p_vec[i_dim]; }

  } // for-loop : i_p

  return EXIT_SUCCESS;

}

