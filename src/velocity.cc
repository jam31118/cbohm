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
//  { return return_with_mesg("Cannot shift to both direction"); }
  { return debug_mesg("Cannot shift to both direction"); }
//  { return return_with_debug_mesg(
//      "Cannot shift to both direction", __FILE__, __LINE__, __func__); }

  //// Evaluate the offset
  *offset 
    = do_shift_to_right * shift_offset_to_right 
    + do_shift_to_left * shift_offset_to_left;

  //// Return if everything is fine
  return EXIT_SUCCESS;

}


template <typename T>
int eval_psi_and_dpsidx_arr(
    double *x_p_arr, T *psi_arr, double *x_arr, 
    int N_s, int N_p, int N_x, const double *x_p_lim, 
    T *psi_p_arr, T *dpsidx_p_arr) {

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


  //// Define some useful variables
  const double delta_x = x_arr[1] - x_arr[0];
//  const int i_x_max = N_x - 1;

  
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

  double *x_p_arr_p, *x_p_arr_p_max = x_p_arr + N_p;
  T *psi_p_arr_p, *dpsidx_p_arr_p;
  double x_p;

  int i_x_s_arr[num_of_stencils];
//  int *i_x_s_arr = new int[num_of_stencils];
  int i_s, *i_x_s_arr_p, *i_x_s_arr_p_max = i_x_s_arr + N_s;
  int i_x_s_at_nls;

//  const int N_s_on_left = (num_of_stencils / 2) - 1;
//  const int N_s_on_right = num_of_stencils - N_s_on_left - 1;
//  int shift_offset_to_right, shift_offset_to_left, shift_offset;
//  bool do_shift_to_right, do_shift_to_left;
  int shift_offset;

  double *power_matrix[num_of_stencils];
//  double **power_matrix = new double*[num_of_stencils];
  double *power_matrix_1d = new double[num_of_stencils*num_of_stencils];
  for (int i_s = 0; i_s < num_of_stencils; ++i_s) {
    power_matrix[i_s] = power_matrix_1d + i_s * num_of_stencils;
  }
  int b_vec_matrix_col_num = 2; // for psi and dpsidx respectively.
  double *b_vec_matrix[b_vec_matrix_col_num];
//  double **b_vec_matrix = new double*[b_vec_matrix_col_num];
  double *b_vec_matrix_1d = new double[num_of_stencils * b_vec_matrix_col_num];
  for (int i_col = 0; i_col < b_vec_matrix_col_num; i_col++) {
    b_vec_matrix[i_col] = b_vec_matrix_1d + i_col * num_of_stencils;
  }


  // An alias for `b_vec_matrix`
  // Note that the `b_vec_matrix` holds the coefficients after gesv() routine
  // This `coef_vec_matrix` will be used after the gesv() routine.
  double **coef_vec_matrix = b_vec_matrix; 
  
  double delta_x_s;

  // The pivot indices that define the permutation matrix P;
  // row i of the matrix was interchanged with row IPIV(i).
  int pivot_indices[num_of_stencils];
//  int *pivot_indices = new int[num_of_stencils];
  int gesv_info;

  T psi_x_s;


  for (
      x_p_arr_p = x_p_arr, 
      psi_p_arr_p = psi_p_arr, dpsidx_p_arr_p = dpsidx_p_arr;
      x_p_arr_p < x_p_arr_p_max; 
      ++x_p_arr_p, ++psi_p_arr_p, ++dpsidx_p_arr_p
      )
  {

    //// Evaluate `x_p` for preventing repetitive 
    x_p = *x_p_arr_p;

    //// Check in-range
    if (x_p_lim[0] > x_p or x_p >= x_p_lim[1]) { continue; }
  

    i_x_s_at_nls = (x_p - x_arr[0]) / delta_x;

    if (eval_shift_offset(N_s,i_x_s_at_nls,N_x,&shift_offset) != EXIT_SUCCESS)
    { return debug_mesg("Failed to evaluate 'shift_offset'"); }
//    shift_offset_to_right = N_s_on_left - i_x_s_at_nls;
//    shift_offset_to_left = (i_x_max - i_x_s_at_nls) - N_s_on_right;
//
//    do_shift_to_right = shift_offset_to_right > 0;
//    do_shift_to_left = shift_offset_to_left < 0;
//    if (do_shift_to_right and do_shift_to_left) {
//      return return_with_mesg("Cannot shift to both direction");
//    }
//
//    shift_offset 
//      = do_shift_to_right * shift_offset_to_right 
//      + do_shift_to_left * shift_offset_to_left;



    for (
        i_s = 0, i_x_s_arr_p = i_x_s_arr;
        i_x_s_arr_p < i_x_s_arr_p_max;
        ++i_x_s_arr_p, ++i_s
        ) 
    {
      *i_x_s_arr_p = i_x_s_at_nls + (i_s - i_nlp) + shift_offset; 
//        + do_shift_to_right * shift_offset_to_right
//        + do_shift_to_left * shift_offset_to_left;

      
      // Construct the power matrix - in column major for FORTRAN LAPACK
      power_matrix[i_s][0] = 1.0; 
      delta_x_s = x_arr[*i_x_s_arr_p] - x_p;
      for (int i_row = 1; i_row < N_s; ++i_row) {
        power_matrix[i_s][i_row] = power_matrix[i_s][i_row-1] * delta_x_s;
      }

      b_vec_matrix[0][i_s] = 0.0;
      b_vec_matrix[1][i_s] = 0.0;

    } // for-loop `i_s`



    b_vec_matrix[0][0] = 1.0; // for psi
    b_vec_matrix[1][1] = 1.0; // for dpsidx


 
#ifdef DEBUG
    printf("x_p = %7.3f\n",x_p);
    for (int i_row = 0; i_row < N_s; i_row++) {
      printf("%7.3f",x_arr[i_x_s_arr[i_row]]);
      for (int i_col = 0; i_col < N_s; i_col++) {
        printf("%7.3f", power_matrix[i_col][i_row]); 
      } printf("%7.3f%7.3f", b_vec_matrix[0][i_row], b_vec_matrix[1][i_row]);
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
    if ( handle_gesv_info(gesv_info) != EXIT_SUCCESS) {
      return return_with_mesg("Failed solving for coeffcients");
    }


#ifdef DEBUG

    printf("power_matrix_1d: ");
    for (int i=0; i < N_s*N_s; i++) {
      printf("%7.3f", power_matrix_1d[i]);
    } printf("\n");
    printf("coef_vec_matrix: ");
    for (int i=0; i < N_s*b_vec_matrix_col_num; i++) {
      printf("%7.3f", coef_vec_matrix[0][i]);
    } printf("\n");

#endif // DEBUG
    

    //// Evaluate the finite-difference-approximated values: psi, dpsidx
    *psi_p_arr_p = 0, *dpsidx_p_arr_p = 0;
    for (int i_s = 0; i_s < N_s; i_s++) {
      psi_x_s = psi_arr[i_x_s_arr[i_s]];
      *psi_p_arr_p += coef_vec_matrix[0][i_s] * psi_x_s;
      *dpsidx_p_arr_p += coef_vec_matrix[1][i_s] * psi_x_s;

#ifdef DEBUG
      printf(
          "psi_x_s: %7.3f,"
          " coef_vec_matrix[0][i_s]: %7.3f,"
          " coef_vec_matrix[1][i_s]: %7.3f,"
          " *psi_p_arr_p: %7.3f, *dpsidx_p_arr_p: %7.3f\n",
          std::real(psi_x_s), 
          coef_vec_matrix[0][i_s], 
          coef_vec_matrix[1][i_s], 
          std::real(*psi_p_arr_p), std::real(*dpsidx_p_arr_p)
      );
#endif // DEBUG

    }
  
  } // for-loop `x_p_arr_p`


  //// Deallocate memory
//  delete [] i_x_s_arr; 

  delete [] power_matrix[0];
//  delete [] power_matrix;
  delete [] b_vec_matrix[0];
//  delete [] b_vec_matrix;


  // Returns if everything is fine.
  return EXIT_SUCCESS;
}


//// Instantiation
template int eval_psi_and_dpsidx_arr<double>(
    double *x_p_arr, double *psi_arr, double *x_arr, 
    int N_s, int N_p, int N_x, const double *x_p_lim, 
    double *psi_p_arr, double *dpsidx_p_arr);

template int eval_psi_and_dpsidx_arr< std::complex<double> >(
    double *x_p_arr, std::complex<double> *psi_arr, double *x_arr, 
    int N_s, int N_p, int N_x, const double *x_p_lim, 
    std::complex<double> *psi_p_arr, std::complex<double> *dpsidx_p_arr);



//// `eval_v_p_arr_for_sph_harm_basis`
int eval_v_p_arr_for_sph_harm_basis(
    const int N_s, const int N_p, const int N_r_dim, 
    const int N_rho, const int N_lm,
    double **r_p_arr, std::complex<double> **psi_in_sph_harm_basis_arr,
    double *rho_arr, int *l_arr, int *m_arr, const double *rho_p_lim, 
    double **v_p_arr)
{

  //// Function Argument
  //
  // [INPUT]
  //
  // `N_s` : the number of stencils used in finite-difference approximation
  //         for estimating psi and dpsidrho at particle's radial coodinate
  // `N_p` : the number of particles
  // `N_r_dim` : the dimension of paritlce's position vector space
  // `N_rho` : the number of radial grid points
  // `N_lm` : the number of spherical harmonics basis
  //
  // `r_p_arr` : 2D array of shape (`N_r_dim`,`N_p`)
  // `psi_in_sph_harm_basis_arr` : 2D array of shape (`N_lm`,`N_rho`)
  // `rho_arr` : 1D array of shape (`N_rho`,)
  // `l_arr` : 1D array of shape (`l_arr`,)
  // `m_arr` : 1D array of shape (`m_arr`,)
  // `rho_p_lim` : 1D array of shape (2,)
  //
  // [OUTPUT]
  //
  // `v_p_arr` : 2D array of shape (`N_r_dim`,`N_p`)
  //


  //// Check arguments
  assert(N_s > 0 and N_p > 0 and N_r_dim > 0 and N_rho > N_s and N_lm > 0);


  //// Define useful variables
  const double rho_p_min = rho_p_lim[0], rho_p_max = rho_p_lim[1];


  //// Aliasing
  std::complex<double> **psi_arr_arr = psi_in_sph_harm_basis_arr;


  //// Allocate data storage
  std::complex<double> *Ylm_p_arr = new std::complex<double>[N_p];
  std::complex<double> *psi_p_arr = new std::complex<double>[N_p];
  std::complex<double> *dpsidrho_p_arr = new std::complex<double>[N_p];

  const size_t type_size = sizeof(std::complex<double>);
  std::complex<double> 
    *denum_p_arr = (std::complex<double> *) std::calloc(N_p, type_size),
    *numer_p_rho_arr = (std::complex<double> *) std::calloc(N_p, type_size),
    *numer_p_theta_arr = (std::complex<double> *) std::calloc(N_p, type_size),
    *numer_p_phi_arr = (std::complex<double> *) std::calloc(N_p, type_size);
  

  double 
    *rho_p_arr = r_p_arr[0], 
    *theta_p_arr = r_p_arr[1],
    *phi_p_arr = r_p_arr[2];

  int return_code = EXIT_FAILURE;
  

  //// Variables to be used in loop
  int i_lm;
  std::complex<double> *psi_arr, **psi_arr_p;
  int *l_p, *m_p;
  std::complex<double> exp_phi, Ylm, Yl1m, psi_p, dpsi_p;
  double l, m, m_power_of_minus_1;
  int li, mi, m_sign, m_abs;

  double *theta_p_p, *phi_p_p, rho_p, theta_p, phi_p;
  int i_p;

  bool out_of_range;

  std::complex<double> *psi_p_p, *dpsi_p_p;


#ifdef DEBUG
  printf("rho_p_arr: \n");
  for (int i_p = 0; i_p < N_p; i_p++) {
    printf("%7.3f",r_p_arr[0][i_p]);
  } printf("\n");
#endif // DEBUG


  for (
      i_lm = 0, psi_arr_p = psi_arr_arr, l_p=l_arr, m_p=m_arr;
      i_lm < N_lm;
      i_lm++, psi_arr_p++, l_p++, m_p++
      ) 
  {

    li = *l_p, mi = *m_p;
    l = (double) li;
    m = (double) mi;
    m_sign = 1 - 2 * (1 - (mi>=0));
    m_power_of_minus_1 = 1.0 - 2.0 * (mi%2);
    m_abs = m_sign * mi;
    psi_arr = *psi_arr_p;

  
    //// Estimate `psi_p` and `dpsidrho_p` using finite difference
    return_code = eval_psi_and_dpsidx_arr< std::complex<double> >(
        rho_p_arr, psi_arr, rho_arr,
        N_s, N_p, N_rho, rho_p_lim,
        psi_p_arr, dpsidrho_p_arr);
    if (return_code != EXIT_SUCCESS) {
      return return_with_mesg("Failed to 'eval_psi_and_dpsidx_arr()'");
    }


#ifdef DEBUG
    printf("psi_arr: \n");
    for (int i_p = 0; i_p < N_p; i_p++) {
      printf(
          "%7.3f%7.3f\n", 
          psi_p_arr[i_p].real(), dpsidrho_p_arr[i_p].real());
    } printf("\n");
#endif // DEBUG


    for (
        i_p = 0, theta_p_p = theta_p_arr, phi_p_p = phi_p_arr,
        psi_p_p = psi_p_arr, dpsi_p_p = dpsidrho_p_arr;
        i_p < N_p;
        i_p++, theta_p_p++, phi_p_p++, psi_p_p++, dpsi_p_p++
        )
    {

      rho_p = rho_p_arr[i_p];

      out_of_range = rho_p_min > rho_p or rho_p >= rho_p_max; 
      if (out_of_range) { continue; }

      theta_p = *theta_p_p, phi_p = *phi_p_p;
      psi_p = *psi_p_p, dpsi_p = *dpsi_p_p;

      if (rho_p == 0 or sin(theta_p) == 0) {
        return return_with_mesg("Zero 'rho_p' or 'sin(theta_p)' occurred");
      }

      exp_phi = std::complex<double>(cos(m_abs*phi_p), sin(m_abs*phi_p));
      Ylm = gsl_sf_legendre_sphPlm(li, m_abs, cos(theta_p)) * exp_phi;
      Yl1m = gsl_sf_legendre_sphPlm(li+1, m_abs, cos(theta_p)) * exp_phi;
      if (mi < 0) {
        Ylm = m_power_of_minus_1 * std::conj(Ylm);
        Yl1m = m_power_of_minus_1 * std::conj(Yl1m);
      }

      denum_p_arr[i_p] += psi_p * Ylm;
      numer_p_rho_arr[i_p] += dpsi_p * Ylm;
      numer_p_theta_arr[i_p] 
        += -cos(theta_p) * l * psi_p * Ylm
        + sqrt( (2.0*l+1.0)/(2.0*l+3.0) * (l+1.0+m) * (l+1.0-m) )
          * psi_p * Yl1m;
      numer_p_phi_arr[i_p] += m * psi_p * Ylm;
      
    } // i_p

  } // i_lm



#ifdef DEBUG


  //// Print real part
  printf("real part: \n");
  printf(
      "%15s%15s%15s%15s\n", 
      "denum","numer_p_rho","numer_p_theta","numer_p_phi");

  for (int i_p = 0; i_p < N_p; i_p++) {
    printf("%15.7f%15.7f%15.7f%15.7f\n",
        std::real(denum_p_arr[i_p]), std::real(numer_p_rho_arr[i_p]),
        std::real(numer_p_theta_arr[i_p]), std::real(numer_p_phi_arr[i_p]));
  } printf("\n");


  //// Print imaginary part
  printf("imag part: \n");
  printf(
      "%15s%15s%15s%15s\n", 
      "denum","numer_p_rho","numer_p_theta","numer_p_phi");

  for (int i_p = 0; i_p < N_p; i_p++) {
    printf("%15.7f%15.7f%15.7f%15.7f\n",
        std::imag(denum_p_arr[i_p]), std::imag(numer_p_rho_arr[i_p]),
        std::imag(numer_p_theta_arr[i_p]), std::imag(numer_p_phi_arr[i_p]));
  } printf("\n");


#endif // DEBUG




#ifdef DEBUG
  printf("v_p_arr: \n");
#endif // DEBUG

  for (int i_p = 0; i_p < N_p; i_p++)
  {
    
    //// Check in-range
    out_of_range = rho_p_min > rho_p_arr[i_p] or rho_p_arr[i_p] >= rho_p_max;
    
    //// Velocity evaluation
    if (out_of_range) { // then velocity is set to be zero

      for (int i_r_dim = 0; i_r_dim < N_r_dim; i_r_dim++) {
        v_p_arr[i_r_dim][i_p] = 0.0;
      }
    } 
    else {
      
      // Check singularity
      if (denum_p_arr[i_p] == 0.0) { 
  
        fprintf(stderr,"[ LOG ] %-25s\n","denum_p_arr: ");
        fprintf(stderr,"%5s (%20s,%20s)\n","i_p","real","imag");
        for (int i_p = 0; i_p < N_p; i_p++) {
          fprintf(stderr,"%5d (%20.15f,%20.15f)\n",
              i_p, denum_p_arr[i_p].real(),denum_p_arr[i_p].imag());
        } fprintf(stderr,"\n");
  
        fprintf(stderr,"[ LOG ] index of particle with zero psi_p: %d\n",i_p);
  
        fprintf(stderr,
            "[ LOG ] r_p_arr\n"
            "        = (%19s,%19s,%19s)\n"
            "        = (%19.15f,%19.15f,%19.15f)\n",
            "rho_p","theta_p","phi_p",
            r_p_arr[0][i_p],r_p_arr[1][i_p],r_p_arr[2][i_p]);
  
        fprintf(stderr,"\n");
  
        fprintf(stderr,
            "[ LOG ] v_p_arr\n"
            "        = (%19s,%19s,%19s)\n"
            "        = (%19.15f,%19.15f,%19.15f)\n",
            "v_rho_p","v_theta_p","v_phi_p",
            v_p_arr[0][i_p],v_p_arr[1][i_p],v_p_arr[2][i_p]);
  
        fprintf(stderr,"\n");
        
        double delta_rho = rho_arr[1] - rho_arr[0];
        int i_rho_nls = (int) (r_p_arr[0][i_p] - rho_arr[0]) / delta_rho;
  
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
  
        return return_with_mesg("Zero 'psi_p' occurred"); 
      }


      // Evaluate velocity if not singular
      v_p_arr[0][i_p] 
        = (numer_p_rho_arr[i_p] / denum_p_arr[i_p]).imag();

      v_p_arr[1][i_p] 
        = (numer_p_theta_arr[i_p] / denum_p_arr[i_p]).imag()
        / (rho_p_arr[i_p] * sin(theta_p_arr[i_p]));

      v_p_arr[2][i_p] 
        = (numer_p_phi_arr[i_p] / denum_p_arr[i_p]).real()
        / (rho_p_arr[i_p] * sin(theta_p_arr[i_p]));

    }

#ifdef DEBUG
    printf("%5d%15.5f%15.5f%15.5f\n", 
        i_p, 
        v_p_arr[0][i_p], v_p_arr[1][i_p], v_p_arr[2][i_p]);
#endif // DEBUG

  }


  //// Deallocate data storage
  delete [] Ylm_p_arr;
  delete [] psi_p_arr;
  delete [] dpsidrho_p_arr;

  std::free(denum_p_arr);
  std::free(numer_p_rho_arr);
  std::free(numer_p_theta_arr);
  std::free(numer_p_phi_arr);


  //// Return from thie program
  return EXIT_SUCCESS;

}

