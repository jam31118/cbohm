#include "../include/velocity.hh"


int handle_gesv_info(int info) {

  if (info != 0) {

    fprintf(stderr, "[ERROR] Unsuccessful exit from 'dgesv_()'\n");

    if (info < 0) {             
      fprintf(stderr, "[ERROR] [info == '%d'] "
          "the '%d'-th argument had an illegal value\n", info, -info);
    } 
    else if (info > 0) {
      fprintf(stderr, "[ERROR] [info == '%d'] singularity happended\n", info);  
    }
    return EXIT_FAILURE;

  }

  return EXIT_SUCCESS;
}


int return_with_mesg(const char *mesg, int return_code) {
  fprintf(stderr,"[ERROR] %s\n", mesg);
  return return_code;
}


//int eval_psi_and_dpsidx_arr_form(
//    double *x_p_arr, std::complex<double> *psi_arr, double *x_arr, 
//    int N_s, int N_p, int N_x, const double *x_p_lim, 
//    std::complex<double> *psi_p_arr, std::complex<double> *dpsidx_p_arr,
//    ); 

//int eval_psi_and_dpsidx_arr(
//    double *x_p_arr, std::complex<double> *psi_arr, double *x_arr, 
//    int N_s, int N_p, int N_x, const double *x_p_lim, 
//    std::complex<double> *psi_p_arr, std::complex<double> *dpsidx_p_arr); 


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
  const int i_x_max = N_x - 1;

  
  //// Aliasing
  const int num_of_stencils = N_s;
  const int num_of_particles = N_p;


  //// Construct array of indices for a set of stencils per particle
  int **i_p_arr_arr = new int*[num_of_stencils];
  i_p_arr_arr[0] = new int[num_of_stencils * num_of_particles];
//  int *i_p_arr, *i_p_arr_max = i_p_arr_arr + num_of_stencils;
  for (int i_s=0; i_s<N_s; ++i_s) {
    i_p_arr_arr[i_s] = i_p_arr_arr[0] + i_s * num_of_particles; 
  }

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


  const int N_s_on_left = (num_of_stencils / 2) - 1;
  const int index_of_nearest_left_stencil = (num_of_stencils / 2) - 1;
  const int i_nlp = index_of_nearest_left_stencil; // aliasing
  const int N_s_on_right = num_of_stencils - N_s_on_left - 1;

  double *x_p_arr_p, *x_p_arr_p_max = x_p_arr + N_p;
  double *psi_p_arr_p, *dpsidx_p_arr_p;
  double x_p;

  int *i_x_s_arr = new int[num_of_stencils];
  int i_s, *i_x_s_arr_p, *i_x_s_arr_p_max = i_x_s_arr + N_s;
  int i_x_s_at_nls;

  int shift_offset_to_right, shift_offset_to_left;
  bool do_shift_to_right, do_shift_to_left;

  double **power_matrix = new double*[num_of_stencils];
  double *power_matrix_1d = new double[num_of_stencils*num_of_stencils];
  for (int i_s = 0; i_s < num_of_stencils; ++i_s) {
    power_matrix[i_s] = power_matrix_1d + i_s * num_of_stencils;
  }
  int b_vec_matrix_col_num = 2; // for psi and dpsidx respectively.
  double **b_vec_matrix = new double*[b_vec_matrix_col_num];
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
  int *pivot_indices = new int[num_of_stencils];
  int gesv_info;



  for (
      x_p_arr_p = x_p_arr, 
      psi_p_arr_p = psi_p_arr, dpsidx_p_arr_p = dpsidx_p_arr;
      x_p_arr_p < x_p_arr_p_max; 
      ++x_p_arr_p, ++psi_p_arr_p, ++dpsidx_p_arr_p
      )
  {

    x_p = *x_p_arr_p;

    //// Check in-range
    assert(x_p_lim[0] <= x_p and x_p < x_p_lim[1]);

  
    i_x_s_at_nls = (x_p - x_arr[0]) / delta_x;

    shift_offset_to_right = N_s_on_left - i_x_s_at_nls;
    shift_offset_to_left = (i_x_max - i_x_s_at_nls) - N_s_on_right;

    do_shift_to_right = shift_offset_to_right > 0;
    do_shift_to_left = shift_offset_to_left < 0;
    if (do_shift_to_right and do_shift_to_left) {
      return return_with_mesg("Cannot shift to both direction");
    }


    for (
        i_s = 0, i_x_s_arr_p = i_x_s_arr;
        i_x_s_arr_p < i_x_s_arr_p_max;
        ++i_x_s_arr_p, ++i_s
        ) 
    {
      *i_x_s_arr_p = i_x_s_at_nls + (i_s - i_nlp) 
        + do_shift_to_right * shift_offset_to_right
        + do_shift_to_left * shift_offset_to_left;

      
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
      double psi_x_s = psi_arr[i_x_s_arr[i_s]];
      *psi_p_arr_p += coef_vec_matrix[0][i_s] * psi_x_s;
      *dpsidx_p_arr_p += coef_vec_matrix[1][i_s] * psi_x_s;

#ifdef DEBUG
      printf(
          "psi_x_s: %7.3f,"
          " coef_vec_matrix[0][i_s]: %7.3f,"
          " coef_vec_matrix[1][i_s]: %7.3f,"
          " *psi_p_arr_p: %7.3f, *dpsidx_p_arr_p: %7.3f\n",
          psi_x_s, 
          coef_vec_matrix[0][i_s], 
          coef_vec_matrix[1][i_s], 
          *psi_p_arr_p, *dpsidx_p_arr_p
      );
#endif // DEBUG

    }
  
  } // for-loop `x_p_arr_p`


  //// Deallocate memory
  delete [] i_x_s_arr; 

  delete [] power_matrix[0];
  delete [] power_matrix;
  delete [] b_vec_matrix[0];
  delete [] b_vec_matrix;


  // Returns if everything is fine.
  return EXIT_SUCCESS;
}



//// Instantiation
template int eval_psi_and_dpsidx_arr<double>(
    double *x_p_arr, double *psi_arr, double *x_arr, 
    int N_s, int N_p, int N_x, const double *x_p_lim, 
    double *psi_p_arr, double *dpsidx_p_arr);



