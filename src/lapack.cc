#include "../include/lapack.hh"
#include "../include/log.hh"

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


