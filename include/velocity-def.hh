#ifndef _VELOCITY_CONSTANT_HH_
#define _VELOCITY_CONSTANT_HH_

#include <complex>

#define CODE_OUT_OF_RANGE 256
#define DIM_R 3

typedef double jac_t[DIM_R][DIM_R];
typedef double vec_t[DIM_R];

typedef std::complex<double> z_t;

typedef z_t jac_z_t[DIM_R][DIM_R];
typedef z_t vec_z_t[DIM_R];

#endif // _VELOCITY_CONSTANT_HH_
