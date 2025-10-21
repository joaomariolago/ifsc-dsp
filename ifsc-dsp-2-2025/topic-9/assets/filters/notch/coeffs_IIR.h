#ifndef INCLUDE_COEFFS_IIR_H_
#define INCLUDE_COEFFS_IIR_H_

#include <stdint.h>

#ifndef NUM_STAGES
#define NUM_STAGES 1
#endif

// --freq 874 --fs 48000 --rho 0.9 --alpha 0.3 --analysis
static const float32_t iirCoeffsF32[NUM_STAGES*5] = {
  +1.0000000000000000e+00f,
  -1.9869254612738059e+00f,
  +1.0000000000000000e+00f,
  +1.8340850411758207e+00f,
  -8.5384615384615392e-01f
};

#endif /* INCLUDE_COEFFS_IIR_H_ */
