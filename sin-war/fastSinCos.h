
#ifndef FAST_SIN_COS_H_
#define FAST_SIN_COS_H_

#include <cmath>
#include "sin_lut.h"


namespace lut {
   const double PI = 3.14159265358979323846;
   const double HALF_PI = 0.5 * PI;
   const double THREE_HALF_PI = 3.0 * HALF_PI;
   const double TWO_PI = 2.0 * PI;
   const double TWO_PI_INV = 1.0 / TWO_PI;

   inline double fast_sin(double x)
   {

      x = fmod(x, TWO_PI); // 0 <= x <= 2pi
      if (x < HALF_PI)
         return sin_lut[static_cast<unsigned>(x / lut_step)];
      else if (x < PI)
         return sin_lut[static_cast<unsigned>((PI - x) / lut_step)];
      else if (x < THREE_HALF_PI)
         return -sin_lut[static_cast<unsigned>((x - PI) / lut_step)];
      else
         return -sin_lut[static_cast<unsigned>((TWO_PI - x) / lut_step)];
   }

   inline double fast_cos(double x)
   {
      return fast_sin(x + HALF_PI);
   }

}
#endif /* FAST_SIN_COS_H_ */
