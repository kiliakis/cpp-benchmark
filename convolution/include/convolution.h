///////////////////////////////////////////////////////////////////////////////
// convolution.h
// =============
// convolution 1D and 2D
//
//  AUTHOR: Song Ho Ahn
// CREATED: 2005-07-18
// UPDATED: 2005-09-08
//
// Copyright (c) 2005 Song Ho Ahn
///////////////////////////////////////////////////////////////////////////////

#ifndef CONVOLUTION_H
#define CONVOLUTION_H


// linear convolution function
// @a: first vector
// @b: second vector
// @return convolution of a and b
static inline std::vector<double> convolution1(const std::vector<double> &a,
      const std::vector<double> &b)
{
   std::vector<double> res;
   res.resize(a.size() + b.size() - 1);

   #pragma omp parallel for
   for (uint i = 0; i < res.size(); ++i) {
      uint i1 = i;
      double temp = 0;
      for (unsigned int j = 0; j < b.size(); ++j) {
         if (i1 >= 0 && i1 < a.size()) {
            temp += a[i1] * b[j];
         }
         i1--;
         res[i] = temp;
      }
   }
   return res;
}


// 1D convolution /////////////////////////////////////////////////////////////
// We assume input and kernel signal start from t=0. (The first element of
// kernel and input signal is at t=0)
// it returns false if parameters are not valid.
// 1D convolution
// We assume input and kernel signal start from t=0.
///////////////////////////////////////////////////////////////////////////////
bool convolution2(double *in, double *out, int dataSize, double *kernel, int kernelSize)
{
   int i, j, k;

   // check validity of params
   if (!in || !out || !kernel) return false;
   if (dataSize <= 0 || kernelSize <= 0) return false;

   // start convolution from out[kernelSize-1] to out[dataSize-1] (last)
   for (i = kernelSize - 1; i < dataSize; ++i) {
      out[i] = 0;                             // init to 0 before accumulate

      for (j = i, k = 0; k < kernelSize; --j, ++k)
         out[i] += in[j] * kernel[k];
   }

   // convolution from out[0] to out[kernelSize-2]
   for (i = 0; i < kernelSize - 1; ++i) {
      out[i] = 0;                             // init to 0 before sum

      for (j = i, k = 0; j >= 0; --j, ++k)
         out[i] += in[j] * kernel[k];
   }

   return true;
}


//convolution algorithm
double *convolution3(double *A, double *B, int lenA, int lenB, int *lenC)
{
   int nconv;
   int i, j, i1;
   double tmp;
   double *C;

   //allocated convolution array
   nconv = lenA + lenB - 1;
   C = (double *) calloc(nconv, sizeof(double));

   //convolution process
   for (i = 0; i < nconv; i++) {
      i1 = i;
      tmp = 0.0;
      for (j = 0; j < lenB; j++) {
         if (i1 >= 0 && i1 < lenA)
            tmp = tmp + (A[i1] * B[j]);

         i1 = i1 - 1;
         C[i] = tmp;
      }
   }

   //get length of convolution array
   (*lenC) = nconv;

   //return convolution array
   return (C);
}


double *convolution4(const double *__restrict__ Signal,
                     const uint SignalLen,
                     const double *__restrict__ Kernel,
                     const uint KernelLen)
{
   const uint size = SignalLen + KernelLen - 1;
   double *__restrict__ Result = (double *) malloc(
                                    size * sizeof(double));

   #pragma omp parallel for
   for (uint n = 0; n < size; ++n) {
      //Result[n] = 0;
      const uint kmin = (n >= KernelLen - 1) ? n - (KernelLen - 1) : 0;
      const uint kmax = (n < SignalLen - 1) ? n : SignalLen - 1;

      for (uint k = kmin; k <= kmax; ++k) {
         Result[n] += Signal[k] * Kernel[n - k];
      }
   }
   return Result;
}

/*
void convolve(Workspace &ws, double *src, double *kernel)
{
   double temp;
   int i, j, k, l;
   int low_k, high_k, low_l, high_l;
   int i_src, j_src;

   if (ws.h_dst <= 0 || ws.w_dst <= 0)
      return;

   switch (ws.mode) {
      case LINEAR_FULL:
         // Full linear convolution of size N + M -1
         for (i = 0 ; i < ws.h_dst ; ++i) {
            low_k = std::max(0, i - ws.h_kernel + 1);
            high_k = std::min(ws.h_src - 1, i);
            for (j = 0 ; j < ws.w_dst ; ++j) {
               low_l = std::max(0, j - ws.w_kernel + 1);
               high_l = std::min(ws.w_src - 1 , j);
               temp = 0.0;
               for (k = low_k ; k <= high_k ; ++k) {
                  for (l = low_l ; l <= high_l ; ++l) {
                     temp += src[k * ws.w_src + l] * kernel[(i - k) * ws.w_kernel + (j - l)];
                  }
               }
               ws.dst[i * ws.w_dst + j] = temp;
            }
         }
         break;
      case LINEAR_SAME:
         // Same linear convolution, of size N
         for (i = 0 ; i < ws.h_dst ; ++i) {
            low_k = std::max(0, i - int((ws.h_kernel - 1.0) / 2.0));
            high_k = std::min(ws.h_src - 1, i + int(ws.h_kernel / 2.0));
            for (j = 0 ; j < ws.w_dst ; ++j) {
               low_l = std::max(0, j - int((ws.w_kernel - 1.0) / 2.0));
               high_l = std::min(ws.w_src - 1, j + int(ws.w_kernel / 2.0));
               temp = 0.0;
               for (k = low_k ; k <= high_k ; ++k) {
                  for (l = low_l ; l <= high_l ; ++l) {
                     temp += src[k * ws.w_src + l] * kernel[(i - k + int(ws.h_kernel / 2.0)) * ws.w_kernel + (j - l + int(ws.w_kernel / 2.0))];
                  }
               }
               ws.dst[i * ws.w_dst + j] = temp;
            }
         }
         break;
      case LINEAR_VALID:
         // Valid linear convolution, of size N - M
         for (i = 0 ; i < ws.h_dst ; ++i) {
            for (j = 0 ; j < ws.w_dst ; ++j) {
               temp = 0.0;
               for (k = i ; k <= i + ws.h_kernel - 1; ++k) {
                  for (l = j ; l <= j + ws.w_kernel - 1 ; ++l) {
                     temp += src[k * ws.w_src + l] * kernel[(i + ws.h_kernel - 1 - k) * ws.w_kernel + (j + ws.w_kernel - 1 - l)];
                  }
               }
               ws.dst[i * ws.w_dst + j] = temp;
            }
         }
         break;
   }
}
*/

#endif
