/*
 * math_functions.h
 *
 *  Created on: Mar 21, 2016
 *      Author: kiliakis
 */

#ifndef INCLUDES_MATH_FUNCTIONS_H_
#define INCLUDES_MATH_FUNCTIONS_H_

#include <cmath>
#include "sin.h"
#include <omp.h>
#include  <cassert>
#include "utilities.h"
#include "configuration.h"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_errno.h>
#include <algorithm>


#ifdef USE_FFTW
#include <fftw3.h>
#else
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_fft_complex.h>
#endif

namespace mymath {

   // Wrapper function for vdt::fast_sin
   static inline ftype fast_sin(ftype x)
   {
      return vdt::fast_sin(x);
   }


   static inline ftype fast_cos(ftype x)
   {
      return vdt::fast_sin(x + M_PI_2);
   }

   // linear convolution function
   // @a: first vector
   // @b: second vector
   // @return convolution of a and b
   template<typename T>
   static inline std::vector<T> convolution(const std::vector<T> &a,
         const std::vector<T> &b)
   {
      std::vector<T> res;
      res.resize(a.size() + b.size() - 1);
      for (unsigned int i = 0; i < res.size(); i++) {
         unsigned int i1 = i;
         T temp = T();
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


   static inline void real_to_complex(const std::vector<ftype> &in,
                                      std::vector<complex_t> &out)
   {
      assert(out.empty());
      out.reserve(in.size());
      for (const auto &real : in)
         out.push_back(complex_t(real, 0));
   }


   static inline void pack_to_complex(const std::vector<ftype> &in,
                                      std::vector<complex_t> &out)
   {
      assert(out.empty());
      assert(in.size() % 2 == 0);
      out.reserve(in.size() / 2);
      for (unsigned int i = 0; i < in.size(); i += 2)
         out.push_back(complex_t(in[i], in[i + 1]));
   }


   static inline void complex_to_real(const std::vector<complex_t> &in,
                                      std::vector<ftype> &out)
   {
      assert(out.empty());
      out.reserve(in.size());
      for (const auto &z : in)
         out.push_back(z.real());
   }

   static inline void unpack_complex(const std::vector<complex_t> &in,
                                     std::vector<ftype> &out)
   {
      assert(out.empty());
      out.reserve(2 * in.size());
      for (const auto &z : in) {
         out.push_back(z.real());
         out.push_back(z.imag());
      }
   }


#ifdef USE_FFTW
   static inline fftw_plan init_fft(const int n,
                                    complex_t *in,
                                    complex_t *out,
                                    const int sign = FFTW_FORWARD,
                                    const unsigned flag = FFTW_ESTIMATE,
                                    const int threads = 1)
   {
      if (threads > 1) {
         fftw_init_threads();
         fftw_plan_with_nthreads(threads);
      }
      fftw_complex *a, *b;
      a = reinterpret_cast<fftw_complex *>(in);
      b = reinterpret_cast<fftw_complex *>(out);
      return fftw_plan_dft_1d(n, a, b, sign, flag);
   }



   static inline fftw_plan init_rfft(const int n,
                                     ftype *in,
                                     complex_t *out,
                                     const unsigned flag = FFTW_ESTIMATE,
                                     const int threads = 1)

   {
      if (threads > 1) {
         fftw_init_threads();
         fftw_plan_with_nthreads(threads);
      }
      fftw_complex *b;
      b = reinterpret_cast<fftw_complex *>(out);
      return fftw_plan_dft_r2c_1d(n, in, b, flag);
   }

   static inline fftw_plan init_irfft(const int n,
                                      complex_t *in,
                                      ftype *out,
                                      const unsigned flag = FFTW_ESTIMATE,
                                      const int threads = 1)
   {
      if (threads > 1) {
         fftw_init_threads();
         fftw_plan_with_nthreads(threads);
      }
      fftw_complex *b;
      b = reinterpret_cast<fftw_complex *>(in);
      return fftw_plan_dft_c2r_1d(n, b, out, flag);

   }

   static inline void run_fft(const fftw_plan &p)
   {
      fftw_execute(p);
   }

   static inline void destroy_fft(fftw_plan &p)
   {
      fftw_destroy_plan(p);
   }


#endif



   // Parameters are like python's numpy.fft.rfft
   // @in:  input data
   // @n:   number of points to use. If n < in.size() then the input is cropped
   //       if n > in.size() then input is padded with zeros
   // @out: the transformed array

   static inline void rfft(f_vector_t &in,
                           complex_vector_t &out,
                           uint n = 0,
                           const uint threads = 1)
   {
      if (n == 0)
         n = in.size();
      else
         in.resize(n);

#ifdef USE_FFTW
      out.resize(n / 2 + 1);
      //out.resize(n);
      auto p = mymath::init_rfft(n, in.data(), out.data(),
                                 FFTW_ESTIMATE, threads);
      mymath::run_fft(p);
      mymath::destroy_fft(p);

#else
      std::cerr << "Use of gsl ffts is depricated\n";


      gsl_fft_real_wavetable *real;
      gsl_fft_real_workspace *work;

      work = gsl_fft_real_workspace_alloc(n);
      real = gsl_fft_real_wavetable_alloc(n);

      gsl_fft_real_transform(in.data(), 1, n, real, work);

      out.clear();
      out.reserve(in.size() / 2);
      // first element is real => imag is zero
      in.insert(in.begin() + 1, 0.0);

      // if n is even => last element is real
      if (n % 2 == 0)
         in.push_back(0.0);

      //pack_to_complex(v, out);
      pack_to_complex(in, out);

      gsl_fft_real_wavetable_free(real);
      gsl_fft_real_workspace_free(work);

#endif
   }


   // Parameters are like python's numpy.fft.fft
   // @in:  input data
   // @n:   number of points to use. If n < in.size() then the input is cropped
   //       if n > in.size() then input is padded with zeros
   // @out: the transformed array

   static inline void fft(complex_vector_t &in,
                          complex_vector_t &out,
                          uint n = 0,
                          const uint threads = 1)
   {
      if (n == 0)
         n = in.size();

#ifdef USE_FFTW
      out.resize(n);

      auto p = mymath::init_fft(n, in.data(), out.data(), FFTW_FORWARD,
                                FFTW_ESTIMATE, threads);
      mymath::run_fft(p);
      mymath::destroy_fft(p);

#else
      std::cerr << "Use of gsl ffts is depricated\n";

      std::vector<ftype> v;
      //v.resize(2 * n, 0);
      unpack_complex(in, v);

      gsl_fft_complex_wavetable *wave;
      gsl_fft_complex_workspace *work;

      wave = gsl_fft_complex_wavetable_alloc(n);
      work = gsl_fft_complex_workspace_alloc(n);

      gsl_fft_complex_forward(v.data(), 1, n, wave, work);

      //printf("ok inside\n");

      out.clear();

      pack_to_complex(v, out);

      out.resize(n, 0);

      gsl_fft_complex_wavetable_free(wave);
      //printf("ok here7\n");

      gsl_fft_complex_workspace_free(work);
      //printf("ok here8\n");

#endif

   }


// Parameters are like python's numpy.fft.ifft
// @in:  input data
// @n:   number of points to use. If n < in.size() then the input is cropped
//       if n > in.size() then input is padded with zeros
// @out: the inverse Fourier transform of input data

   static inline void ifft(complex_vector_t &in,
                           complex_vector_t &out,
                           uint n = 0,
                           const uint threads = 1)
   {
      if (n == 0)
         n = in.size();

#ifdef USE_FFTW
      out.resize(n);
      auto p = mymath::init_fft(n, in.data(), out.data(), FFTW_BACKWARD,
                                FFTW_ESTIMATE, threads);
      mymath::run_fft(p);
      std::transform(out.begin(), out.end(), out.begin(),
                     std::bind2nd(std::divides<complex_t>(), n));
      mymath::destroy_fft(p);

#else
      std::cerr << "Use of gsl ffts is depricated\n";

      std::vector<ftype> v;
      //v.resize(2 * n, 0);

      unpack_complex(in, v);

      gsl_fft_complex_wavetable *wave;
      gsl_fft_complex_workspace *work;

      work = gsl_fft_complex_workspace_alloc(n);
      wave = gsl_fft_complex_wavetable_alloc(n);

      gsl_fft_complex_inverse(v.data(), 1, n, wave, work);


      out.clear();

      pack_to_complex(v, out);

      out.resize(n, 0);

      gsl_fft_complex_wavetable_free(wave);
      gsl_fft_complex_workspace_free(work);
#endif
   }


// Inverse of rfft
// @in: input vector which must be the result of a rfft
// @out: irfft of input, always real
// Missing n: size of output
// TODO fix this one!!
   static inline void irfft(complex_vector_t in,
                            f_vector_t &out,
                            uint n = 0,
                            const uint threads = 1)
   {
      n = (n == 0) ? 2 * (in.size() - 1) : n;

#ifdef USE_FFTW
      out.resize(n);

      auto p = mymath::init_irfft(n, in.data(), out.data(),
                                  FFTW_ESTIMATE, threads);
      mymath::run_fft(p);
      std::transform(out.begin(), out.end(), out.begin(),
                     std::bind2nd(std::divides<ftype>(), n));
      mymath::destroy_fft(p);


#else
      std::cerr << "Use of gsl ffts is depricated\n";

      assert(in.size() > 1);

      uint last = in.size() - 2;

      if (n == 2 * in.size() - 1) {
         last = in.size() - 1;
      } else if (n == 2 * (in.size() - 1)) {
         ;
      } else {
         std::cerr << "[mymath::ifft] Size not supported!\n"
                   << "[mymath::ifft] default size: " << n
                   << " will be used\n";
      }

      for (uint i = last; i > 0; --i) {
         in.push_back(std::conj(in[i]));
      }

      complex_vector_t temp;
      mymath::ifft(in, n, temp);

      mymath::complex_to_real(temp, out);
#endif

   }


// Same as python's numpy.fft.rfftfreq
// @ n: window length
// @ d (optional) : Sample spacing
// @return: A vector of length (n div 2) + 1 of the sample frequencies
   static inline f_vector_t rfftfreq(const uint n, const ftype d = 1.0)
   {
      f_vector_t v(n / 2 + 1);
      const ftype factor = 1.0 / (d * n);
      #pragma omp parallel for
      for (uint i = 0; i < v.size(); ++i) {
         v[i] = i * factor;
      }
      return std::move(v);
   }


// Parameters are like python's np.interp
// @x: x-coordinates of the interpolated values
// @xp: The x-coords of the data points
// @fp: the y-coords of the data points
// @y: the interpolated values, same shape as x
// @left: value to return for x < xp[0]
// @right: value to return for x > xp[last]
   static inline void lin_interp(const std::vector<ftype> &x, const std::vector<ftype> &xp,
                                 const std::vector<ftype> &fp, std::vector<ftype> &y,
                                 const ftype left = 0, const ftype right = 0)
   {
      assert(y.empty());

      gsl_interp *interp =
         gsl_interp_alloc(gsl_interp_linear, xp.size());

      gsl_interp_init(interp, &xp[0], &fp[0], xp.size());

      gsl_interp_accel *acc = gsl_interp_accel_alloc();

      for (uint i = 0; i < x.size(); ++i) {
         double val;
         if (x[i] < interp->xmin) {
            val = left;
         } else if (x[i] > interp->xmax) {
            val = right;
         } else {
            val = gsl_interp_eval(interp, &xp[0],
                                  &fp[0], x[i],
                                  acc);
         }
         y.push_back(val);
      }

      gsl_interp_free(interp);
      gsl_interp_accel_free(acc);



   }


// Function to implement integration of f(x) over the interval
// [a,b] using the trapezoid rule with nsub subdivisions.
   template<typename T>
   static inline ftype *cum_trapezoid(const T *f, const T deltaX, const int nsub)
   {
      // initialize the partial sum to be f(a)+f(b) and
      // deltaX to be the step size using nsub subdivisions
      ftype *psum = new ftype[nsub];
      psum[0] = 0;

      // increment the partial sum
      //#pragma omp parallel for
      for (int i = 1; i < nsub; ++i)
         psum[i] = psum[i - 1] + (f[i] + f[i - 1]) * (deltaX / 2.0);

      return psum;

   }

   template<typename T>
   static inline ftype trapezoid(const T *f, const T *deltaX, const int nsub)
   {
      // initialize the partial sum to be f(a)+f(b) and
      // deltaX to be the step size using nsub subdivisions

      ftype psum = 0.0;
      // increment the partial sum
      #pragma omp parallel for reduction(+ : psum)
      for (int index = 1; index < nsub; ++index) {
         psum += (f[index] + f[index - 1])
                 * (deltaX[index] - deltaX[index - 1]);
      }

      // return approximation
      return psum / 2;

   }

   template<typename T>
   static inline ftype trapezoid(const T *f, const T deltaX, const int nsub)
   {
      // initialize the partial sum to be f(a)+f(b) and
      // deltaX to be the step size using nsub subdivisions
      ftype psum = f[0] + f[nsub - 1]; //f(a)+f(b);
      //ftype deltaX = (b-a)/nsub;

      // increment the partial sum
      #pragma omp parallel for reduction(+ : psum)
      for (int index = 1; index < nsub - 1; ++index) {
         psum += 2 * f[index];
      }

      // multiply the sum by the constant deltaX/2.0
      psum = (deltaX / 2) * psum;

      // return approximation
      return psum;

   }
   template<typename T>
   static inline int min(T *a, int size, int step = 1)
   {
      int p = 0;
      T min = a[0];
      //#pragma omp parallel for  shared(p) reduction(min : min)
      for (int i = 1; i < size; i += step) {
         if (a[i] < min) {
            min = a[i];
            p = i;
         }
      }
      return p;

   }

   template<typename T>
   static inline int max(T *a, int size, int step = 1)
   {
      int p = 0;
      T max = a[0];
      //#pragma omp parallel for shared(p) reduction(max : max)
      for (int i = 1; i < size; i += step) {
         if (a[i] > max) {
            max = a[i];
            p = i;
         }
      }
      return p;

   }

   static inline void linspace(ftype *a, const ftype start, const ftype end,
                               const int n, const int keep_from = 0)
   {
      const ftype step = (end - start) / (n - 1);
      //ftype value = start;
      #pragma omp parallel for
      for (int i = 0; i < n; ++i) {
         if (i >= keep_from)
            a[i - keep_from] = start + i * step;
         //value += step;
      }
   }

   template<typename T>
   static inline ftype mean(const T data[], const int n)
   {
      ftype m = 0.0;
      #pragma omp parallel for reduction(+ : m)
      for (int i = 0; i < n; ++i) {
         m += data[i];
      }
      return m / n;
   }

   template<typename T>
   static inline ftype standard_deviation(const T data[], const int n,
                                          const ftype mean)
   {
      ftype sum_deviation = 0.0;
      #pragma omp parallel for reduction(+ : sum_deviation)
      for (int i = 0; i < n; ++i)
         sum_deviation += (data[i] - mean) * (data[i] - mean);
      return std::sqrt(sum_deviation / n);
   }

   template<typename T>
   static inline ftype standard_deviation(const T data[], const int n)
   {
      const ftype mean = mymath::mean(data, n);
      ftype sum_deviation = 0.0;
      #pragma omp parallel for reduction(+ : sum_deviation)
      for (int i = 0; i < n; ++i)
         sum_deviation += (data[i] - mean) * (data[i] - mean);
      return std::sqrt(sum_deviation / n);
   }

}

#endif /* INCLUDES_MATH_FUNCTIONS_H_ */
