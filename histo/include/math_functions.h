/*
 * math_functions.h
 *
 *  Created on: Mar 21, 2016
 *      Author: kiliakis
 */

#ifndef INCLUDES_MATH_FUNCTIONS_H_
#define INCLUDES_MATH_FUNCTIONS_H_

#include <algorithm>
#include <utilities.h>
#include <cmath>
#include <omp.h>

namespace mymath {


    template<typename T>
    static inline uint min(T *a,
                           uint size,
                           uint step = 1)
    {
        uint p = 0;
        T min = a[0];
        //#pragma omp parallel for  shared(p) reduction(min : min)
        for (uint i = 1; i < size; i += step) {
            if (a[i] < min) {
                min = a[i];
                p = i;
            }
        }
        return p;

    }

    template<typename T>
    static inline uint max(T *a,
                           uint size,
                           uint step = 1)
    {
        uint p = 0;
        T max = a[0];
        //#pragma omp parallel for shared(p) reduction(max : max)
        for (uint i = 1; i < size; i += step) {
            if (a[i] > max) {
                max = a[i];
                p = i;
            }
        }
        return p;

    }


    static inline void linspace(ftype *a,
                                const ftype start,
                                const ftype end,
                                const uint n,
                                const uint keep_from = 0)
    {
        const ftype step = (end - start) / (n - 1);
        //ftype value = start;
        //#pragma omp parallel for
        for (uint i = 0; i < n; ++i) {
            if (i >= keep_from)
                a[i - keep_from] = start + i * step;
            //value += step;
        }
    }


    template<typename T>
    static inline std::vector<T> arange(T start,
                                        T stop,
                                        T step = 1)
    {
        std::vector<T> values;
        for (T value = start; value < stop; value += step)
            values.push_back(value);
        return values;
    }


    template<typename T>
    static inline ftype mean(const T data[],
                             const int n)
    {
        ftype m = 0.0;
        #pragma omp parallel for reduction(+ : m)
        for (int i = 0; i < n; ++i) {
            m += data[i];
        }
        return m / n;
    }


    template<typename T>
    static inline ftype standard_deviation(const T data[],
                                           const int n,
                                           const ftype mean)
    {
        ftype sum_deviation = 0.0;
        #pragma omp parallel for reduction(+ : sum_deviation)
        for (int i = 0; i < n; ++i)
            sum_deviation += (data[i] - mean) * (data[i] - mean);
        return std::sqrt(sum_deviation / n);
    }


    template<typename T>
    static inline ftype standard_deviation(const T data[],
                                           const int n)
    {
        const ftype mean = mymath::mean(data, n);
        ftype sum_deviation = 0.0;
        #pragma omp parallel for reduction(+ : sum_deviation)
        for (int i = 0; i < n; ++i)
            sum_deviation += (data[i] - mean) * (data[i] - mean);
        return std::sqrt(sum_deviation / n);
    }


    template <typename T>
    int sign(T val)
    {
        return (T(0) < val) - (val < T(0));
    }

}
#endif /* INCLUDES_MATH_FUNCTIONS_H_ */
