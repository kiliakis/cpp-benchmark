#include <iostream>
#include <chrono>
#include <vector>
#include <cmath>
#include "sin.h"

const unsigned N = 20000000;

const double B = 4.0/M_PI;
const double C = -4.0/(M_PI*M_PI);
const double P = 0.225;

double inline fastSin( double x ){
    x = fmod(x + M_PI, M_PI * 2) - M_PI; // restrict x so that -M_PI < x < M_PI

    const double y = B * x + C * x * std::abs(x);


    return P * (y * std::abs(y) - y) + y; 
}

const double PI = M_PI;
const double HALF_PI = 0.5 * PI;
const double TWO_PI = 2.0 * PI;
const double TWO_PI_INV = 1.0 / TWO_PI;
 
const double a0 = 1.0;
const double a2 = 2.0 / PI - 12.0 / (PI * PI);
const double a3 = 16.0 / (PI * PI * PI) - 4.0 / (PI * PI);

inline double  Hill(double x)
{
  const double xx = x * x;
  const double xxx = xx * x;
 
  return a0 + a2 * xx + a3 * xxx;
}
 
inline double fastestSin(double x)
{
  // wrap x within [0, TWO_PI)
  const double a = x * TWO_PI_INV;
  x -= static_cast<int>(a) * TWO_PI;
  if (x < 0.0)
    x += TWO_PI;
 
  // 4 pieces of hills
  if (x < HALF_PI)
    return Hill(HALF_PI - x);
  else if (x < PI)
    return Hill(x - HALF_PI);
  else if (x < 3.0 * HALF_PI)
    return -Hill(3.0 * HALF_PI - x);
  else
    return -Hill(x - 3.0 * HALF_PI);
}

int main()
{
   
   std::vector<double> v(N);
   auto start = std::chrono::system_clock::now();
   auto end = std::chrono::system_clock::now();
   auto sum = 0.0;
   std::chrono::duration<double> elapsed = end - start;
   /*
   auto start = std::chrono::system_clock::now();
   for(unsigned i = 0; i < v.size(); i++)
      v[i] = sin(i);

   auto end = std::chrono::system_clock::now();
   std::chrono::duration<double> elapsed = end - start;
   
   auto sum = 0.0;
   for(const auto &n: v)
      sum +=n;
   std::cout << "---- sin ----\n";
   std::cout << "sum : " << sum << std::endl;
   std::cout << "elapsed time : " << elapsed.count() << std::endl;
   */

   //v.clear();
   start = std::chrono::system_clock::now();
   for(unsigned i = 0; i < v.size(); i++)
      v[i] = std::sin(i);

   end = std::chrono::system_clock::now();
   elapsed = end - start;
   
   sum = 0.0;
   for(const auto &n: v)
      sum +=n;
   std::cout << "---- std::sin ----\n";
   std::cout << "sum : " << sum << std::endl;
   std::cout << "elapsed time : " << elapsed.count() << std::endl;

   start = std::chrono::system_clock::now();
   
   for(unsigned i = 0; i < v.size(); i++)
      v[i] = vdt::fast_sin(i);

   end = std::chrono::system_clock::now();
   elapsed = end - start;
   
   sum = 0.0;
   for(const auto &n: v)
      sum +=n;
   std::cout << "---- vdt::fast_sin ----\n";
   std::cout << "sum : " << sum << std::endl;
   std::cout << "elapsed time : " << elapsed.count() << std::endl;

   start = std::chrono::system_clock::now();
   
   for(unsigned i = 0; i < v.size(); i++)
      v[i] = fastSin(i);

   end = std::chrono::system_clock::now();
   elapsed = end - start;
   
   sum = 0.0;
   for(const auto &n: v)
      sum +=n;
   std::cout << "---- fastSin ----\n";
   std::cout << "sum : " << sum << std::endl;
   std::cout << "elapsed time : " << elapsed.count() << std::endl;

   start = std::chrono::system_clock::now();
   
   for(unsigned i = 0; i < v.size(); i++)
      v[i] = fastestSin(i);

   end = std::chrono::system_clock::now();
   elapsed = end - start;
   
   sum = 0.0;
   for(const auto &n: v)
      sum +=n;
   std::cout << "---- fastestSin ----\n";
   std::cout << "sum : " << sum << std::endl;
   std::cout << "elapsed time : " << elapsed.count() << std::endl;



   return 0;
}
