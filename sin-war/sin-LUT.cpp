#include <iostream>
#include <chrono>
#include <vector>
#include <cmath>
#include "sin.h"
#include "fastSinCos.h"

const unsigned N = 20000000;


int main()
{
   std::vector<double> v(N);
   std::vector<double> res(N);
   for (unsigned i = 0; i < v.size(); i++) {
      float r2 = static_cast <float>(rand()) / (static_cast <float>(RAND_MAX / 100.0));
      v[i] = r2;
   }

   auto start = std::chrono::system_clock::now();
   auto end = std::chrono::system_clock::now();
   auto sum = 0.0;
   std::chrono::duration<double> elapsed = end - start;



   // ------------- std::sin()
   start = std::chrono::system_clock::now();
   for (unsigned i = 0; i < v.size(); i++)
      res[i] = std::sin(v[i]);

   end = std::chrono::system_clock::now();
   elapsed = end - start;

   sum = 0.0;
   for (const auto &n : res)
      sum += n;

   std::cout << "---- std::sin ----\n";
   std::cout << "sum : " << sum << std::endl;
   std::cout << "elapsed time : " << elapsed.count() << std::endl;
   
   // -------------- fastSin

   start = std::chrono::system_clock::now();

   for (unsigned i = 0; i < v.size(); i++)
      res[i] = fastSin(v[i]);

   end = std::chrono::system_clock::now();
   elapsed = end - start;

   sum = 0.0;
   for (const auto &n : res)
      sum += n;
   std::cout << "---- fastSin ----\n";
   std::cout << "sum : " << sum << std::endl;
   std::cout << "elapsed time : " << elapsed.count() << std::endl;


   // -------------- vdt::fast_sin
   start = std::chrono::system_clock::now();

   for (unsigned i = 0; i < v.size(); i++)
      res[i] = vdt::fast_sin(v[i]);

   end = std::chrono::system_clock::now();
   elapsed = end - start;

   sum = 0.0;
   for (const auto &n : res)
      sum += n;
   std::cout << "---- vdt::fast_sin ----\n";
   std::cout << "sum : " << sum << std::endl;
   std::cout << "elapsed time : " << elapsed.count() << std::endl;


   return 0;
}
