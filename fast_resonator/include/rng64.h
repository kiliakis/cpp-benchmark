/**********************************************************************

  2D3V PIC-MCC CODE '2D Arc-PIC'
  (C) H. Timko, 2010-2011

  rng64.h
  Header file for rng64.cpp

***********************************************************************/
#ifndef RNG64_H_
#define RNG64_H_

extern   unsigned long long  seed;

void  initrand(unsigned long long *seed, int, int, int, int);

double  Random(unsigned long long *seed);

void   newrands(unsigned long long *seed, int *s);

#endif /* RNG64_H_ */