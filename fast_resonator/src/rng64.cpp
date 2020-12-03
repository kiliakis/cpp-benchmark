/**********************************************************************

  2D3V PIC-MCC CODE '2D Arc-PIC'
  (C) H. Timko, 2010-2011

  rng64.cpp
  Random number genertor after Friedrich Hertweck

***********************************************************************/

#include <stdlib.h>
#include <stdio.h>



struct  random64 {
    unsigned long long   seed, one;
    int           sequ, mask, shift;
};


extern struct random64  r = { 0, 0x3ff0000000000000LLU, -1, -1, -1 };



void  RNGERR(int fehler)
{
    /*   printf ("%4d> Error found by RNG:  ", _my_pe() ); */

    switch (fehler) {
        case  1:
            printf("second initialization not permitted\n");
            break;

        case  2:
            printf("value of ctrl is not 1 or 2\n");
            break;

        case  3:
            printf("value of  kr  must be at least 24\n");
            break;

        case  4:
            printf("value of  ks  exceeds 24\n");
            break;

        case  5:
            printf("sum  kr + kp + ks  exceeds 64\n");
            break;

        case  6:
            printf("dimension of nCube exceeds  kp\n");
            break;

        case  7:
            printf("number of permitted sequences exhausted\n");
            break;
    }

    /*   printf ("%4d> *****  Program terminated  *****\n", _my_pe() ); */
    exit(1000 + fehler);
}




void  initrand(unsigned long long      *seed, int ctrl, int kr, int kp, int ks)
{

    unsigned long long   t, t1;
    int                  n, pid, dim;
    // int                  proc, host;

    if (r.shift > 0)  RNGERR(1);
    if ((ctrl < 1) || (ctrl > 2))  RNGERR(2);
    if (kr < 24)  RNGERR(3);
    if (ks > 24)  RNGERR(4);
    if ((kr + kp + ks) > 64)  RNGERR(5);

    // whoami (&pid, &proc, &host, &dim);
    pid = dim = 0;  // only 1 PE ;-)
    if (dim > kp)  RNGERR(6);

    r.one = 0x3ff0000000000000LLU;
    r.sequ = 0;
    r.mask = (1 << ks) - 1;

    if (ctrl == 1) {
        r.shift = kr + kp;
        n = kr;
    } else {
        r.shift = kr;
        n = kr + ks;
    }

    t = pid;  t = t << n;
    t1 = 1;  t1 = (t1 << (64 - kr)) - 1;
    r.seed = (*seed & t1) | t;

    *seed = r.seed;
}




void   newrands(unsigned long long *seed, int *s)
{

    unsigned long long   t;

    r.sequ = (r.sequ + 1) & r.mask;
    if (r.sequ == 0)  RNGERR(7);

    t = r.sequ;  t = t << r.shift;
    *seed = r.seed | t;
    *s = r.sequ;
}




double  Random(unsigned long long      *seed)
{

    unsigned long long    a = 6364136223846793005LLU;
    unsigned long long    x;

    *seed = (*seed) * a + 1;

    x = (*seed >> 12) | r.one;
    return (*(double *)&x - 1.0);
}




/*=============================================================================
*       Parallelized Linear Congruential Random Number Generator
*                          for the nCube 2
*
*  Copyright:  Friedrich Hertweck
*              Max Planck Institute of Plasma Physics
*               D-8046 Garching, Germany
*==============================================================================
*
*  The package consists of three functions:
*
*------------------------------------------------------------------------------
*
*  FORTRAN:   call INITRAND (seed, ctrl, kr, kp, ks)
*  C:         initrand (&seed, ctrl, kr, kp, ks)
*
*  inputs:
*       seed = global initial seed (should be less than 2**kr)
*       ctrl = 1:  pid first, sequ second
*              2:  sequ first, pid second
*       kr   = number of bits of straight run of rundom numbers
*       kp   = number of bits for pid (must be >= dim of nCube)
*       ks   = number of bits for sequence numbers
*
*  result:
*       seed =      current seed , as described below
*
*  errors (program will stop with one of the messages):
*       1   second initialization not permitted
*       2   input value of ctrl not 1 or 2
*       3   input value of kr is < 24
*       4   input value of ks is > 24
*       5   the sum of  kr + kp + ks exceeds 64
*       6   dimension of nCube exceeds kp
*
*  Initializes the generator. May only be called once. The sequences to
*  be produced have a maximum length of 2**kr; the first sequence has
*  sequence id = 0. A check is made that  kr + kp + ks <= 64, so that each
*  sequence on each processor is different. If the initial seed is >= 2**kr,
*  it is truncated without warning to kr bits; this truncated value, denoted
*  by seed0, will be used further on as the global seed for constructing
*  initial seeds for sequences.
*
*  When ctrl = 1, the initial seeds produced by INITRAND and NEWRANDS are
*
*      seed = sequ * 2**(kr + kp) + pid * 2**kr + seed0,
*
*  when ctrl = 2, the seeds are
*
*      seed = pid * 2**(kr + ks) + sequ * 2**kr + seed0.
*
*  NOTE:  When the initialization returns an error, the function
*         RANDOM will always return zero as "random numbers"!
*
*------------------------------------------------------------------------------
*
*  FORTRAN:    call  NEWRANDS (seed, sequ)
*  C:          newrands (&seed, &sequ)
*
*  inputs:
*       none
*
*  results:
*       seed =    seed for a new sequence
*       sequ =    id of the new sequence
*
*  errors (program will stop with the message):
*       7   number of sequences exhausted
*
*  Generate the initial seed for the next sequence, the id number of
*  which is returned in sequ.
*
*------------------------------------------------------------------------------
*
*  FORTRAN:    r = RANDOM (seed)
*  C:          r = random (&seed)
*
*  inputs:
*       seed =    current seed
*
*  results:
*       seed =    next seed
*       r    =    next random number in the range 0.0 <= r < 1.0
*
*============================================================================*/
