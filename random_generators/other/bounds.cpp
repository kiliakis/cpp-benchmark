/**********************************************************************

  2D3V PIC-MCC CODE '2D Arc-PIC'
  (C) H. Timko, 2010-2011

  bounds.cpp
  2D boundary conditions

***********************************************************************/

#include <math.h>
#include  <stdio.h> 

#include  "pic.h"
#include  "dim.h"
#include  "mydef.h" 

#define   XTRN extern
#include  "rng64.h"
#define   RAND  Random(&seed)

#include  "var.h"
   


/***********************
 Remove particles
***********************/

void  remove_particles( Particle pa[], int *np, double  zmin, double  zmax, double rmax, 
			int *n_wall1, int *n_wall2, int *n_inf )
{
  int n, n_lost = 0;
  double delta;

  for ( n=0; n<*np; n++ )
    {
     
      if ( pa[n].p.r >= rmax ) 	
	{
	  /*
	  delta = pa[n].p.r - rmax;
	  pa[n].p.r -= 2.*delta;
	  pa[n].p.vr *= -1.;
	  */
	  *n_inf +=1;
	  n_lost++;
	  continue; 
	}
      else if ( pa[n].p.z >= zmax ) 		
	{
	  *n_wall2 += 1;
	  n_lost++;
	  continue; 
	}
      else if ( pa[n].p.z < zmin ) 	
	{
	  *n_wall1 +=1;
	  n_lost++;
	  continue; 
	}


      pa[n-n_lost].p = pa[n].p;
    }

  *np -= n_lost;

}   


void  remove_arc_e( Particle pa[], int *np, double  zmin, double  zmax, double rmax, 
		    int *n_wall1, int *n_wall2, int *n_inf, int *I_wall1, int *I_wall2, int *I_C, int *I_A )
{
  int jr, n, n_lost = 0;
  double delta;

  for ( n=0; n<*np; n++ )
    {
     
      if ( pa[n].p.r >= rmax ) 	
	{
	  *n_inf +=1;
	  n_lost++;
	  continue; 
	}
      else if ( pa[n].p.z >= zmax ) 
	{
	  jr = int (pa[n].p.r);
	  *n_wall2 += 1;
	  *I_wall2 += 2*jr+1;
	  *I_A += 1;
	  n_lost++;
	  continue; 
	}
      else if ( pa[n].p.z < zmin ) 	
	{
	  jr = int (pa[n].p.r);
	  *n_wall1 +=1;
	  *I_wall1 += 2*jr+1;
	  *I_C -= 1;
	  n_lost++;
	  continue; 
	}


      pa[n-n_lost].p = pa[n].p;
    }

  *np -= n_lost;

}   


void  remove_arc_n( Particle pa[], int *np, double zmin, double zmax, double rmax, 
		    int *n_wall1, int *n_wall2, int *n_inf,   
		    int *ind1, Sput s1[], int *ind2, Sput s2[], double cs_n, double Te )
{

  int n;
  int n_lost = 0; 

  double nrg_w1, eps_1, sn_1, nrg_w2, eps_2, sn_2;

  // Yield: fractional and integer part
  double p; 
  int Y;



  for ( n=0; n<*np; n++ )
    {

      if ( pa[n].p.r >= rmax )
	{
	  *n_inf +=1;
	  n_lost++;
	  continue; 
	}


      else if ( pa[n].p.z >= zmax ) 		
	{	
 	
	  *n_wall2  += 1; 

	  nrg_w2 = SQU(pa[n].p.vz) + SQU(pa[n].p.vr) + SQU(pa[n].p.vt); 
	  nrg_w2 *= Te/(2*SQU(cs_n));
       
	  if ( nrg_w2 > 23.383 ) // in eV
	    {
	      eps_2 = (4.45394e-06)*nrg_w2; 
	      sn_2 = 8205*3.441*sqrt(eps_2)*log(eps_2 + 2.718)/(1 + 6.355*sqrt(eps_2) + eps_2*(6.882*sqrt(eps_2) - 1.708));
	      p = 0.042*(0.2525/3.49)*(sn_2/(1 + 0.010124*0.1573*pow(eps_2,0.3)))*pow(1 - sqrt(23.383/nrg_w2),2.5);

	      // Fractional part of p handled probabilistically
	      Y = int ( p );
	      p -= Y;
	      if ( RAND <= p )
		Y++;
		
	      // NEW: register r-coordinates of bombarding particles
	      if ( Y > 0 )
		{
		  s2[*ind2].Y = Y;
		  s2[*ind2].r = pa[n].p.r;
		  (*ind2)++;
		}

	    }

	  n_lost++; 
	  continue; 
	}


      else if ( pa[n].p.z < zmin ) 	
	{
	  
	  *n_wall1  += 1;
 
	  nrg_w1 = SQU(pa[n].p.vz) + SQU(pa[n].p.vr) + SQU(pa[n].p.vt);
	  nrg_w1 *= Te/(2*SQU(cs_n));
	  
	  if ( nrg_w1 > 23.383 )
	    {
	      // Reduced energy
	      eps_1 = (4.45394e-06)*nrg_w1;
	      
	      // Nuclear stopping cross section
	      sn_1 = 8205*3.441*sqrt(eps_1)*log(eps_1 + 2.718)/(1 + 6.355*sqrt(eps_1) + eps_1*(6.882*sqrt(eps_1) - 1.708));
	      
	      // Sputtering yield 
	      p = 0.042*(0.2525/3.49)*(sn_1/(1 + 0.010124*0.1573*pow(eps_1,0.3)))*pow(1 - sqrt(23.383/nrg_w1),2.5);

	      // Fractional part of p handled probabilistically
	      Y = int ( p );
	      p -= Y;
	      if ( RAND <= p )
		Y++;
		
	      // Register r-coordinates of bombarding particles
	      if ( Y > 0 )
		{
		  s1[*ind1].Y = Y;
		  s1[*ind1].r = pa[n].p.r;
		  (*ind1)++;
		}

	    }
	  
	  n_lost++; 
	  continue; 
	}
      pa[n-n_lost].p = pa[n].p; 
    }

  *np -= n_lost; 

}


void  remove_arc_i( Particle pa[], int *np, double zmin, double zmax, double rmax, 
		    int *n_wall1, int *n_wall2, int *n_inf, int *I_wall1, int *I_wall2, int *I_C, int *I_A,
		    double cs_i, double Te, double n0, double NDebye, 
		    Sput s1_temp[], int *ind1, Sput s1[], int *ind2, Sput s2[], 
		    int *inde, Sput se[], double SEY, int nr, double dr, int ion_step )
{

  int n, n_lost = 0;

  // check threshold of sputtering yield, only for ions, only at the cathode
  // 2D: flux will now depend on the area!! taken into account further below
  double threshold = 1.e7*(NDebye*ion_step*Omega_pe)/(6.7193e-12*n0*sqrt(Te));  // threshold = 1e7 A/cm^2
  double current [nr];

  double nrg_w1, eps_1, sn_1, nrg_w2, eps_2, sn_2;
  
  // Yield: fractional and integer part
  double p; 
  int Y, Y_sum = 0;
  int i, j, jr;
  int k = 0;
  
  // SEY integer, fractional part, arrays
  int SEY_i = int ( SEY );
  double SEY_f = SEY - SEY_i;

  double r, sigma;

  // Initialise array
  for ( i=0; i<nr; i++ )
    current[i] = 0.;


  for ( n=0; n<*np; n++ )
    {

      if ( pa[n].p.r >= rmax )
	{
	  *n_inf +=1;
	  n_lost++;
	  continue; 
	}
      

      else if ( pa[n].p.z >= zmax ) 		
	{	
	  jr = int (pa[n].p.r);
	  *n_wall2 += 1;
	  *I_wall2 += 2*jr+1;
	  *I_A -= 1;
	  
	  nrg_w2 = SQU(pa[n].p.vz) + SQU(pa[n].p.vr) + SQU(pa[n].p.vt);
	  nrg_w2 *= Te/(2*SQU(cs_i));
	  
	  if ( nrg_w2 > 23.383 ) // in eV
	    {
	      eps_2 = (4.45394e-06)*nrg_w2; 
	      sn_2 = 8205*3.441*sqrt(eps_2)*log(eps_2 + 2.718)/(1 + 6.355*sqrt(eps_2) + eps_2*(6.882*sqrt(eps_2) - 1.708));
	      p = 0.042*(0.2525/3.49)*(sn_2/(1 + 0.010124*0.1573*pow(eps_2,0.3)))*pow(1 - sqrt(23.383/nrg_w2),2.5);
	      
	      // Fractional part of p handled probabilistically
	      Y = int ( p );
	      p -= Y;
	      if ( RAND <= p )
		Y++;
		
	      // Register r-coordinates of bombarding particles
	      if ( Y > 0 )
		{
		  s2[*ind2].Y = Y; // Same particle weight as neutrals SN = 1
		  s2[*ind2].r = pa[n].p.r;
		  (*ind2)++;
		}

	    }
	  
	  n_lost++;
	  continue; 
	}

      else if ( pa[n].p.z < zmin ) 	
	{
	  jr = int (pa[n].p.r);
	  *n_wall1 += 1;
	  *I_wall1 += 2*jr+1;
	  *I_C += 1;

	  nrg_w1 = SQU(pa[n].p.vz) + SQU(pa[n].p.vr) + SQU(pa[n].p.vt);
	  nrg_w1 *= Te/(2*SQU(cs_i));

	  if ( nrg_w1 > 23.383 )
	    {
	      eps_1 = (4.45394e-06)*nrg_w1;
	      sn_1 = 8205*3.441*sqrt(eps_1)*log(eps_1 + 2.718)/(1 + 6.355*sqrt(eps_1) + eps_1*(6.882*sqrt(eps_1) - 1.708));
	      p = 0.042*(0.2525/3.49)*(sn_1/(1 + 0.010124*0.1573*pow(eps_1,0.3)))*pow(1 - sqrt(23.383/nrg_w1),2.5);

	      // Register fluxes going through each cell for enhanced Y
	      for ( i=0; i<nr; i++ )
		{
		  if ( (i <= pa[n].p.r) && (pa[n].p.r < (i+1)) )
		    current[i] += 1.; // to be rescaled!
		}

	      // Fractional part of p handled probabilistically
	      Y = int ( p );
	      p -= Y;
	      if ( RAND <= p )
		Y++;
		
	      // Register r-coordinates of bombarding particles
	      if ( Y > 0 )
		{
		  s1_temp[k].Y = Y; // Same particle weight as neutrals, SN = 1
		  s1_temp[k].r = pa[n].p.r;
		  k++;
		  Y_sum += Y;
		}
	    }

	  // SEY = 0.5 = constant, only from ions hitting the cathode
	  // SEY with registering r-coordinates
	  // Set reasonable threshold for incident ion energy (e.g. 100 eV)
	  if ( nrg_w1 > 100 )
	    {
	      if ( RAND <= SEY_f )
		{
		  se[*inde].Y = SEY_i + 1;
		  se[*inde].r = pa[n].p.r;
		  (*inde)++;
		}
	      else if ( SEY_i > 0 )
		{
		  se[*inde].Y = SEY_i;
		  se[*inde].r = pa[n].p.r;
		  (*inde)++;
		}
	    }
	    

	  n_lost++;
	  continue; 
	}

      pa[n-n_lost].p = pa[n].p;

    }
  *np -= n_lost; 


  // Enhanced yield?
  int check_enh = 0;

  for ( i=0; i<nr; i++ )
    {
      // Rescale current with area of cell / area Ldb^2
      current[i] /= PI * (2*i + 1) * SQU(dr);

      if ( current[i] >= threshold )
	{
	    check_enh = 1;
	    sigma = i+1;
	}
    }
  
  if ( (check_enh == 0) && (k > 0) ) // CORR! 11.1.2011
    {
      // Yamamura-Tawara fitting for Cu -> Cu
      for ( j=0; j<k; j++ )
	{
	  s1[*ind1].Y = s1_temp[j].Y;
	  s1[*ind1].r = s1_temp[j].r;
	  (*ind1)++;
	}
    }
  else if ( check_enh == 1 )
    {
      // Enhanced yield from MD
      // Generate Gaussian distribution for r-coord's 
      for ( j=0; j<1000; j++ ) // 1000/SN
	{
	  do { r = sigma * sqrt(-2.*log(RAND+1.e-20)); }
	  while ( r >= nr );

	  s1[*ind1].Y = 1;
	  s1[*ind1].r = r;
	  (*ind1)++;
	}

      // Print control message
      if ( Y_sum > 1000 ) // 1000/SN
	{
	  printf("*** UNDERESTIMATED *** Yamamura sputtering yield Y_sum ( %d ) > enhanced sputtering yield ( %d )! \n", Y_sum, 1000 );
	  fflush( stdout );
	}

    }

  // empty arrays
  for (i=0; i<NSY; i++)
    {
      s1_temp[i].r = 0.;
      s1_temp[i].Y = 0;
    }

   

}



void  reflect_particles( Particle pa[], int np, double  zmin, double zmax, double rmax )
{
  int n;
  double delta;

  for ( n=0; n<np; n++ )
    {
     
      if ( pa[n].p.r >= rmax ) 	
	{
	  delta = pa[n].p.r - rmax;
	  pa[n].p.r -= 2.*delta;
	  //pa[n].p.r = 2.*rmax - pa[n].p.r;
	  pa[n].p.vr *= -1.;
	  continue; 
	}
      else if ( pa[n].p.z >= zmax ) 		
	{
	  delta = pa[n].p.z - zmax;
	  pa[n].p.z -= 2.*delta;
	  //pa[n].p.z = 2.*zmax - pa[n].p.z;
	  pa[n].p.vz *= -1.;
	  continue; 
	}
      else if ( pa[n].p.z < zmin ) 	
	{
	  delta = zmin - pa[n].p.z;
	  pa[n].p.z += 2.*delta;
	  //pa[n].p.z = 2.*zmin - pa[n].p.z;
	  pa[n].p.vz *= -1.;
	  continue; 
	}

    }
}   


void  periodic_BC( Particle pa[], int np, double  zmin, double zmax, double rmax )
{
  int n;
  double delta; 

  for ( n=0; n<np; n++ )
    {
	
      // Reflection at infinity
      if ( pa[n].p.r >= rmax ) 	
	{
	  delta = pa[n].p.r - rmax;
	  pa[n].p.r -= 2.*delta;
	  pa[n].p.vr *= -1.;
	  continue; 
	}

      else if ( pa[n].p.z >= zmax ) 		
	{
	  pa[n].p.z += zmin - zmax;
	  continue; 
	}
      else if ( pa[n].p.z < zmin ) 	
	{
	  pa[n].p.z += zmax - zmin;
	  continue; 
	}
     

    }
}   


void  remove_simply_n( Particle pa[], int *np, double zmin, double zmax, double rmax, int *n_inf )
{
  int n, n_lost = 0;

  for ( n=0; n<*np; n++ )
    {
      // Remove
      if ( pa[n].p.r >= rmax ) 	
	{
	  *n_inf +=1;
	  n_lost++;
	  continue; 
	}
      // Reflect
      else if ( pa[n].p.z >= zmax ) 		
	{
	  pa[n].p.z = 2.*zmax - pa[n].p.z;
	  pa[n].p.vz *= -1.;
	  continue; 
	}
      // Reflect
      else if ( pa[n].p.z < zmin ) 	
	{
	  pa[n].p.z = 2.*zmin - pa[n].p.z;
	  pa[n].p.vz *= -1.;
	  continue; 
	}

      pa[n-n_lost].p = pa[n].p;

    }
  *np -= n_lost;
}   



void  remove_simply_i( Particle pa[], int *np, double zmin, double zmax, double rmax, 
		       int *sput1, int *sput2, int *n_wall1, int *n_wall2, int *n_inf )
{
  int n, n_lost = 0;

  for ( n=0; n<*np; n++ )
    {
      // Remove     
      if ( pa[n].p.r >= rmax ) 	
	{
	  *n_inf +=1;
	  n_lost++;
	  continue; 
	}
      // Sputter
      else if ( pa[n].p.z >= zmax ) 		
	{
	  *sput2 += 1;
	  *n_wall2 += 1;
	  n_lost++;
	  continue; 
	}
      // Sputter
      else if ( pa[n].p.z < zmin ) 	
	{
	  *sput1 += 1;
	  *n_wall1 +=1;
	  n_lost++;
	  continue; 
	}

      pa[n-n_lost].p = pa[n].p;

    }
  *np -= n_lost;
}   



void  remove_simply_e( Particle pa[], int *np, double zmin, double zmax, double rmax,
		       int *sput1, int *sput2, int *n_wall1, int *n_wall2, int *n_inf )
{
  int n, n_lost = 0;
  double tmp1 = 0, tmp2 = 0;

  for ( n=0; n<*np; n++ )
    {
      // Remove
      if ( pa[n].p.r >= rmax ) 	
	{
	  *n_inf +=1;
	  n_lost++;
	  continue; 
	}
      // Sputter
      else if ( pa[n].p.z >= zmax ) 		
	{
	  tmp2 += .01;
	  *n_wall2 += 1;
	  n_lost++;
	  continue; 
	}
      // Sputter
      else if ( pa[n].p.z < zmin ) 	
	{
	  tmp1 += .01;
	  *n_wall1 +=1;
	  n_lost++;
	  continue; 
	}

      pa[n-n_lost].p = pa[n].p;

    }
  *np -= n_lost;

  int Y1 = int ( tmp1 );
  tmp1 -= Y1;
  if ( RAND <= tmp1 )
    Y1++;
  *sput1 += Y1;

  int Y2 = int ( tmp2 );
  tmp2 -= Y2;
  if ( RAND <= tmp2 )
    Y2++;
  *sput2 += Y2;

}   



/***********************
 Initialise particles
***********************/

void  init_uniformly( Particle pa[], int *np, double vt, double rmax_inj, double sign, 
		      double ndb, double dz, double nz, double zmin, double zmax, double rmax )
{

  //double rmax_inj = 10; 
  double r1,r2;
  int i,m;

  m = (int) (PI*dz*dz*dz*rmax_inj*rmax_inj*nz*ndb);
  //printf("Initialise %d electrons and ions \n",nr_e);
  
  for (i=0; i<m; i++)
    {
      // OLD SCHEME  
      do { r1 = sqrt(-2.*log(RAND+1.e-20)); }
      while( r1 > 5. );
      r2 = RAND * TWOPI;
      pa[i].p.vz = r1*sin(r2)*vt;
      pa[i].p.vt = r1*cos(r2)*vt;
      
      do { r1 = sqrt(-2.*log(RAND+1.e-20)); } 
      while( r1 > 5. );
      pa[i].p.vr = r1*cos(RAND*TWOPI)*vt;
      
      pa[i].p.z = zmin + zmax*RAND;
      pa[i].p.r = rmax_inj*sqrt(RAND);
      pa[i].p.m = 1;
      
      // NEW SCHEME
      //r1 = sqrt(-2.*log(RAND+1.e-20)); 
      //r2 = RAND * TWOPI;
      //pa[i].p.vr = r1*cos(r2)*vt;
      //pa[i].p.vt = r1*sin(r2)*vt;
      //pa[i].p.vz = sign*sqrt(-2.*log(RAND+1.e-20))*vt; // from cathode -> anode, if +1
      
      //pa[i].p.z = zmin + pa[i].p.vz;
      //pa[i].p.r = rmax_inj*sqrt(RAND);
    }

  *np += m;
}


/* 
**********************
INJECT PARTICLES
**********************
*/

void  inject_arc_e( Particle pa[], int *np, int step, double vt, double zmin, double zmax, int nr,
		    double rem_th, double rem, int *emitted, int *inde, Sput se[], int *erosion, int *check, 
		    int *el_inwards_emis1, int *el_inwards_emisf, int *el_inwards_sey1, 
		    int *el_inwards_w1, int *I_wall1, int *I_C,
		    double *B, double B_f, double const E[], double Te, double n0, double Debye )
{
 int n, n_lost = 0;
 int i, k, jr, jj;
 int totY = 0;
 double r1, r2;
 double j, Eloc, field;
 int tmp;
 double jFN[nr];
 

 // FN AT FIELD EMITTER
 field = E[0];
 if (field < 0.) 
   {
     // Fowler-Nordheim with Wang-Loew approximation
     // W-L: v(y) ~ 0.956 - 1.062*(3.7947e-5)^2*Eloc/(4.5*4.5)
     // beta=dynamic, work fct=4.5eV, j in A/cm^2, E in V/m
     // E = - field * 1e-7 * 2* Te * dz / (sqrt(552635 * Te / n0) * SQU( Omega_pe )); // in GV/m
     // j = 1e-4 *1.54e6 / 4.5 * SQU (Eloc) * exp (-6.8309e3 * pow(4.5,1.5) * (0.956 - 7.552e-6 * Eloc)/ (Eloc)); // in A/cm^2

     Eloc = - 2.69036254e-10*dz/SQU(Omega_pe)*sqrt(Te*n0)*field*(*B); // rescale the field to GV/m
     // Protect against numerical fluctuations
     if ( Eloc > 12. )
       Eloc = 12.;

     j = 4.7133e9 * SQU(Eloc) * exp(-62.338/Eloc); // in A/cm^2
     j *= (Debye*step*Omega_pe*SQU(rem_th*dz))/(6.7193e-12*n0*sqrt(Te)); //dimless, 2D: rem^2
   }

 else { j=0.; }

 // FN OUTSIDE THE FIELD EMITTER
 for ( jj=0; jj<nr; jj++ )
   {
     field = E[jj+1];
     if (field < 0.) 
       {
	 Eloc = - 2.69036254e-10*dz/SQU(Omega_pe)*sqrt(Te*n0)*field*B_f; // B=B_f, rescale the field to GV/m

	 if ( Eloc > 12. )
	   Eloc = 12.;

	 jFN[jj] = 4.7133e9 * SQU(Eloc) * exp(-62.338/Eloc); // in A/cm^2

	 if ( rem <= jj )
	   jFN[jj] *= (Debye*step*Omega_pe*PI*SQU(dz)*(2*jj+1))/(6.7193e-12*n0*sqrt(Te)); //dimless
	 else if ( (rem > jj) && (rem < jj+1) )
	   jFN[jj] *= Debye*step*Omega_pe*PI*SQU(dz)*( SQU(jj+1)-SQU(rem) )/(6.7193e-12*n0*sqrt(Te)); //dimless
	 else 
	   jFN[jj] = 0.;
       }

     else { jFN[jj]=0.; }
   }

 //2D: use Rem (=emission radius) instead of 20nm
 if (*check == 0) 
   {
     // for 20 nm: B -= erosion * (n0 * pow(552635*Te/n0,1.5) / Debye) * 1.9121e-6 / 4; 
     *B -= (*erosion) * n0 / ( Debye * pow(rem_th*dz,3)) * 3.824131e-24; 
     *erosion = 0;
     if (*B < B_f) { *B = B_f; }
   }
 //else 
 //  *B = 1;
   


 // Check indeces
 if ( *inde > NSY ) 
   {
     printf ("*** ERROR *** in sputtering yield index: inde ( %d ) ! \n",*inde);
     fflush( stdout );
   }

 // injection only from cathode
 
 // SEY: inject with incident r-coordinates, empty SEY variables!
 for ( k=0; k<*inde; k++ )
   {
     for ( i=0; i<se[k].Y; i++ )
       {
	 do { r1 = sqrt(-2.*log(RAND+1.e-20)); }
	 while( r1 > 5. );
	 //r2 = RAND * PI;
	 
	 //pa[*np+totY].p.vz = r1*sin(r2)*vt;
	 //pa[*np+totY].p.vt = r1*cos(r2)*vt;
	 
	 // Gaussian scheme
	 r2 = RAND * TWOPI;
	 pa[*np+totY].p.vr = r1*cos(r2)*vt;
	 pa[*np+totY].p.vt = r1*sin(r2)*vt;

	 do { r1 = sqrt(-2.*log(RAND+1.e-20)); } 
	 while( r1 > 5. );
	 
	 //pa[*np+totY].p.vr = r1*cos(RAND*TWOPI)*vt;
	 
	 // Gaussian scheme
	 pa[*np+totY].p.vz = r1*vt; 
	 
	 pa[*np+totY].p.z = zmin + pa[*np+totY].p.vz;
	 pa[*np+totY].p.r = se[k].r + pa[*np+totY].p.vr;
	 if ( pa[*np+totY].p.r < 0 )
	   pa[*np+totY].p.r = 1.e-20;
	 
	 pa[*np+totY].p.m = 1;
	 
	 jr = int (pa[*np+totY].p.r);
	 *I_wall1 += 2*jr+1;
	 totY++;
       }
   }
 *np += totY;

 *el_inwards_sey1 += totY; 
 *el_inwards_w1 += totY; 
 *I_C += totY;

 // empty arrays
 *inde = 0;
 for (i=0; i<NSY; i++)
   {
     se[i].r = 0.;
     se[i].Y = 0;
   }


 // FN INJECTION FROM FIELD EMITTER
 // FE: inject with flat distribution over the emission radius
 tmp = int (j);
 j -= tmp;
 if ( RAND <= j )
   tmp++;

   
 for  ( k=0; k<tmp; k++ )
   { 
     do { r1 = sqrt(-2.*log(RAND+1.e-20)); }
     while( r1 > 5. );
     //r2 = RAND * PI;
          	
     //pa[*np+k].p.vz = r1*sin(r2)*vt;
     //pa[*np+k].p.vt = r1*cos(r2)*vt;

     // Gaussian scheme
     r2 = RAND * TWOPI;
     pa[*np+k].p.vr = r1*cos(r2)*vt;
     pa[*np+k].p.vt = r1*sin(r2)*vt;
     
     do { r1 = sqrt(-2.*log(RAND+1.e-20)); } 
     while( r1 > 5. );
     //pa[*np+k].p.vr = r1*cos(RAND*TWOPI)*vt;
     
     // Gaussian scheme
     pa[*np+k].p.vz = r1*vt; 

     pa[*np+k].p.z = zmin + pa[*np+k].p.vz;
     pa[*np+k].p.r = rem * sqrt(RAND) + pa[*np+k].p.vr; // 15.2.2011

     if ( pa[*np+k].p.r < 0 )
       pa[*np+k].p.r = 1.e-20;
     pa[*np+k].p.m = 1;
    
     jr = int (pa[*np+k].p.r);
     *I_wall1 += 2*jr+1;
   }
 *np += tmp;
       
 *el_inwards_emis1 += tmp; 
 *el_inwards_w1 += tmp; 
 *emitted += tmp;
 *I_C += tmp;


 // FN INJECTION OUTSIDE THE FIELD EMITTER
 for ( jj=0; jj<nr; jj++ )
   {
     if ( jFN[jj] > 0 )
       {
	 // THEN INJECT 
	 tmp = int (jFN[jj]);
	 jFN[jj] -= tmp;
	 if ( RAND <= jFN[jj] )
	   tmp++;
	 
	 for  ( k=0; k<tmp; k++ )
	   { 
	     do { r1 = sqrt(-2.*log(RAND+1.e-20)); }
	     while( r1 > 5. );
	     //r2 = RAND * PI;
	     
	     //pa[*np+k].p.vz = r1*sin(r2)*vt;
	     //pa[*np+k].p.vt = r1*cos(r2)*vt;

	     // Gaussian scheme
	     r2 = RAND * TWOPI;
	     pa[*np+k].p.vr = r1*cos(r2)*vt; 
	     pa[*np+k].p.vt = r1*sin(r2)*vt; 
	       
	     do { r1 = sqrt(-2.*log(RAND+1.e-20)); } 
	     while( r1 > 5. );
	     //pa[*np+k].p.vr = r1*cos(RAND*TWOPI)*vt;
	     
	     // Gaussian scheme
	     pa[*np+k].p.vz = r1*vt; 
	     
	     pa[*np+k].p.z = zmin + pa[*np+k].p.vz;
	     if ( rem <= jj )
	       pa[*np+k].p.r = jj + RAND + pa[*np+k].p.vr;
	     else if ( (rem > jj) && (rem < jj+1) )
	       pa[*np+k].p.r = rem + (jj + 1 - rem) * RAND + pa[*np+k].p.vr;

	     if ( pa[*np+k].p.r < 0 )
	       pa[*np+k].p.r = 1.e-20;
	     pa[*np+k].p.m = 1;
	     
	     jr = int (pa[*np+k].p.r);
	     *I_wall1 += 2*jr+1;
	   }
	 *np += tmp;
	 
	 *el_inwards_emisf += tmp; 
	 *el_inwards_w1 += tmp; 
	 *emitted += tmp;
	 *I_C += tmp;
       }
   }

}


void  inject_arc_n( Particle  pa[], int *np, double vt, double zmin, double zmax, 
		    double rem_th, double rem, double ratio,
		    int *ind1, Sput s1[], int *ind2, Sput s2[], int *Cu_inwards_w1, int *Cu_inwards_w2, 
		    int *Cu_inwards_sput1, int *Cu_inwards_sput2, int *Cu_inwards_evap1, 
		    int *emitted, int *erosion )
{

 int n, n_lost = 0;
 int tot1 = 0, tot2 = 0;
 int m, el_eroded, i, k;
 double r1, r2;

 // evaporation of neutrals
 double tmp = ratio * (*emitted);
 int n2inject_evap = int ( tmp );
 tmp -= n2inject_evap;
 if ( RAND <= tmp )
   n2inject_evap++;
 *emitted = 0;


 // first  wall (x= 0)   
 
 // neutral evaporation
 m = n2inject_evap;
 
 for  ( k=0; k<m; k++ )
   { 
     do { r1 = sqrt(-2.*log(RAND+1.e-20)); }
     while( r1 > 5. );
     //r2 = RAND * PI;
     
     //pa[*np+k].p.vz = r1*sin(r2)*vt;
     //pa[*np+k].p.vt = r1*cos(r2)*vt/100.; // suppress perp components

     // Gaussian scheme
     r2 = RAND * TWOPI;
     pa[*np+k].p.vr = r1*cos(r2)*vt; // do not suppress
     pa[*np+k].p.vt = r1*sin(r2)*vt; // do not suppress

     do { r1 = sqrt(-2.*log(RAND+1.e-20)); } 
     while( r1 > 5. );
     //pa[*np+k].p.vr = r1*cos(RAND*TWOPI)*vt/100.; // suppress perp components
     
     // Gaussian scheme
     pa[*np+k].p.vz = r1*vt; 
     
     pa[*np+k].p.z = zmin + pa[*np+k].p.vz;
     pa[*np+k].p.r = rem * sqrt(RAND) + pa[*np+k].p.vr; // 15.2.2011

     if ( pa[*np+k].p.r < 0 )
       pa[*np+k].p.r = 1.e-20;
     pa[*np+k].p.m = 1; //SN;
     
   }
 *np += m;
	 
 // Check indeces
 if ( (*ind1 > NSY) || (*ind2 > NSY) ) 
   {
     printf ("*** ERROR *** in sputtering yield indeces: ind1 ( %d ), ind2 ( %d ) ! \n",*ind1,*ind2);
     fflush( stdout );
   }


 // sputtering
 m = *ind1;

 for ( k=0; k<m; k++ )
   {
     for ( i=0; i<s1[k].Y; i++ )
       {
	 do { r1 = sqrt(-2.*log(RAND+1.e-20)); }
	 while( r1 > 5. );
	 //r2 = RAND * PI;

	 //pa[*np+tot1].p.vz = r1*sin(r2)*vt;
	 //pa[*np+tot1].p.vt = r1*cos(r2)*vt;

	 // Gaussian scheme
	 r2 = RAND * TWOPI;
	 pa[*np+tot1].p.vr = r1*cos(r2)*vt; 
	 pa[*np+tot1].p.vt = r1*sin(r2)*vt; 

	 do { r1 = sqrt(-2.*log(RAND+1.e-20)); } 
	 while( r1 > 5. );
	 
	 //pa[*np+tot1].p.vr = r1*cos(RAND*TWOPI)*vt;
	 
	 // Gaussian scheme
	 pa[*np+tot1].p.vz = r1*vt; 
	 
	 pa[*np+tot1].p.z = zmin + pa[*np+tot1].p.vz;
	 pa[*np+tot1].p.r = s1[k].r + pa[*np+tot1].p.vr; 

	 if ( pa[*np+tot1].p.r < 0 )
	   pa[*np+tot1].p.r = 1.e-20;
	 pa[*np+tot1].p.m = 1; //SN;

	 // Sum up total yield for outputting
	 tot1++;
       }
   }
 *np += tot1;

 // update currents
 *Cu_inwards_w1 += n2inject_evap + tot1;
 *Cu_inwards_sput1 += tot1;
 *Cu_inwards_evap1 += n2inject_evap;
 
 *erosion += n2inject_evap; // erosion due to evaporation, not sputtering
 
 
	 
 // second wall
 
 // only sputtering
 m = *ind2;

 for ( k=0; k<m; k++ )
   {
     for ( i=0; i<s2[k].Y; i++ )
       {
	 do { r1 = sqrt(-2.*log(RAND+1.e-20)); }
	 while( r1 > 5. );
	 //r2 = -RAND * PI;

	 //pa[*np+tot2].p.vz = r1*sin(r2)*vt;
	 //pa[*np+tot2].p.vt = r1*cos(r2)*vt;

	 // Gaussian scheme
	 r2 = RAND * TWOPI;
	 pa[*np+tot2].p.vr = r1*cos(r2)*vt; 
	 pa[*np+tot2].p.vt = r1*sin(r2)*vt; 

	 do { r1 = sqrt(-2.*log(RAND+1.e-20)); } 
	 while( r1 > 5. );
	     
	 //pa[*np+tot2].p.vr = r1*cos(RAND*TWOPI)*vt;
             
	 // Gaussian scheme
	 pa[*np+tot2].p.vz = -1.*r1*vt; 
             
	 pa[*np+tot2].p.z = zmax + pa[*np+tot2].p.vz;
	 pa[*np+tot2].p.r = s2[k].r + pa[*np+tot2].p.vr; 

	 if ( pa[*np+tot2].p.r < 0 )
	   pa[*np+tot2].p.r = 1.e-20;
	 pa[*np+tot2].p.m = 1; //SN;

	 // Sum up total yield for outputting
	 tot2++;
       }
   }
 *np += tot2;
 
 // update currents
 *Cu_inwards_sput2 += tot2;
 *Cu_inwards_w2 += tot2;	

 // empty arrays
 *ind1 = 0;
 *ind2 = 0;
 for (i=0; i<NSY; i++)
   {
     s1[i].r = 0., s2[i].r = 0.;
     s1[i].Y = 0, s2[i].Y = 0;
   }


}


void  inject_at_centre( Particle pa[], int *np, int n2inject, double vt, 
			double zmin, double zmax, double rmax, int *p_inwards_c )
{
  int m, k;
  double r1, r2;

  // Inject to inner few cells 
  double percent = 1.;  // how many percent of the cells will be injected into
  int j = 5;            // innermost j cells in r-dir
  double lsys = zmax - zmin;
  double offset = (1 - percent)*lsys/2.;
  double range = percent * lsys;

  if ( j > rmax )
    printf ("Error in injection scheme, system too small! \n");

               
         
         m = n2inject;                                             
       
	 for  (k=0; k < m; k++)
            { 
	      do { r1 = sqrt(-2.*log(RAND+1.e-20)); }
	      while( r1 > 5. );
	      r2 = RAND * TWOPI;
	      pa[*np+k].p.vz = r1*sin(r2)*vt;
	      pa[*np+k].p.vt = r1*cos(r2)*vt;
	      //pa[*np+k].p.vt = r1*cos(r2)*vt/10.; // energycons tests
	      
	      do { r1 = sqrt(-2.*log(RAND+1.e-20)); } 
	      while( r1 > 5. );
	      pa[*np+k].p.vr = r1*cos(RAND*TWOPI)*vt;
	      //pa[*np+k].p.vr = r1*cos(RAND*TWOPI)*vt/1000.; // enegycons tests
	      
	      pa[*np+k].p.z = offset + range*RAND;
	      pa[*np+k].p.r = j*RAND;
	      //pa[*np+k].p.r = 50*RAND; // enegycons tests
	      pa[*np+k].p.m =  1;
            }
            
	 *np += m;       

	 // currents
	 *p_inwards_c += m;
}


void  inject_uniformly( Particle pa[], int *np, int n2inject, double vt, 
			double zmin, double zmax, double rmax, int *p_inwards_c )
{
  int m, k;
  double r1, r2;

  // Inject to inner few cells with uniform density
  double percent = 1.;  // how many percent of the cells will be injected into
  int j = 5;            // innermost j cells in r-dir
  double lsys = zmax - zmin;
  double offset = (1 - percent)*lsys/2.;
  double range = percent * lsys;

  if ( j > rmax )
    printf ("Error in injection scheme, system too small! \n");

               
         
         m = n2inject;                                             
       
	 for  (k=0; k < m; k++)
            { 
	      do { r1 = sqrt(-2.*log(RAND+1.e-20)); }
	      while( r1 > 5. );
	      r2 = RAND * TWOPI;
	      pa[*np+k].p.vz = r1*sin(r2)*vt;
	      pa[*np+k].p.vt = r1*cos(r2)*vt;
	      
	      do { r1 = sqrt(-2.*log(RAND+1.e-20)); } 
	      while( r1 > 5. );
	      pa[*np+k].p.vr = r1*cos(RAND*TWOPI)*vt;
	      
	      pa[*np+k].p.z = offset + range*RAND;
	      pa[*np+k].p.r = j*sqrt(RAND);
	      pa[*np+k].p.m =  1;
            }
            
	 *np += m;       

	 // currents
	 *p_inwards_c += m;
}
    

void  inject_from_cathode( Particle pa[], int *np, int n2inject, double vt, 
			   double zmin, double zmax, double rmax, int *p_inwards_w1 )
{
  int m, k;
  double r1, r2;
            
         
         m = n2inject;                                             
       
	 for  (k=0; k < m; k++)
            { 
	      do { r1 = sqrt(-2.*log(RAND+1.e-20)); }
	      while( r1 > 5. );
	      r2 = RAND * PI;
	      pa[*np+k].p.vz = r1*sin(r2)*vt;
	      //pa[*np+k].p.vt = r1*cos(r2)*vt; 
	      pa[*np+k].p.vt = r1*cos(r2)*vt/10.; // suppress perp. comp.
	      
	      do { r1 = sqrt(-2.*log(RAND+1.e-20)); } 
	      while( r1 > 5. );
	      //pa[*np+k].p.vr = r1*cos(RAND*TWOPI)*vt;
	      pa[*np+k].p.vr = r1*cos(RAND*TWOPI)*vt/10.; // suppress perp. comp.
	      
	      pa[*np+k].p.z = zmin + pa[*np+k].p.vz;
	      //pa[*np+k].p.r = fabs (pa[*np+k].p.vr); // pointlike source
	      pa[*np+k].p.r = RAND; // uniform source

	      pa[*np+k].p.m =  1;
            }
            
	 *np += m;       

	 // currents
	 *p_inwards_w1 += m;
}


void  inject_from_anode( Particle pa[], int *np, int n2inject, double vt, 
			 double zmin, double zmax, double rmax, int *p_inwards_w2 )
{
  int m, k;
  double r1, r2;
            
         
         m = n2inject;                                             
       
	 for  (k=0; k < m; k++)
            { 
	      do { r1 = sqrt(-2.*log(RAND+1.e-20)); }
	      while( r1 > 5. );
	      r2 = -RAND * PI;
	      pa[*np+k].p.vz = r1*sin(r2)*vt;
	      pa[*np+k].p.vt = r1*cos(r2)*vt; 
	      //pa[*np+k].p.vt = r1*cos(r2)*vt/10.; // energycons tests
	      
	      do { r1 = sqrt(-2.*log(RAND+1.e-20)); } 
	      while( r1 > 5. );
	      pa[*np+k].p.vr = r1*cos(RAND*TWOPI)*vt;
	      //pa[*np+k].p.vr = r1*cos(RAND*TWOPI)*vt/1000.; // enegycons tests
	      
	      pa[*np+k].p.z = zmax + pa[*np+k].p.vz;
	      pa[*np+k].p.r = fabs (pa[*np+k].p.vr); 

	      pa[*np+k].p.m =  1;
            }
            
	 *np += m;       

	 // currents
	 *p_inwards_w2 += m;
}


void  inject_simply_n( Particle pa[], int *np, int n2inject, int *sput1, int *sput2, double vt, 
		       double zmin, double zmax, double rmax, int *p_inwards_w1, int *p_inwards_w2 )
{
  int m, k;
  double r1, r2;
            
         
         // Evaporation & sputtering
         m = n2inject + *sput1;                                             
       
	 for  (k=0; k < m; k++)
            { 
	      do { r1 = sqrt(-2.*log(RAND+1.e-20)); }
	      while( r1 > 5. );
	      r2 = RAND * PI;
	      pa[*np+k].p.vz = r1*sin(r2)*vt;
	      //pa[*np+k].p.vt = r1*cos(r2)*vt; 
	      pa[*np+k].p.vt = r1*cos(r2)*vt/10.; // suppress perp. comp.
	      
	      do { r1 = sqrt(-2.*log(RAND+1.e-20)); } 
	      while( r1 > 5. );
	      //pa[*np+k].p.vr = r1*cos(RAND*TWOPI)*vt;
	      pa[*np+k].p.vr = r1*cos(RAND*TWOPI)*vt/10.; // suppress perp. comp.
	      
	      pa[*np+k].p.z = zmin + pa[*np+k].p.vz;
	      pa[*np+k].p.r = rmax*sqrt(RAND); // uniform source

	      pa[*np+k].p.m =  1;
            }
            
	 *np += m;       

	 // currents
	 *p_inwards_w1 += m;
	 *sput1 = 0;

         // Only sputtering
         m = *sput2;                                             
       
	 for  (k=0; k < m; k++)
            { 
	      do { r1 = sqrt(-2.*log(RAND+1.e-20)); }
	      while( r1 > 5. );
	      r2 = -RAND * PI;
	      pa[*np+k].p.vz = r1*sin(r2)*vt;
	      pa[*np+k].p.vt = r1*cos(r2)*vt/10.; // suppress perp. comp.
	      
	      do { r1 = sqrt(-2.*log(RAND+1.e-20)); } 
	      while( r1 > 5. );
	      pa[*np+k].p.vr = r1*cos(RAND*TWOPI)*vt/10.; // suppress perp. comp.
	      
	      pa[*np+k].p.z = zmax + pa[*np+k].p.vz;
	      pa[*np+k].p.r = rmax*sqrt(RAND); //fabs (pa[*np+k].p.vr); 

	      pa[*np+k].p.m =  1;
            }
            
	 *np += m;       

	 // currents
	 *p_inwards_w2 += m;
	 *sput2 = 0;
}


void  inject_simply_e( Particle pa[], int *np, int n2inject, double vt, 
		       double zmin, double zmax, double rmax, int *p_inwards_w1 )
{
  int m, k;
  double r1, r2;
            
         
         // Emission only
         m = n2inject;                                             
       
	 for  (k=0; k < m; k++)
            { 
	      do { r1 = sqrt(-2.*log(RAND+1.e-20)); }
	      while( r1 > 5. );
	      r2 = RAND * PI;
	      pa[*np+k].p.vz = r1*sin(r2)*vt;
	      pa[*np+k].p.vt = r1*cos(r2)*vt; 
	      
	      do { r1 = sqrt(-2.*log(RAND+1.e-20)); } 
	      while( r1 > 5. );
	      pa[*np+k].p.vr = r1*cos(RAND*TWOPI)*vt;
	      
	      pa[*np+k].p.z = zmin + pa[*np+k].p.vz;
	      pa[*np+k].p.r = rmax*sqrt(RAND); // uniform source

	      pa[*np+k].p.m =  1;
            }
            
	 *np += m;       

	 // currents
	 *p_inwards_w1 += m;


}


