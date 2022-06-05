/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Disk-Planet interaction problem.

  Simulate the interaction of a planet embedded in a disk as described
  in section 3.4 of Mignone et al., A&A (2012) 545, A152.
  This test is a nice benchmark for the FARGO module and the
  \c ROTATING_FRAME switch.
  For testing-purposes  no viscosity is used here.
  The initial condition consists of a locally isothermal configuration
  with temperature profile \f$\propto T^{-1}\f$  yielding a disk vertical
  height to radius of \c 0.05.
  The gravitational potential due to the presence of the star and the
  planet is defined in BodyForcePotential() function.

  The conventions used throught the implementation are the following:
 
  - \c r  = spherical radius
  - \c R  = cylindrical radius
  - \c z  = cylindrical height
  - \c th = meridional angle
  
  The test can be carried out in polar (2D or 3D) or spherical (3D)
  coordinates and the following parameters determine the initial configuration:
   
  -# <tt>g_inputParam[Mstar]</tt>: controls the star mass (in solar masses)
  -# <tt>g_inputParam[Mdisk]</tt>: controls the disk mass (in solar masses)
  -# <tt>g_inputParam[Mplanet]</tt>: sets the planet mass (in earth masses)
  -# <tt>g_inputParam[Viscosity]</tt>: sets the amount of viscosity


  Computation can be carried in the rotating or in the observer's frame
  of reference (\c ROTATING_FRAME to \c YES or \c NO, respectively).
  In particular:

  - Configurations #01 and #02 are in 2D polar coordinates without and with
    FARGO, in the rotating frame.
  - Configurations #03, #04 and #05 are in spherical 3D coordinates with
    and without FARGO in the rotating frame
  - Configurations #06 and #07 are in 2D polar coordinates but in the
    observer's frame
  - Configuration #08 employs static AMR  (grid levels are spatial
    dependent but not dependent on time) in the rotating frame.
    

  \image html hd_disk_planet.08.png "Density map for configuration #08 using AMR" width=1cm

  \author A. Mignone (mignone@to.infn.it)
  \date   Jul 19, 2019

  \b References:
     - "A Conservative orbital advection scheme for simulations
        of magnetized shear flows with the PLUTO Code"
        Mignone et al, A&A (2012)
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#define MIN_DENSITY 1e-8

static void NormalizeDensity (const Data *d, Grid *g);
#if ROTATING_FRAME == NO
 #define g_OmegaZ  0.0
#endif

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  double r, th, R, z, H, OmegaK, cs, RHO0, sigma, c1, c2, deltar;
  double phi;
  double scrh;
  double ex;
  double pert;
  double dphi, dp;

 #if EOS == IDEAL
  g_gamma = 5./3.;
  #endif


  #if ROTATING_FRAME == YES
  g_OmegaZ  = sqrt(1.0 + g_inputParam[Mplanet]/g_inputParam[Mstar]*CONST_Mearth/CONST_Msun);
  g_OmegaZ *= 2.0*CONST_PI;
  #endif
  
#if GEOMETRY == POLAR
  R  = x1;
  #if DIMENSIONS == 2
  z  = 0.0;
  r  = R;
  th = 0.5*CONST_PI;
  #else
  z  = x3;
  r  = sqrt(R*R + z*z);
  th = atan2(R,z);
  #endif
#elif GEOMETRY == SPHERICAL-
  r  = x1;
  th = x2;
  R  = r*sin(th);
  z  = r*cos(th);
#endif
 
  H      = 0.05*R;
  OmegaK = 1/(R*sqrt(R));  

  cs     = 0.05; /*0.05*/
  RHO0   = 1.0;
  c1     = 1.0;
  c2     = 0.01;
  deltar = 0.03;
  
  sigma   = RHO0/r;
  scrh    = (r-1.0)/deltar;
  ex      = exp(-0.5*scrh*scrh) ;
  pert    = c1*(1+c2*cos(3*x2));
  
  us[RHO] = sigma;
  us[RHO] += sigma*c1*(1+c2*cos(3*x2))*exp(-0.5*scrh*scrh);

  dphi = 1./(r*r);
  dp = cs*cs*( -(sigma/r)* -(sigma/r)*(1+pert*ex)-sigma*(pert*ex/deltar)  );

  us[VX1] = us[VX3] = 0.0;
  us[VX2] = sqrt(r*(dphi+(dp/us[RHO])));



  /*us[iVPHI] = R*(OmegaK - g_OmegaZ);*/
  #ifdef FARGO
  us[FARGO_W] = R*(OmegaK - g_OmegaZ);
  us[iVPHI]   = 0.0;
  #endif  
  #if EOS == IDEAL
  us[PRS] = us[RHO]*cs*cs;
  #elif EOS == ISOTHERMAL
  g_isoSoundSpeed = 0.05;
  #endif

#if DUST_FLUID == YES
  us[RHO_D] = 1.e-4;
  us[VX1_D] = us[VX2_D] = us[VX3_D] = 0.0;
  us[VX2_D] = us[iVPHI];
#endif

#if PHYSICS == MHD
  us[BX1] = us[BX3] = 0;
  us[BX2] = 0.01776/r;    /* B_phi = B0/r beta = 50*/
#endif
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{

}
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/* 
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double *x1, *x2, *x3, R, OmegaK, v[256];
  static int do_once = 1;
  
  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  #if DIMENSIONS == 3
  if (side == 0){
    if (do_once){
      NormalizeDensity(d, grid);
      do_once = 0;
    }
  }
  #endif

  if (side == X1_BEG){
    X1_BEG_LOOP(k,j,i){
      NVAR_LOOP(nv) d->Vc[nv][k][j][i] = d->Vc[nv][k][j][2*IBEG - i - 1];
      d->Vc[VX1][k][j][i] *= -1.0;
      #if GEOMETRY == POLAR
      R = x1[i];
      #elif GEOMETRY == SPHERICAL
      R = x1[i]*sin(x2[j]);
      #endif
      OmegaK = 2.0*CONST_PI/(R*sqrt(R));
      
      d->Vc[iVPHI][k][j][i] = R*(OmegaK - g_OmegaZ);
      #ifdef FARGO
      d->Vc[iVPHI][k][j][i] = 0.0;
//d->Vc[iVPHI][k][j][i] = d->Vc[iVPHI][k][j][IBEG];
      #endif      
      #if DUST_FLUID == YES      
      d->Vc[VX2_D][k][j][i] = d->Vc[iVPHI][k][j][i];
      #endif      
    }
  }

  if (side == X1_END){
    X1_END_LOOP(k,j,i){
      NVAR_LOOP(nv)  d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IEND];
      #if GEOMETRY == POLAR
      R = x1[i];
      #elif GEOMETRY == SPHERICAL
      R = x1[i]*sin(x2[j]);
      d->Vc[iVR][k][j][i]  = 0.0;
      d->Vc[iVTH][k][j][i] = 0.0;
      #endif
      OmegaK = 2.0*CONST_PI/(R*sqrt(R));
      d->Vc[iVPHI][k][j][i] = R*(OmegaK - g_OmegaZ);
      #ifdef FARGO
      d->Vc[iVPHI][k][j][i] = 0.0;
      #endif      
      #if DUST_FLUID == YES      
      d->Vc[VX2_D][k][j][i] = d->Vc[iVPHI][k][j][i];
      #endif      
    }
  }
}

/* ************************************************************** */
void NormalizeDensity (const Data *d, Grid *grid)
/*
 *
 * Normalize density and pressure as   rho -> K*rho, where
 *
 *   K = M/(\sum rho*dV)
 *
 **************************************************************** */
{
  
  
}

#if (BODY_FORCE & VECTOR)
/* ************************************************************************ */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
  g[IDIR] = -g_inputParam[Mstar]*CONST_G/(x1*x1);
  g[JDIR] = 0.0; 
  g[KDIR] = 0.0;
}
#endif

#if (BODY_FORCE & POTENTIAL)
/* ************************************************************************ */
double BodyForcePotential(double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
  double d, R, r, z, th, x, y, phiplanet, rsm;
  double xp, yp, t, phi;

#if GEOMETRY == POLAR
  R  = x1;
  #if DIMENSIONS == 2
  z  = 0.0;
  r  = R;
  th = 0.5*CONST_PI;
  #else
  z  = x3;
  r  = sqrt(R*R + z*z);
  th = atan2(R,z);
  #endif
  
  x  = R*cos(x2);
  y  = R*sin(x2);

#elif (GEOMETRY == SPHERICAL)
  r  = x1;
  th = x2;
  R = r*sin(th);
  z = r*cos(th);
  x = r*sin(th)*cos(x3);
  y = r*sin(th)*sin(x3);
#endif

/* ---------------------------------------------
             planet position
   --------------------------------------------- */

#if ROTATING_FRAME == NO
  
  phi  = - 1./x1;

  return phi;

#endif
}
#endif
