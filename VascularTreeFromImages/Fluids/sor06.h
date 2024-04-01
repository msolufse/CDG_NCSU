/***************************************************************************/
/*                                                                         */
/*  Program: sor06.h                                                       */
/*  Version: 2.0                                                           */
/*  By: Mette Olufsen, Math-Tech                                           */
/*  Date: 14. Jan. 1997                                                    */
/*                                                                         */
/*  This header file defines the global parameters                         */
/*                                                                         */
/***************************************************************************/

// $Id: sor06.h,v 1.8 2010-10-20 15:38:28 mette Exp $

#ifndef _SOR06_H
#define _SOR06_H

#include <cmath>
#include "tools.h"

int    nbrves; // Number of vessels in the tree.

int    tmstps = 8192, // The number of timesteps per period, must match number of points in Qin
       plts   = 512; // Number of plots per period.

const char *CO_filename = "Qin.dat"; // Inflow data

int max_D = 4; // For passing bifurcation pointer


const int WALL_MODEL = 2; // Define which wall model to use.

double conv   = 1333.220,  // Conversion from mmHg to SI-units.
       rho    = 1.055,    // Density of blood [g/cm^3].
       mu     = 0.032,   // Viscosity of blood [g/cm/s].
       mu_pl  = mu,  // Viscosity of blood [g/cm/s].
       nu     = mu/rho,

       Tper   = 0.9, // The period of one heart beat [s].

       Lr     = 1.0, // Characteristic radius of the vessels in the tree [cm].
       Lr2    = sq(Lr), // The squared radius [cm2].
       Lr3    = cu(Lr), // The radius to the third power [cm^3].
       g      = 981.0,  // The gravitational force [cm/s^2].
       q      = 10.0*Lr2, // The characteristic flow [cm^3/s].
       Fr2    = sq(q)/g/pow(Lr,5), // The squared Froudes number.
       Re     = q*rho/mu/Lr, // Reynolds number.
       bound_thick = sqrt(nu*Tper/(2.0*M_PI))/Lr, // Boundary layer thickness 
       Period = Tper*q/Lr3,   // The dimensionless period.
       k      = Period/tmstps, // Length of a timestep.
       Deltat = Period/plts,  // Interval between each point plot.
       p0     = 0.0*rho/g/Lr*conv, // Ensures a certain diastolic pressure.
       
    //  for stenosis
    inert_perm = 0.10,
    ALF = 1.0e-4,

       *fjac[24],   // Work space used by bound_bif.
       *fjac_ST[4], fvec_ST, //for Lax-Wendroff in bound_right
       xr, f, df;  // Work space used by bound_right.

#endif
