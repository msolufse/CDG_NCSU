
/***************************************************************************/
/*                                                                         */
/* The junction.h header program                                            */
/*  Version: 1.0                                                           */
/*  Date: 24 March 2020                                                    */
/*                                                                         */
/*  Primary Authors: M.J. Colebank & J. Mackenzie                          */
/*  Key Contributers: M.U. Qureshi & M.S. Olufsen                          */
/*  Institute: North Carolina State University, Raleigh, NC, USA           */
/*  Funded by: NSF-DMS # 1615820 & AHA #19PRE34380459                      */
/*                                                                         */
/***************************************************************************/


#ifndef _JUNCTION_H
#define _JUNCTION_H

#include <cmath>
#include "tools.h"

// Updates bifurcation conditions. Uses daughter vessels, and should
// only be called when such exist.

void junction   (double theta, double gamma, Tube *Arteries[], int parent, int WM);

//void bound_monf (double theta, double gamma, Tube *Arteries[], int parent, int D1, int WM);

void bound_screen_linesearch (double theta, double gamma, Tube *Arteries[], int parent, int D1, int WM);
void bound_screen_linesearch_fvec (double *fvec, double *xb, double theta, double gamma, double *k1, double *k2, double *k3, double *k4, Tube *PV,Tube *B1, int WM);
void bound_screen_linesearch_fjac (double *fjac[], double *xb, double theta, double gamma, Tube *PV,Tube *B1, int WM);


void bound_bif  (double theta, double gamma, Tube *Arteries[], int parent, int D1, int D2, int WM);
void bound_bif_fvec (double *fvec, double *xb, double theta, double gamma, double *k1, double *k2, double *k3, double *k4, Tube *PV,Tube *B1, Tube *B2, int WM);
void bound_bif_fjac (double *fjac[], double *xb, double theta, double gamma, Tube *PV,Tube *B1, Tube *B2, int WM);


void bound_sten_bif  (double theta, double gamma, Tube *Arteries[], int parent, int D1, int D2, int WM);
void bound_sten_fvec (double *fvec, double *xb, double theta, double gamma, double *k1, double *k2, double *k3, double *k4, double dp1, double dp2, double dp3, Tube *PV,Tube *B1, Tube *B2, int WM);
void bound_sten_fjac (double *fjac[], double *xb, double theta, double gamma, double dp1, double dp2, double dp3, Tube *PV, Tube *B1, Tube *B2, int WM);


void bound_trif (double theta, double gamma, Tube *Arteries[], int parent, int D1, int D2, int D3, int WM);
void bound_trif_fvec ();
void bound_trif_fjac ();

// In order to ensure a more efficient execution of the program the following
// functions is made as in-line functions.

// A function returning the Friction of the system. The definition of this
// function is given according to the derivation in the mathematical model.
// The constant cst, determines the amount of damping in the system.
inline double F (double Q, double A)
{
    double tmp1 = -2.0*sqrt(M_PI)*Q; 
    double tmp2 = bound_thick*Re*sqrt(A);        
    double tmp3 = tmp1/tmp2;

  return(tmp3);
}

inline double dFdQ (double A)
{
    return((-2.0*sqrt(M_PI))/(bound_thick*Re*sqrt(A))); 
}

inline double dFdA (double Q, double A)
{
    double tmp1 = Q*sqrt(M_PI);
    double tmp2 = bound_thick*Re*sqrt(cu(A));
    return(tmp1/tmp2); 
}


#endif
