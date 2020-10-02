/* The sor06.C main program */

// $Id: sor06.C,v 1.10 2005/10/14 18:05:59 heine Exp $

#include "sor06.h"
#include "tools.h"
#include "arteries.h"

extern "C"  void impedance_init_driver_(int *tmstps);

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>

using namespace std;

// The vessel-network is specified and initialized, and the flow and
// pressures are to be determed. All global constants must be defined
// in the header file. That is the number of dots per cm, and the number
// of vessels.
int main(int argc, char *argv[])
{
  double tstart, tend, finaltime;

  double rm  = 0.04;
    
    double f1   = 0;//1.99925e07;// 8.99925e07;
    double f2   = 1;//-25.5267;
//    double f3   = 0.1*8.7401569e05;// 1.95*8.65e05;
//
    //==========adjustment in resistances (WK parameters) for the control========

    double Len, rtop, rbot, f3, R1, R2, CT, L;
    int HB, cycles, id, nv;
        
    if (argc != 13) //argv[0] is the name of the program, here sor06
    {
        printf("Not enough input arguments, noargc %d and they are %s\n", argc, argv[0]);
        return 1;
    }
    Len = atof(argv[1]);
    rtop = atof(argv[2]);
    rbot = atof(argv[3]);
    f3 = atof(argv[4]);
    R1 = atof(argv[5]);
    R2 = atof(argv[6]);
    CT = atof(argv[7]);
    L  = atof(argv[8]);
    nv = atof(argv[9]);
    HB = atof(argv[10]);
    cycles = atof(argv[11]);
    id = atof(argv[12]);
    

    char namepu1 [20];
    
    sprintf(namepu1, "pu1_%d.2d", id);
    
    FILE *fpu1 = fopen (namepu1, "w");

    
  // Workspace used by bound_right

  // Workspace used by bound_bif
  for(int i=0; i<18; i++) fjac[i] = new double[18];

//  clock_t c1 = clock();        // Only used when timing the program.
  nbrves    = nv;             // The number of vessels in the network.
  tstart    = 0.0;            // Starting time.
  finaltime = HB*Period;       // Final end-time during a simulation.
  tend      = (HB-cycles)*Period;       // Timestep before the first plot-point
                              // is reached.

  // The number of vessels in the network is given when the governing array of
  // vessels is declared.

  impedance_init_driver_(&tmstps);

  Tube   *Arteries[nbrves];                    // Array of blood vessels.
    
    
    // Initialization of the Arteries.

   
    // Definition of Class Tube: (Length, topradius, botradius, *LeftDaughter, *RightDaughter, rmin, points, init, K, f1, f2, f3, R1, R2,  CT, LT);

    
    // =========================NETWORK MITCHELS CONSENSUS=================================================================================
    
    
//    Arteries[0] = new Tube( 0.410, 0.047, 0.047*(1-tap), 0, 0, rm, 40, 1, 0,f1,f2,f3, R1, R2, CT, 0);
    Arteries[0] = new Tube( Len, rtop, rbot, 0, 0, rm, 4, 1, 0,f1,f2,f3, R1, R2, CT, L);

    
    
      // In the next three statements the simulations are performed until
  // tstart = tend. That is this is without making any output during this
  // first period of time. If one needs output during this period, these three
  // lines should be commented out, and the entire simulation done within the
  // forthcomming while-loop.

  // Solves the equations until time equals tend.
  solver (Arteries, tstart, tend, k);
  tstart = tend;
  tend = tend + Deltat;

  // fprintf (stdout,"saves Q0\n");
  // Arteries[ 0] -> printQ0(fq0);

//  fprintf (stdout,"plots start\n");

  // The loop is continued until the final time
  // is reached. If one wants to make a plot of
  // the solution versus x, tend is set to final-
  // time in the above declaration.
  while (tend <= finaltime)
  {
    for (int j=0; j<nbrves; j++)
    {
      int ArtjN = Arteries[j]->N;
      for (int i=0; i<ArtjN; i++)
      {
        Arteries[j]->Qprv[i+1] = Arteries[j]->Qnew[i+1];
        Arteries[j]->Aprv[i+1] = Arteries[j]->Anew[i+1];
      }
    }

    // Solves the equations until time equals tend.
    solver (Arteries, tstart, tend, k);
//    fprintf (stdout,".");

    // A 2D plot of P(x_fixed,t) is made. The resulting 2D graph is due to
    // the while loop, since for each call of printPt only one point is set.
      Arteries[ 0] -> printPxt (fpu1, tend, 0);

    // The time within each print is increased.
    tstart = tend;
    tend   = tend + Deltat; // The current ending time is increased by Deltat.
  }
//  fprintf(stdout,"\n");

  // The following statements is used when timing the simulation.
//  fprintf(stdout,"nbrves = %d, Lax, ", nbrves);
//  clock_t c2 = clock(); // FIXME clock() may wrap after about 72 min.
//  int tsec = (int) ((double) (c2-c1)/CLOCKS_PER_SEC + 0.5);
//  fprintf(stdout,"cpu-time %d:%02d\n", tsec / 60, tsec % 60);
//  fprintf(stdout,"\n");

  // In order to termate the program correctly the vessel network and hence
  // all the vessels and their workspace are deleted.
  for (int i=0; i<nbrves; i++) delete Arteries[i];

  // Matrices and arrays are deleted
  for (int i=0; i<18; i++) delete[] fjac[i];

    
  return 0;
}
