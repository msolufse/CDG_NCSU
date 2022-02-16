/* The sor06.C main program */

// $Id: sor06.C,v 1.13 2010-10-20 15:38:28 mette Exp $

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
    int num_pts = 40;
    // Declare parameters: MJC 9/2019
    double f1, f2, f3, fs1, fs2, fs3;
    double alpha,beta, lrr1,lrr2, rm;
    int total_vessels, total_terminal, total_conn;
    
    if (argc < 14) {
        fprintf(stdout,"Not enough entries: Only %d inputs. Exiting now.\n",argc);
        return 1;
    }
    ////// Load in fluids parameters
    f1    = atof(argv[1]);
    f2    = atof(argv[2]);
    f3    = atof(argv[3]);
    fs1   = atof(argv[4]);
    fs2   = atof(argv[5]);
    fs3   = atof(argv[6]);
    alpha = atof(argv[7]);
    beta  = atof(argv[8]);
    lrr1   = atof(argv[9]);
    lrr2   = atof(argv[10]);
    rm    = atof(argv[11]);
    int cycles = 1; // NOTE: Can make this a parameter to pass in

    /////// Load in network parameters
    total_vessels  = atoi(argv[12]);
    total_terminal = atoi(argv[13]);
    total_conn     = total_vessels-total_terminal;
    nbrves         = total_vessels;
    
    
//     printf("PARAMETERS IN SOR06.C:\n f1: %lf\n f2: %lf\n f3: %lf\n fs1: %lf\n fs2: %lf\n fs3: %lf\n alpha: %lf\n beta: %lf\n lrr: %lf\n rm: %lf\n tot_ves %d\n tot_term %d\n tot_conn %d\n",f1,f2,f3,fs1,fs2,fs3,alpha,beta,lrr,rm, total_vessels,total_terminal,total_conn);
    
    
//    const char *namepu1  = "pu1.2d";
//    FILE *fpu1 = fopen (namepu1, "w");
//
//    const char *namepu2  = "pu2.2d";
//    FILE *fpu2 = fopen (namepu2, "w");//
//
//    const char *namepu3  = "pu3.2d";
//    FILE *fpu3 = fopen (namepu3, "w");
    
    // Workspace used by bound_right
    
    // Workspace used by bound_bif
    for(int i=0; i<18; i++) fjac[i] = new double[18];
    
//  clock_t c1 = clock();        // Only used when timing the program.
    tstart    = 0.0;            // Starting time.
    
    // The number of vessels in the network is given when the governing array of
    // vessels is declared.
    
    impedance_init_driver_(&tmstps);
    Tube   *Arteries[nbrves];                    // Array of blood vessels.
    
    
    
    
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
////////////////// Load in network ////////////////////////////////////////
    int conn_rows = (total_vessels-1)/2;
    int conn_cols = 3;
    int connectivity_matrix[conn_rows][conn_cols];
    int terminal_vessels[total_terminal];
    FILE *conn;
    conn = fopen("connectivity.txt","rt");
    int parent, daughter1, daughter2, r_in;
    int conn_id = 0;
    
    // Check to see if we have the connecitivity file
    if (conn == NULL)
    {
        fprintf(stdout,"Error: Connectivity File Does Not Exist \n");
        return 1;
    }
    while ((r_in = fscanf(conn, "%d %d %d", &parent, &daughter1, &daughter2)) != EOF)
    {
        connectivity_matrix[conn_id][0] = parent;
        connectivity_matrix[conn_id][1] = daughter1;
        connectivity_matrix[conn_id][2] = daughter2;
        conn_id++;
    }
    fclose(conn);
    
    
    /////////////// Load in terminal vessels
    
    FILE *terminal_file;
    terminal_file = fopen("terminal_vessels.txt","rt");
    for (int i=0; i<total_terminal; i++){
        fscanf(terminal_file, "%d", &terminal_vessels[i]);
    }
    fclose(terminal_file);
    
    
    /////////////// Load in Dimensions
    FILE *dim_file;
    dim_file = fopen("Dimensions.txt","rt");
    double dimensions_matrix[total_vessels][2];
    for (int dim_ID = 0; dim_ID < total_vessels; dim_ID++){
        fscanf(dim_file, "%5lf %5lf", &dimensions_matrix[dim_ID][0],&dimensions_matrix[dim_ID][1]);
    }
    fclose(dim_file);
    
    
    /* Initialization of the Arteries.
    // Definition of Class Tube:
//     Tube (double Length,
//               double topradius, double botradius,
//               Tube *LeftDaughter, Tube *RightDaughter,
//               double rmin_in, double points, int init, double K,
//               double f1, double f2, double f3, double fs1, double fs2, double fs3, double alpha_in, double beta_in, double lrr_in, double terminal_ID):
     Note: In this version, we pass the terminal ID to save impdance files from the ST*/
    /*                 FORWARD VERSION             */
     int curr_d1, curr_d2;
        int term_id = total_terminal-1;
        int bc_id = total_terminal-1;
        conn_id = total_conn-1;
        curr_d1 = 0;
        curr_d2 = 0;
        
    ///// A hard code for a single vessel
        if (total_vessels == 1){
            Arteries[0] = new Tube( dimensions_matrix[0][0], dimensions_matrix[0][1], dimensions_matrix[0][1], 0, 0, rm, num_pts, 1,0,f1,f2,f3,fs1,fs2,fs3, alpha, beta, lrr1,lrr2,0);
        }
        else if (total_vessels > 1){
            for (int i=total_vessels-1; i>=0; i--) {
                //// Define the connecitivity between the vessels
                if (i == connectivity_matrix[conn_id][0])
                {
                    curr_d1     = connectivity_matrix[conn_id][1];
                    curr_d2     = connectivity_matrix[conn_id][2];
                    conn_id--;
                }
                //first vessels
                if (i==0) {
    //                 fprintf(stdout,"At first vessel with dimensions %lf %lf  and daughters %d  %d\n",dimensions_matrix[i][0],dimensions_matrix[i][1], curr_d1, curr_d2);
                    Arteries[i] = new Tube( dimensions_matrix[0][0], dimensions_matrix[0][1], dimensions_matrix[0][1], Arteries[ curr_d1], Arteries[ curr_d2], rm, num_pts, 1,0, f1,f2,f3,fs1,fs2,fs3, alpha, beta, lrr1,lrr2,0);
                }
                else{
                    if (i == terminal_vessels[term_id] && term_id>=0){
                        Arteries[i] = new Tube( dimensions_matrix[i][0], dimensions_matrix[i][1], dimensions_matrix[i][1], 0, 0, rm, num_pts, 0, 0,f1,f2,f3,fs1,fs2,fs3, alpha, beta, lrr1,lrr2,term_id);
                        term_id--;
                        bc_id--;
                    }
                    else
                    {
    //                     fprintf(stdout,"At vessel %d with dimensions %lf %lf  and daughters %d  %d\n",i,dimensions_matrix[i][0],    dimensions_matrix[i][1], curr_d1, curr_d2);
                        Arteries[i] = new Tube( dimensions_matrix[i][0], dimensions_matrix[i][1], dimensions_matrix[i][1], Arteries[ curr_d1], Arteries[ curr_d2], rm, num_pts, 0, 0, f1,f2,f3,0,0,0,0,0,0,0,0);
                    }
                }
                
            }
            
        }
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    
    
    
    // Solves the equations until time equals tend.///////////////////////////
    /* ADDED BY Mitchel J. Colebank
     * Rather than specifying the number of cycles as an input to the function,
     * we want to test to see if the solution has converged. If so, we should exit.*/
    
    int period_counter = 1; // Count the number of periods you have solved for
    double norm_sol = 1e+6;
    double sol_tol  = 1e+1;
    printf("NORM_SOL: %f\n",norm_sol);
    double sol_p1[tmstps],sol_p2[tmstps];
    tend      = Deltat;
    
    
    // SOLVE THE MODEL ONCE
    // Note: Only want to test the pressure at the inlet
    int sol_ID = 0;
    while (tend<=period_counter*Period)
    {
    solver (Arteries, tstart, tend, k); 
    sol_p1[sol_ID] = Arteries[0]->P(0,Arteries[0]->Anew[0]); // for printing
    sol_p1[sol_ID] *= rho*g*Lr/conv;
    tstart = tend;
    tend   = tend + Deltat; // The current ending time is increased by Deltat.
    sol_ID++;
    }
    
    
    // LOOP FOR CONVERGENCE
    double sse;
    int max_count = 30;
    while (norm_sol>=sol_tol && period_counter<=max_count)
    {
        sol_ID = 0;
        sse    = 0;
        period_counter++;
        while (tend<=period_counter*Period)
        {
            solver (Arteries, tstart, tend, k);
            sol_p2[sol_ID] = Arteries[0]->P(0,Arteries[0]->Anew[0]); // for printing
            sol_p2[sol_ID] *= rho*g*Lr/conv;
            sse = sse+ sq(sol_p1[sol_ID]-sol_p2[sol_ID]);
            tstart = tend;
            tend   = tend + Deltat; // The current ending time is increased by Deltat.
            sol_ID++;
        }
        norm_sol = sse;
        memcpy (sol_p1, sol_p2, sizeof(sol_p2));
        printf("NORM_SOL:%f\n",norm_sol);
    }
    printf("num_cylces:%d\n",period_counter);
    
//   tend = tend + Deltat;

    //////////////////////////////////////////////
  // The loop is continued until the final time
  // is reached. If one wants to make a plot of
  // the solution versus x, tend is set to final-
  // time in the above declaration.
    
  period_counter++;
//   finaltime = period_counter*Period;
    
    const char *namepuALL  = "pu_ALL.2d";
    FILE *fpALL = fopen (namepuALL, "w");
    
  finaltime = (period_counter+(cycles-1))*Period;
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
	
	
//          Arteries[ 0] -> printPxt (fpu1, tend, 0);
//        Arteries[ connectivity_matrix[0][1]] -> printPxt (fpu2, tend, 0);
//        Arteries[ connectivity_matrix[0][2]] -> printPxt (fpu3, tend, 0);
      
         for (int save_id=0; save_id<nbrves; save_id++)
         {
             Arteries[ save_id] -> printAllt (fpALL, tend, 0);
         }
        for (int save_id=0; save_id<nbrves; save_id++)
        {
            Arteries[ save_id] -> printAllt (fpALL, tend, Arteries[save_id]->N);
        }
        // The time within each print is increased.
        tstart = tend;
        tend   = tend + Deltat; // The current ending time is increased by Deltat.
  }
  
    
    
    // In order to termate the program correctly the vessel network and hence
    // all the vessels and their workspace are deleted.
    for (int i=0; i<nbrves; i++) delete Arteries[i];
    
    // Matrices and arrays are deleted
    for (int i=0; i<18; i++) delete[] fjac[i];
         
         fclose(fpALL);

    
//    const char *namepu1  = "pu1.2d";
//    FILE *fpu1 = fopen (namepu1, "w");
//
//    const char *namepu2  = "pu2.2d";
//    FILE *fpu2 = fopen (namepu2, "w");//
//
//    const char *namepu3  = "pu3.2d";
//    FILE *fpu3 = fopen (namepu3, "w");
    
    
//    fclose(fpu1);
//    fclose(fpu2);
//    fclose(fpu3);
    
    return 0;
}
