/***************************************************************************/
/*                                                                         */
/* The junction.C main program                                             */
/*  Version: 1.0                                                           */
/*  Date: 24 March 2020                                                    */
/*                                                                         */
/*  Primary Authors: M.J. Colebank & J. Mackenzie                          */
/*  Key Contributers: M.U. Qureshi & M.S. Olufsen                          */
/*  Institute: North Carolina State University, Raleigh, NC, USA           */
/*  Funded by: NSF-DMS # 1615820 & AHA #19PRE34380459                      */
/*                                                                         */
/*  In contrast to previous versions, this code calls a generalized        */
/*  junction condition, which can be a single vessel, bifurcation, or n-   */
/*  furcation, and consructs the appropriate Jacobian.                     */
/***************************************************************************/


#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <algorithm>

#include "tools.h"
#include "arteries.h"
#include "junction.h"

using namespace std;

extern int max_D;


void junction (double theta, double gamma, Tube *Arteries[], int parent, int WM, double t, double Period)
{
    int D1,D2,D3;//,D4;
    D1 = Arteries[parent]->daughters[max_D*parent + 1];
    D2 = Arteries[parent]->daughters[max_D*parent + 2];
    D3 = Arteries[parent]->daughters[max_D*parent + 3];
    

    if (D2==0 ) {
        bound_screen_lnsearch(theta, gamma, Arteries, parent,D1,WM);
    }
    else if (D3==0){
//        fprintf(stdout,"Call bif with D1:%d D2:%d D3:%d\n",D1,D2,D3);
        bound_bif(theta, gamma, Arteries, parent,D1,D2,WM);
    }
    else if (D3==-1){
            bound_sten_bif(theta, gamma, Arteries, parent,D1,D2,WM);
        }

    else
    {
//        fprintf(stdout,"Call trif with D1:%d D2:%d D3:%d\n",D1,D2,D3);
        bound_trif(theta, gamma, Arteries, parent,D1,D2,D3,WM);
    }
    // Can add quadfurcation here if necessary.
}

// FROM GITLAB
/////////////////////////////////////////////////////////////////
///* ADDED BY MJC: This uses a newton-line search method rather than just Newton-Raphson
// */
void bound_screen_lnsearch(double theta, double gamma,Tube *Arteries[], int parent, int D1, int WM)
{
    // MJC: Try to define the daughters here
    Tube* PV = Arteries[ parent]; // Parent vessel
    Tube* B1 = Arteries[ D1];      // Branch 1
    int N = PV->N;
    int j = 1;
    int ok = false;
    const int ntrial = 1500000;
    

    double k1[2];
    double k2[2];
    double k3[2];
    double k4[2];

        // Define the screen parameters
        double k_darcy,k_forch,Ls;//,reynolds;
        // For now, just set values based on literature
        // MJC: for flow through a web-like stenoses, set Kt=0
        k_darcy = PV->sten_factor;//0.9; // Darcy-parameter/permeability
        k_forch = 0.10/sqrt(k_darcy);//0.05;                  // Forchhemier parameter (inertial permeability)
        Ls = PV->sten_length; //Length of the screen




    k1[0] = PV->Qold[N]  + theta*(PV->R2h[N-1]) + gamma*(PV->S2h[N-1]);
    k1[1] = B1->Qold[0] - theta*(B1->R2h[0]) + gamma*(B1->S2h[0]);

    k2[0] = PV->Aold[N] + theta*(PV->R1h[N-1]);
    k2[1] = B1->Aold[0] - theta*(B1->R1h[0]);

    k3[0]   = PV->Qh[N-1]/2.0;
    k3[1]   = B1->Qh[0]/2.0;

    k4[0]   = PV->Ah[N-1]/2.0;
    k4[1]   = B1->Ah[0]/2.0;

    double xb[12];

    // The approximative initial guesses are applied.
    xb[ 0] =  PV->Qh[N-1];                      //Initial guess for Q1_xb n+1
    xb[ 1] = (PV->Qold[N-1] + PV->Qold[N])/2.0;       //Initial guess for Q1_xb^n+0.5
    xb[ 2] =  PV->Qold[N];                      //Initial guess for Q1_xb+0.5 n+0.5
    xb[ 3] =  B1->Qh[0];                    //Initial guess for Q2_xb n+1
    xb[ 4] = (B1->Qold[0] + B1->Qold[1])/2.0; //Initial guess for Q2_xb n+0.5
    xb[ 5] =  B1->Qold[0];                  //Initial guess for Q2_xb+0.5 n+0.5
    xb[ 6] =  PV->Ah[N-1];                      //Initial guess for A1_xb n+1
    xb[ 7] = (PV->Aold[N-1] + PV->Aold[N])/2.0;       //Initial guess for A1_xb^n+0.5
    xb[ 8] =  PV->Aold[N];                      //Initial guess for A1_xb+0.5 n+0.5
    xb[ 9] =  B1->Ah[0];                    //Initial guess for A2_xb n+1
    xb[10] = (B1->Aold[0] + B1->Aold[1])/2.0; //Initial guess for A2_xb n+0.5
    xb[11] =  B1->Aold[0];                  //Initial guess for A2_xb+0.5 n+0.5


    // ADD NEWTON PARAMETERS w/ LINE SEARCH HERE
    double tolx = 1e-8;
    double tol_res = 1e-5;
    double tol_min = 1e-12;
    double stpmaxMAX = 100.0;
    double alam, alamin, a, b, disc, f2, rhs1, rhs2, stpmax, slope, sum, temp, test, tmplam;
    double alam2 = 0.0;
    double fvec[12], pln[12], grad[12], xold[12], fold[12];
    int indx[12];
    double d;
    double ALF = 1.0e-4;
    double SSE_old,SSE_new;
    
    
    sum = 0.0;
    stpmax = 0.0;
    // Calulcate stpmax for linesearch
    for (int i = 0; i<12; i++)
        sum+=xb[i]*xb[i];
    stpmax = stpmaxMAX*fmax(sqrt(sum),double(12));
    
//    fprintf(stdout,"stpmax: %lf\n",stpmax);
    
    
    double momentum = 0.9;//0.9;
    double x_mom[12];
    for (int i = 0; i<12; i++) fold[i] = 0.0; // Initialize
    // The residuals (fvec), and the Jacobian is determined, and if possible
    // the system of equations is solved.
    while (j <= ntrial && ok==false) // Find the zero
    {

         // Residual vector

        // Characteristic Q residual at n+1
        fvec[0]  = k1[0]  - xb[0] -
        theta*(sq(xb[2])/xb[8] + PV->Bh(N,xb[8],WM)) +
        gamma*(F(xb[2],xb[8])+PV->dBdx1h(N,xb[8],WM));

        fvec[1]  = k1[1]  - xb[3] +
        theta*(sq(xb[5])/xb[11] + B1->Bh(-1,xb[11],WM)) +
        gamma*(F(xb[5],xb[11])  + B1->dBdx1h(-1,xb[11],WM));

        // Characteristic A residual at n+1
        fvec[2]  = - theta*xb[2] - xb[6]  + k2[0];
        fvec[3]  =   theta*xb[5] - xb[9]  + k2[1];

        // Flow residuals at n+1/2 (ghost points)
        fvec[4]  = - xb[ 1] + xb[ 2]/2.0  + k3[0];
        fvec[5]  = - xb[ 4] + xb[ 5]/2.0  + k3[1];

        // Area residuals at n+1/2 (ghost points)
        fvec[6]  = - xb[ 7] + xb[ 8]/2.0  + k4[0];
        fvec[7]  = - xb[10] + xb[11]/2.0  + k4[1];

        // Flow conservation residuals (n+1/2 and n+1)
        fvec[8]  = - xb[ 1] + xb[ 4];
        fvec[9]  = - xb[ 0] + xb[ 3];

        // Define the pressure loss terms.
        double u1 = xb[1]/xb[7];
        double u2 = xb[0]/xb[6];
         double delta_PA =  Ls*(mu*u1/k_darcy + rho*sq(u1)*k_forch);//square
         double delta_PB =  Ls*(mu*u2/k_darcy + rho*sq(u2)*k_forch); //square
//        double delta_PA =  Ls*(mu*u1/k_darcy + rho*fabs(u1)*u1*k_forch); //abs
//        double delta_PB =  Ls*(mu*u2/k_darcy + rho*fabs(u2)*u2*k_forch); //abs



            fvec[10] = - PV->P(N,xb[7],WM) + B1->P(0,xb[10],WM) +  delta_PA;
            fvec[11] = - PV->P(N,xb[6],WM) + B1->P(0,xb[9],WM)  +  delta_PB;

//        fprintf(stdout,"ves: %d u1:%lf u2:%lf\n",parent,u1,u2);
//        fprintf(stdout,"Qold: %lf %lf Qhalf: %lf %lf\n",xb[0],xb[3],xb[1],xb[4]);
//        fprintf(stdout,"Aold: %lf %lf Ahalf: %lf %lf\n",xb[6],xb[9],xb[7],xb[10]);
//        fprintf(stdout,"square: dp:  %lf %lf\n",delta_PA,delta_PB);
//                fprintf(stdout,"abs: dp:  %lf %lf\n",Ls*(mu*u1/k_darcy + rho*fabs(u1)*u1*k_forch),Ls*(mu*u2/k_darcy + rho*fabs(u2)*u2*k_forch));



        double chi[10];

        // Here are the residuals for the characteristic matching for flow
        chi[0] = -2.0*theta*xb[ 2]/xb[8] + gamma*dFdQ(xb[8]);
        chi[2] =  2.0*theta*xb[ 5]/xb[11] + gamma*dFdQ(xb[11]);

        // Here are the residuals for the area characteristic matching
        chi[1] = theta*(sq(xb[2]/xb[8]) - PV->dBdAh(N,xb[8],WM)) +
                   gamma*(dFdA(xb[2],xb[8]) + PV->d2BdAdxh(N,xb[8],WM));

        chi[3] = theta*( -sq(xb[5]/xb[11]) + B1->dBdAh(-1,xb[11],WM)) +
                  gamma*(dFdA(xb[5],xb[11]) + B1->d2BdAdxh(-1,xb[11],WM));

        // Here is pressure conservation (n+1/2)
         chi[4]  = -PV->dPdA(N,xb[7],WM) - Ls*(mu*u1/(xb[7]*k_darcy) + 2.0*rho*sq(u1)*k_forch/xb[7]); //square
//        chi[4]  = -PV->dPdA(N,xb[7],WM) - Ls*(mu*xb[1]/(sq(xb[7])*k_darcy) +
//                    2.0*rho*k_forch*(fabs(xb[7]/xb[1])*cu(xb[7])/(sq(xb[1])*sq(xb[1])))); //abs

        chi[5]  =  B1->dPdA(0,xb[10],WM);

        // Here is pressure conservation (n+1)
         chi[6] = -PV->dPdA(N,xb[6],WM) - Ls*(mu*u2/(xb[6]*k_darcy) + 2.0*rho*sq(u2)*k_forch/xb[6]); //square
//        chi[6]  = -PV->dPdA(N,xb[6],WM) - Ls*(mu*xb[0]/(sq(xb[6])*k_darcy) +
//                    2.0*rho*k_forch*(fabs(xb[6]/xb[0])*cu(xb[6])/(sq(xb[0])*sq(xb[0])))); //abs

        chi[7] =  B1->dPdA(0,xb[9],WM);

        // Additional terms that arrive because of the pressure loss (depends on parent flow)
         chi[8] = Ls*(mu/(xb[7]*k_darcy) + 2.0*rho*u1*k_forch/xb[7]); //square
         chi[9] = Ls*(mu/(xb[6]*k_darcy) + 2.0*rho*u2*k_forch/xb[6]); //square
//        chi[8] = Ls*(mu/(xb[7]*k_darcy) + 2.0*rho*k_forch*(fabs(xb[7])*abs(xb[1])/cu(xb[7]))); //abs
//        chi[9] = Ls*(mu/(xb[6]*k_darcy) + 2.0*rho*k_forch*(fabs(xb[6])*abs(xb[0])/cu(xb[0]))); //abs


//        fprintf(stdout,"square chi6 term %lf abs chi 6 term %lf\n",
//               sq(xb[0])*cu(xb[6]),(abs(xb[6]/xb[0])*cu(xb[6])/(sq(xb[0])*sq(xb[0]))));
//        fprintf(stdout,"Ls:%lf, Kdarcy: %lf, Kforch: %lf rho: %lf mu: %lf\n",Ls,k_darcy,k_forch,rho,mu);
//        fprintf(stdout,"xb0:%lf, xb1: %lf, xb6: %lf xb6: %lf\n",xb[0],xb[1],xb[6],xb[7]);
//        fprintf(stdout,"1-over xb0:%lf, xb1: %lf, xb6: %lf xb6: %lf\n",1.0/xb[0],1.0/xb[1],1.0/xb[6],1.0/xb[7]);

        for (int row = 0; row < 12; row++)
            for (int col = 0; col < 12; col++)
                fjac[row][col] = 0.0;

        // The Jacobian.
        // MJC: Ordered by column
        // Column 0:
        fjac[ 0][ 0] = -1.0;
        fjac[ 0][ 2] = chi[0];
        fjac[ 0][ 8] = chi[1];


        fjac[ 1][ 3] = -1.0;
        fjac[ 1][ 5] = chi[2];
        fjac[ 1][11] = chi[3];


        fjac[ 2][2] = -theta;
        fjac[ 2][6] = -1.0;

        fjac[ 3][5] =  theta;
        fjac[ 3][9] = -1.0;


        fjac[ 4][1] = -1.0;
        fjac[ 4][2] =  0.5;


        fjac[ 5][4] = -1.0;
        fjac[ 5][5] =  0.5;


        fjac[ 6][ 7] = -1.0;
        fjac[ 6][ 8] = 0.5;


        fjac[7][10] = -1.0;
        fjac[7][11] = 0.5;


        fjac[8][1] = -1.0;
        fjac[8][4] =  1.0;


        fjac[9][0]  = -1.0;
        fjac[ 9][3] =  1.0;

        fjac[10][ 1] = chi[8]; //Loss
        fjac[10][ 7] = chi[4];
        fjac[10][10] = chi[5];

        fjac[11][ 0] = chi[9]; //Loss
        fjac[11][ 6] = chi[6]; // Might need an if statement for negative flow here
        fjac[11][ 9] = chi[7]; // Might need an if statement for negative flow here

        ///////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////
    // Line-search procedures
        
        // First, construct a gradient and old x/f values
        for (int i=0; i<12; i++)
        {
            sum = 0.0;
            for (int k=0; k<12; k++)
                sum += fjac[k][i]*fvec[k];
            grad[i] = sum;
            pln[i] = -fvec[i]; // Assign Newton direction
            x_mom[i] = xold[i]; // k-1 iterate
            xold[i] = xb[i];   // k iterate
            fold[i] = fvec[i];
//            fprintf(stdout,"grad: %lf p: %lf xold: %lf fold: %lf\n",grad[i],pln[i], xold[i],fold[i]);
        }
        
        
        
        ludcmp (fjac, 12, indx, &d);
        lubksb (fjac, 12, indx, pln);

        sum = 0.0;
        for (int i=0; i<12; i++)
        {
            sum+=pln[i]*pln[i];
        }
        // Scale step if too large
        if (sum>stpmax) {
            for (int i=0; i<12; i++)
            {
                pln[i]*=stpmax/sum;
            }
        }
        slope = 0.0;
        // Construct the slope for the gradient
        for (int i=0; i<12 ;i++)
        {
            slope += grad[i]*pln[i];
//            fprintf(stdout,"g: %lf p: %lf slope: %lf\n",grad[i],pln[i], slope);
        }
        if (slope >= 0.0)
        {
            fprintf(stdout,"Roundoff problem in lnsrch.\n");
            exit(1);
        }
        
        
        test=0.0;
        for (int i=0;i< 12; i++) {
            temp=abs(pln[i])/fmax(abs(xold[i]),1.0);
            if (temp > test) test=temp;
        }
        alamin=tolx/test;
        alam=1.0;
        f2 = 0.0;
        
        // Construct initial SSE
        SSE_old = 0.0;
        for (int i=0; i<12; i++) SSE_old+=fvec[i]*fvec[i];
        SSE_old*=0.5;
        
        
        // Now iterate until the Armijo condition is satisfied

        bool armijo_flag = false;
        while (!armijo_flag) {
            // Take an initial full gradient step
            if (j>5)
            {
                for (int i=0; i<12; i++)
                {
                    xb[i] = xold[i] + alam*pln[i] + momentum*(xold[i]-x_mom[i]);
                    //                fprintf(stdout,"ves: %d   j:  %d i: %d  ALAM: %5.12lf fvec: %5.12lf xb: %5.12lf    xold: %5.12lf     pln: %5.12lf    mom:   %5.12lf  \n",parent,j,i,alam,fvec[i],xb[i],xold[i],pln[i],momentum*(xold[i]-x_mom[i]));
                }
            }
            else
            {
                for (int i=0; i<12; i++)
                    xb[i] = xold[i] + alam*pln[i];
            }
            // Now recompute the residual values
            fvec[0]  = k1[0]  - xb[0] -
            theta*(sq(xb[2])/xb[8] + PV->Bh(N,xb[8],WM)) +
            gamma*(F(xb[2],xb[8])+PV->dBdx1h(N,xb[8],WM));

            fvec[1]  = k1[1]  - xb[3] +
            theta*(sq(xb[5])/xb[11] + B1->Bh(-1,xb[11],WM)) +
            gamma*(F(xb[5],xb[11])  + B1->dBdx1h(-1,xb[11],WM));

            // Characteristic A residual at n+1
            fvec[2]  = - theta*xb[2] - xb[6]  + k2[0];
            fvec[3]  =   theta*xb[5] - xb[9]  + k2[1];

            // Flow residuals at n+1/2 (ghost points)
            fvec[4]  = - xb[ 1] + xb[ 2]/2.0  + k3[0];
            fvec[5]  = - xb[ 4] + xb[ 5]/2.0  + k3[1];

            // Area residuals at n+1/2 (ghost points)
            fvec[6]  = - xb[ 7] + xb[ 8]/2.0  + k4[0];
            fvec[7]  = - xb[10] + xb[11]/2.0  + k4[1];

            // Flow conservation residuals (n+1/2 and n+1)
            fvec[8]  = - xb[ 1] + xb[ 4];
            fvec[9]  = - xb[ 0] + xb[ 3];

            // Define the pressure loss terms.
            double u1 = xb[1]/xb[7];
            double u2 = xb[0]/xb[6];
             double delta_PA =  Ls*(mu*u1/k_darcy + rho*sq(u1)*k_forch);//square
             double delta_PB =  Ls*(mu*u2/k_darcy + rho*sq(u2)*k_forch); //square
    //        double delta_PA =  Ls*(mu*u1/k_darcy + rho*fabs(u1)*u1*k_forch); //abs
    //        double delta_PB =  Ls*(mu*u2/k_darcy + rho*fabs(u2)*u2*k_forch); //abs
                fvec[10] = - PV->P(N,xb[7],WM) + B1->P(0,xb[10],WM) +  delta_PA;
                fvec[11] = - PV->P(N,xb[6],WM) + B1->P(0,xb[9],WM)  +  delta_PB;
            
            
            // Calculate SSE
            SSE_new = 0.0;
            for (int i=0; i<12; i++) SSE_new+=fvec[i]*fvec[i];
            SSE_new*=0.5;
            if (alam<alamin) {
                for (int i=0; i<12; i++) xb[i] = xold[i];
                armijo_flag = true;
//                break;
            }
            else if (SSE_new <= SSE_old +ALF*alam*slope)
            {
                armijo_flag=true;
//                fprintf(stdout,"Armijo true: %lf %lf %lf\n",alam,SSE_new,SSE_old);
            }
            else
            {
                if (alam == 1.0)
                {
                    tmplam = -slope/(2.0*(SSE_new-SSE_old-slope));
                }
                else
                {
                    rhs1 = SSE_new-SSE_old-alam*slope;
                    rhs2 = f2-SSE_old-alam2*slope;
                    a    = (rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
                    b    = (-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
                    if (a == 0.0) tmplam = -slope/(2.0*b);
                    else {
                        disc=b*b-3.0*a*slope;
                        if (disc < 0.0) tmplam=0.5*alam;
                        else if (b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
                        else tmplam=-slope/(b+sqrt(disc));
                        }
                        if (tmplam>0.5*alam)
                            tmplam=0.5*alam;
                }
//                if (isnan(tmplam))
//                 {
//                     printf(" -  tmplam not a number\n");
//                     abort();
//                 }
                alam2=alam;
                f2 = SSE_new;
                alam=fmax(tmplam,0.1*alam);
            }
        
        }
        
        // Now that we have run through the linsearch, fvec and xb have been updated.
        // Now test for convergence and update
        
        
        // Test for residual
        test = 0.0;
        for (int i=0; i<12; i++) {
            if (fabs(fvec[i]) > test) {
                test = fabs(fvec[i]);
            }
        }
        ok = (test < tol_res);
        
        // Test for gradient going to zero
        if (ok==false)
                {
                    test=0.0;
//                    den = fmax(SSE_new,0.5*12.0);
                    for (int i=0; i<12; i++)
                    {
        //                temp = fabs(grad[i])*fmax(fabs(xb[i]),1.0)/den;
                        temp = grad[i]*grad[i]; // test for gradient convergence
                        if (temp > test) test = temp;
                    }
                    ok = (test<tol_min);
                    if (ok==true) fprintf(stdout,"Gradient passed: %12.12lf\n",test);
                }
        
        // Test for converged x value
        if (ok==false)
        {
        test = 0.0;
        for (int i=0; i<12; i++) {
            if (fabs(xb[i]-xold[i]) > test) {
                test = fabs(xb[i]-xold[i]);
            }
        }
        ok = fmax(ok, (test < tolx));
        }
        
        if (j>ntrial-5)
        {
            fprintf(stdout,"xb:");
            for (int i=0; i<12; i++)
                fprintf(stdout,"   %5.15lf    ",xb[i]);
            fprintf(stdout,"\n");
        }
        
        j = j+1;
        }
        
        ///////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////

//        fprintf(stdout,"AREA DIFFERENCE 2: %lf \n \n",xb[6]-xb[9]);



    // Solutions is applied, and right boundary is updated.
    PV->Anew[N] = xb[ 6];
    PV->Qnew[N] = xb[ 0];
    B1->Anew[0] = xb[ 9];
    B1->Qnew[0] = xb[ 3];
    
//    if (alam<1.0) fprintf(stdout,"Armijo: %lf\n",alam);
    
    
//   fprintf(stdout,"Delta P:%lf\n",(PV->P(PV->N,xb[6],WM)-B1->P(0,xb[9],WM))*rho*g*Lr/conv);
//   fprintf(stdout,"Delta A:%lf\n",(xb[9]-xb[6])/xb[6]);
//    fprintf(stdout,"Delta Q:%lf\n",(xb[3]-xb[0])/xb[0]);
    if (j >=ntrial){
        error ("arteries.C","Root not found in the screen");
        fprintf(stdout,"Vessel: %d\n",parent);
        for (int i=0; i<12; i++) fprintf(stdout,"  %5.5lf  ",fvec[i]);
        fprintf(stdout,"\n \n \n");
    }
}




// /////////////////////////////////////////////////////////////////
// ///* ADDED BY MJC: This uses a newton-line search method rather than just Newton-Raphson
// // */
// // FROM CODE THAT ISN"T WORKING
// void bound_screen_lnsearch(double theta, double gamma,Tube *Arteries[], int parent, int D1, int WM)
// {
//     // MJC: Try to define the daughters here
//     Tube* PV = Arteries[ parent]; // Parent vessel
//     Tube* B1 = Arteries[ D1];      // Branch 1
//     int N = PV->N;
//     int j = 1;
//     int ok = false;
//     const int ntrial = 2000000;
//     
// 
//     double k1[2];
//     double k2[2];
//     double k3[2];
//     double k4[2];
// 
//         // Define the screen parameters
//         double k_darcy,k_forch,Ls;//,reynolds;
//         // For now, just set values based on literature
//         // MJC: for flow through a web-like stenoses, set Kt=0
//         k_darcy = PV->sten_factor;//0.9; // Darcy-parameter/permeability
//         k_forch = 0.10/sqrt(k_darcy);//0.05;                  // Forchhemier parameter (inertial permeability)
//         Ls = PV->sten_length; //Length of the screen
// 
// 
// 
// 
//     k1[0] = PV->Qold[N]  + theta*(PV->R2h[N-1]) + gamma*(PV->S2h[N-1]);
//     k1[1] = B1->Qold[0] - theta*(B1->R2h[0]) + gamma*(B1->S2h[0]);
// 
//     k2[0] = PV->Aold[N] + theta*(PV->R1h[N-1]);
//     k2[1] = B1->Aold[0] - theta*(B1->R1h[0]);
// 
//     k3[0]   = PV->Qh[N-1]/2.0;
//     k3[1]   = B1->Qh[0]/2.0;
// 
//     k4[0]   = PV->Ah[N-1]/2.0;
//     k4[1]   = B1->Ah[0]/2.0;
// 
//     double xb[12];
// 
//     // The approximative initial guesses are applied.
//     xb[ 0] =  PV->Qh[N-1];                      //Initial guess for Q1_xb n+1
//     xb[ 1] = (PV->Qold[N-1] + PV->Qold[N])/2.0;       //Initial guess for Q1_xb^n+0.5
//     xb[ 2] =  PV->Qold[N];                      //Initial guess for Q1_xb+0.5 n+0.5
//     xb[ 3] =  B1->Qh[0];                    //Initial guess for Q2_xb n+1
//     xb[ 4] = (B1->Qold[0] + B1->Qold[1])/2.0; //Initial guess for Q2_xb n+0.5
//     xb[ 5] =  B1->Qold[0];                  //Initial guess for Q2_xb+0.5 n+0.5
//     xb[ 6] =  PV->Ah[N-1];                      //Initial guess for A1_xb n+1
//     xb[ 7] = (PV->Aold[N-1] + PV->Aold[N])/2.0;       //Initial guess for A1_xb^n+0.5
//     xb[ 8] =  PV->Aold[N];                      //Initial guess for A1_xb+0.5 n+0.5
//     xb[ 9] =  B1->Ah[0];                    //Initial guess for A2_xb n+1
//     xb[10] = (B1->Aold[0] + B1->Aold[1])/2.0; //Initial guess for A2_xb n+0.5
//     xb[11] =  B1->Aold[0];                  //Initial guess for A2_xb+0.5 n+0.5
// 
// 
//     // ADD NEWTON PARAMETERS w/ LINE SEARCH HERE
//     double tolx = 1e-6;    // Criteria for x
//     double tol_res = 1e-6; // Criteria for residual
//     double tol_min = 1e-6; // Criteria for gradient
//     double stpmaxMAX = 500.0;
//     double alam, alamin, a, b, disc, f2, rhs1, rhs2, stpmax, slope, sum, temp, test, tmplam;
//     double alam2 = 0.0;
//     double den;
//     double fvec[12], pln[12], grad[12], xold[12], fold[12];
//     int indx[12];
//     double d;
//     double ALF = 1.0e-4;
//     double SSE_old,SSE_new;
//     
//     
//     sum = 0.0;
//     stpmax = 0.0;
//     // Calulcate stpmax for linesearch
//     for (int i = 0; i<12; i++)
//         sum+=xb[i]*xb[i];
//     stpmax = stpmaxMAX*fmax(sqrt(sum),double(12));
//     
// //    fprintf(stdout,"stpmax: %lf\n",stpmax);
//     
//     
//     
//     for (int i = 0; i<12; i++) fold[i] = 0.0; // Initialize
//     
//     // MJC: Testing purposes:
//     double xinit[12];
//     for (int i = 0; i<12; i++) xinit[i] = xb[i];
//     // The residuals (fvec), and the Jacobian is determined, and if possible
//     // the system of equations is solved.
//     while (j <= ntrial && ok==false) // Find the zero
//     {
// 
//          // Residual vector
// 
//         // Characteristic Q residual at n+1
//         fvec[0]  = k1[0]  - xb[0] -
//         theta*(sq(xb[2])/xb[8] + PV->Bh(N,xb[8],WM)) +
//         gamma*(F(xb[2],xb[8])+PV->dBdx1h(N,xb[8],WM));
// 
//         fvec[1]  = k1[1]  - xb[3] +
//         theta*(sq(xb[5])/xb[11] + B1->Bh(-1,xb[11],WM)) +
//         gamma*(F(xb[5],xb[11])  + B1->dBdx1h(-1,xb[11],WM));
// 
//         // Characteristic A residual at n+1
//         fvec[2]  = - theta*xb[2] - xb[6]  + k2[0];
//         fvec[3]  =   theta*xb[5] - xb[9]  + k2[1];
// 
//         // Flow residuals at n+1/2 (ghost points)
//         fvec[4]  = - xb[ 1] + xb[ 2]/2.0  + k3[0];
//         fvec[5]  = - xb[ 4] + xb[ 5]/2.0  + k3[1];
// 
//         // Area residuals at n+1/2 (ghost points)
//         fvec[6]  = - xb[ 7] + xb[ 8]/2.0  + k4[0];
//         fvec[7]  = - xb[10] + xb[11]/2.0  + k4[1];
// 
//         // Flow conservation residuals (n+1/2 and n+1)
//         fvec[8]  = - xb[ 1] + xb[ 4];
//         fvec[9]  = - xb[ 0] + xb[ 3];
// 
//         // Define the pressure loss terms.
//         double u1 = fabs(xb[1]/xb[7]);
//         double u2 = fabs(xb[0]/xb[6]);
//         double backflow_flag;
//         if (xb[0]>=0.0)
//             backflow_flag = -1.0;
//         else
//             backflow_flag = 1.0;
//         
//          double delta_PA =  backflow_flag*Ls*(mu*u1/k_darcy + rho*sq(u1)*k_forch);//square
//          double delta_PB =  backflow_flag*Ls*(mu*u2/k_darcy + rho*sq(u2)*k_forch); //square
// 
// //        double delta_PA =  Ls*(mu*u1/k_darcy + rho*fabs(u1)*u1*k_forch); //abs
// //        double delta_PB =  Ls*(mu*u2/k_darcy + rho*fabs(u2)*u2*k_forch); //abs
// 
// 
// 
//             fvec[10] = - PV->P(N,xb[7],WM) + B1->P(0,xb[10],WM) +  delta_PA;
//             fvec[11] = - PV->P(N,xb[6],WM) + B1->P(0,xb[9],WM)  +  delta_PB;
// 
// //        fprintf(stdout,"ves: %d u1:%lf u2:%lf\n",parent,u1,u2);
// //        fprintf(stdout,"Qold: %lf %lf Qhalf: %lf %lf\n",xb[0],xb[3],xb[1],xb[4]);
// //        fprintf(stdout,"Aold: %lf %lf Ahalf: %lf %lf\n",xb[6],xb[9],xb[7],xb[10]);
// //        fprintf(stdout,"square: dp:  %lf %lf\n",delta_PA,delta_PB);
// //                fprintf(stdout,"abs: dp:  %lf %lf\n",Ls*(mu*u1/k_darcy + rho*fabs(u1)*u1*k_forch),Ls*(mu*u2/k_darcy + rho*fabs(u2)*u2*k_forch));
// 
// 
// 
//         double chi[10];
// 
//         // Here are the residuals for the characteristic matching for flow
//         chi[0] = -2.0*theta*xb[ 2]/xb[8] + gamma*dFdQ(xb[8]);
//         chi[2] =  2.0*theta*xb[ 5]/xb[11] + gamma*dFdQ(xb[11]);
// 
//         // Here are the residuals for the area characteristic matching
//         chi[1] = theta*(sq(xb[2]/xb[8]) - PV->dBdAh(N,xb[8],WM)) +
//                    gamma*(dFdA(xb[2],xb[8]) + PV->d2BdAdxh(N,xb[8],WM));
// 
//         chi[3] = theta*( -sq(xb[5]/xb[11]) + B1->dBdAh(-1,xb[11],WM)) +
//                   gamma*(dFdA(xb[5],xb[11]) + B1->d2BdAdxh(-1,xb[11],WM));
// 
//         // Here is pressure conservation (n+1/2)
//          chi[4]  = -PV->dPdA(N,xb[7],WM) - backflow_flag*Ls*(mu*u1/(xb[7]*k_darcy) + 2.0*rho*sq(u1)*k_forch/xb[7]); //square
// //        chi[4]  = -PV->dPdA(N,xb[7],WM) - Ls*((mu/k_darcy)*xb[1]/(sq(xb[7])) +
// //        rho*k_forch*(2.0*cu(xb[1])/(xb[7]*cu(xb[7])))*(fabs(xb[1])/fabs(xb[7]))); //abs V2
// 
//         chi[5]  =  B1->dPdA(0,xb[10],WM);
// 
//         // Here is pressure conservation (n+1)
//          chi[6] = -PV->dPdA(N,xb[6],WM) - backflow_flag*Ls*(mu*u2/(xb[6]*k_darcy) + 2.0*rho*sq(u2)*k_forch/xb[6]); //square
// //        chi[6]  = -PV->dPdA(N,xb[6],WM) - Ls*((mu/k_darcy)*xb[0]/(sq(xb[6])) +
// //        rho*k_forch*(2.0*cu(xb[0])/(xb[6]*cu(xb[6])))*(fabs(xb[0])/fabs(xb[6]))); //abs V2
// 
//         chi[7] =  B1->dPdA(0,xb[9],WM);
// 
//         // Additional terms that arrive because of the pressure loss (depends on parent flow)
//         // Note that sign is accounted for here, since we use xb
//          chi[8] = Ls*(mu/(xb[7]*k_darcy) + 2.0*rho*xb[1]*k_forch/sq(xb[7])); //square
//          chi[9] = Ls*(mu/(xb[6]*k_darcy) + 2.0*rho*xb[0]*k_forch/sq(xb[6])); //square
// //        chi[8] = Ls*( (mu/k_darcy)/xb[7] + rho*k_forch*2.0*(fabs(xb[1])*fabs(xb[7])/cu(xb[7]) )); //abs V2
// //        chi[9] = Ls*( (mu/k_darcy)/xb[6] + rho*k_forch*2.0*(fabs(xb[0])*fabs(xb[6])/cu(xb[6]) )); //abs V2
// 
// 
// //        fprintf(stdout,"square chi6 term %lf abs chi 6 term %lf\n",
// //               sq(xb[0])*cu(xb[6]),(abs(xb[6]/xb[0])*cu(xb[6])/(sq(xb[0])*sq(xb[0]))));
// //        fprintf(stdout,"Ls:%lf, Kdarcy: %lf, Kforch: %lf rho: %lf mu: %lf\n",Ls,k_darcy,k_forch,rho,mu);
// //        fprintf(stdout,"xb0:%lf, xb1: %lf, xb6: %lf xb6: %lf\n",xb[0],xb[1],xb[6],xb[7]);
// //        fprintf(stdout,"1-over xb0:%lf, xb1: %lf, xb6: %lf xb6: %lf\n",1.0/xb[0],1.0/xb[1],1.0/xb[6],1.0/xb[7]);
// 
//         for (int row = 0; row < 12; row++)
//             for (int col = 0; col < 12; col++)
//                 fjac[row][col] = 0.0;
// 
//         // The Jacobian.
//         // MJC: Ordered by column
//         // Column 0:
//         fjac[ 0][ 0] = -1.0;
//         fjac[ 0][ 2] = chi[0];
//         fjac[ 0][ 8] = chi[1];
// 
// 
//         fjac[ 1][ 3] = -1.0;
//         fjac[ 1][ 5] = chi[2];
//         fjac[ 1][11] = chi[3];
// 
// 
//         fjac[ 2][2] = -theta;
//         fjac[ 2][6] = -1.0;
// 
//         fjac[ 3][5] =  theta;
//         fjac[ 3][9] = -1.0;
// 
// 
//         fjac[ 4][1] = -1.0;
//         fjac[ 4][2] =  0.5;
// 
// 
//         fjac[ 5][4] = -1.0;
//         fjac[ 5][5] =  0.5;
// 
// 
//         fjac[ 6][ 7] = -1.0;
//         fjac[ 6][ 8] = 0.5;
// 
// 
//         fjac[7][10] = -1.0;
//         fjac[7][11] = 0.5;
// 
// 
//         fjac[8][1] = -1.0;
//         fjac[8][4] =  1.0;
// 
// 
//         fjac[9][0]  = -1.0;
//         fjac[ 9][3] =  1.0;
// 
//         fjac[10][ 1] = chi[8]; //Loss
//         fjac[10][ 7] = chi[4];
//         fjac[10][10] = chi[5];
// 
//         fjac[11][ 0] = chi[9]; //Loss
//         fjac[11][ 6] = chi[6]; // Might need an if statement for negative flow here
//         fjac[11][ 9] = chi[7]; // Might need an if statement for negative flow here
//         
// //        for (int i=0; i<12; i++) {
// //            if (fabs(xb[i])>5.0) {
// //                fprintf(stdout,"j: %d i: %d xb: %lf\n",j,i,xb[i]);
// //            }
// //            if (fabs(fvec[i])>5.0) {
// //                fprintf(stdout,"j: %d i: %d res: %lf\n",j,i,fvec[i]);
// //            }
// //        }
// 
//         ///////////////////////////////////////////////////////////////////////////////////
//         ///////////////////////////////////////////////////////////////////////////////////
//         ///////////////////////////////////////////////////////////////////////////////////
//     // Line-search procedures
//         
//         // First, construct a gradient and old x/f values
//         for (int i=0; i<12; i++)
//         {
//             sum = 0.0;
//             for (int k=0; k<12; k++)
//                 sum += fjac[k][i]*fvec[k];
//             grad[i] = sum;
//             pln[i] = -fvec[i]; // Assign Newton direction
//             xold[i] = xb[i];
//             fold[i] = fvec[i];
// //            fprintf(stdout,"grad: %lf p: %lf xold: %lf fold: %lf\n",grad[i],pln[i], xold[i],fold[i]);
//         }
//         
//         
//         
//         ludcmp (fjac, 12, indx, &d);
//         lubksb (fjac, 12, indx, pln);
// 
//         sum = 0.0;
//         for (int i=0; i<12; i++)
//         {
//             sum+=pln[i]*pln[i];
//         }
//         // Scale step if too large
//         if (sum>stpmax) {
// //            fprintf(stdout,"stpmax too large: %lf\n",sum);
//             for (int i=0; i<12; i++)
//             {
//                 pln[i]*=stpmax/sum;
// //                fprintf(stdout,"p %d: %5.5lf",i,pln[i]);
//             }
// //            fprintf(stdout,"\n");
//         }
//         slope = 0.0;
//         // Construct the slope for the gradient
//         for (int i=0; i<12 ;i++)
//         {
//             slope += grad[i]*pln[i];
// //            fprintf(stdout,"g: %lf p: %lf slope: %lf\n",grad[i],pln[i], slope);
//         }
//         if (slope >= 0.0)
//         {
//             fprintf(stdout,"Roundoff problem in lnsrch.\n");
//             exit(1);
//         }
//         
//         
//         test=0.0;
//         for (int i=0;i< 12; i++) {
//             temp=abs(pln[i])/fmax(abs(xold[i]),1.0);
//             if (temp > test) test=temp;
//         }
//         alamin=tolx/test;
//         alam=1.0;
//         f2 = 0.0;
//         
//         // Construct initial SSE
//         SSE_old = 0.0;
//         for (int i=0; i<12; i++) SSE_old+=fvec[i]*fvec[i];
//         SSE_old*=0.5;
//         
//         
//         // Now iterate until the Armijo condition is satisfied
//         bool armijo_flag = false;
//         while (!armijo_flag) {
//             // Take an initial full gradient step
//             for (int i=0; i<12; i++)
//             {
//                 xb[i] = xold[i] + alam*pln[i];
//             }
//             // Now recompute the residual values
//             fvec[0]  = k1[0]  - xb[0] -
//             theta*(sq(xb[2])/xb[8] + PV->Bh(N,xb[8],WM)) +
//             gamma*(F(xb[2],xb[8])+PV->dBdx1h(N,xb[8],WM));
// 
//             fvec[1]  = k1[1]  - xb[3] +
//             theta*(sq(xb[5])/xb[11] + B1->Bh(-1,xb[11],WM)) +
//             gamma*(F(xb[5],xb[11])  + B1->dBdx1h(-1,xb[11],WM));
// 
//             // Characteristic A residual at n+1
//             fvec[2]  = - theta*xb[2] - xb[6]  + k2[0];
//             fvec[3]  =   theta*xb[5] - xb[9]  + k2[1];
// 
//             // Flow residuals at n+1/2 (ghost points)
//             fvec[4]  = - xb[ 1] + xb[ 2]/2.0  + k3[0];
//             fvec[5]  = - xb[ 4] + xb[ 5]/2.0  + k3[1];
// 
//             // Area residuals at n+1/2 (ghost points)
//             fvec[6]  = - xb[ 7] + xb[ 8]/2.0  + k4[0];
//             fvec[7]  = - xb[10] + xb[11]/2.0  + k4[1];
// 
//             // Flow conservation residuals (n+1/2 and n+1)
//             fvec[8]  = - xb[ 1] + xb[ 4];
//             fvec[9]  = - xb[ 0] + xb[ 3];
// 
//             // Define the pressure loss terms.
//             double u1 = fabs(xb[1]/xb[7]);
//             double u2 = fabs(xb[0]/xb[6]);
//             
//             // Account for backflow
//             if (xb[0]>=0.0)
//                 backflow_flag = -1.0;
//             else
//                 backflow_flag = 1.0;
//             
//              double delta_PA =  backflow_flag*Ls*(mu*u1/k_darcy + rho*sq(u1)*k_forch);//square
//              double delta_PB =  backflow_flag*Ls*(mu*u2/k_darcy + rho*sq(u2)*k_forch); //square
//             
// //            double delta_PA =  Ls*(mu*u1/k_darcy + rho*fabs(u1)*u1*k_forch); //abs
// //            double delta_PB =  Ls*(mu*u2/k_darcy + rho*fabs(u2)*u2*k_forch); //abs
//             
//             fvec[10] = - PV->P(N,xb[7],WM) + B1->P(0,xb[10],WM) +  delta_PA;
//             fvec[11] = - PV->P(N,xb[6],WM) + B1->P(0,xb[9],WM)  +  delta_PB;
//             
// //            if (xb[0]>0.0)
// //            {
// //                fvec[10] = - PV->P(N,xb[7],WM) + B1->P(0,xb[10],WM) +  delta_PA;
// //                fvec[11] = - PV->P(N,xb[6],WM) + B1->P(0,xb[9],WM)  +  delta_PB;
// //            } else
// //            {
// //                fvec[10] = - PV->P(N,xb[7],WM) + B1->P(0,xb[10],WM) -  delta_PA;
// //                fvec[11] = - PV->P(N,xb[6],WM) + B1->P(0,xb[9],WM)  -  delta_PB;
// //            }
//             
//             // Calculate SSE
//             SSE_new = 0.0;
//             for (int i=0; i<12; i++) SSE_new+=fvec[i]*fvec[i];
//             SSE_new*=0.5;
//             if (alam<alamin) {
//                 for (int i=0; i<12; i++) xb[i] = xold[i];
//                 armijo_flag = true;
//             }
//             else if (SSE_new <= SSE_old +ALF*alam*slope)
//             {
//                 armijo_flag=true;
// //                fprintf(stdout,"Armijo true: %lf %lf %lf\n",alam,SSE_new,SSE_old);
//             }
//             else
//             {
//                 if (alam == 1.0)
//                 {
//                     tmplam = -slope/(2.0*(SSE_new-SSE_old-slope));
//                 }
//                 else
//                 {
//                     rhs1 = SSE_new-SSE_old-alam*slope;
//                     rhs2 = f2-SSE_old-alam2*slope;
//                     a    = (rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
//                     b    = (-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
//                     if (a == 0.0) tmplam = -slope/(2.0*b);
//                     else {
//                         disc=b*b-3.0*a*slope;
//                         if (disc < 0.0) tmplam=0.5*alam;
//                         else if (b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
//                         else tmplam=-slope/(b+sqrt(disc));
//                         }
//                         if (tmplam>0.5*alam)
//                         {
// //                            fprintf(stdout,"tmplam: %12.5lf,  tmplam_new: %12.5lf\n",tmplam,0.5*alam);
//                             tmplam=0.5*alam;
//                         }
//                 }
// //                if (isnan(tmplam))
// //                 {
// //                     printf(" -  tmplam not a number\n");
// //                     abort();
// //                 }
//                 alam2=alam;
//                 f2 = SSE_new;
//                 alam=fmax(tmplam,0.1*alam);
//             }
//         
//         }
//         
//         // Now that we have run through the linsearch, fvec and xb have been updated.
//         // Now test for convergence and update
//         
//         
//         // Test for residual
//         test = 0.0;
//         for (int i=0; i<12; i++) {
//             if (fabs(fvec[i]) > test) {
//                 test = fabs(fvec[i]);
//             }
//         }
//         ok = (test < tol_res);
//         
//         // Test for zero gradient (spurious convergence)
//         if (ok==false)
//         {
//             test=0.0;
//             den = fmax(SSE_new,0.5*12.0);
//             for (int i=0; i<12; i++)
//             {
//                 temp = fabs(grad[i])*fmax(fabs(xb[i]),1.0)/den;
//                 if (temp > test) test = temp;
//             }
//             ok = (test<tol_min);
//         }
//         
//         // Test for converged x value
//         if (ok==false)
//             {
//             test = 0.0;
//             for (int i=0; i<12; i++) {
//                 if (fabs(xb[i]-xold[i]) > test) {
//                     test = fabs(xb[i]-xold[i]);
//                 }
//             }
//             ok =  (test < tolx);
//         }
//         
//         j = j+1;
//         
//         if (j>ntrial-20) {
//             fprintf(stdout,"j: %d  alam: %12.5lf\n",j,alam);
//             for (int i=0; i<12; i++) fprintf(stdout,"  %12.5lf  ",xb[i]);
//             fprintf(stdout,"\n");
//             for (int i=0; i<12; i++) fprintf(stdout,"  %12.5lf  ",fvec[i]);
//             fprintf(stdout,"\n");
//             }
//         }
//         
//         ///////////////////////////////////////////////////////////////////////////////////
//         ///////////////////////////////////////////////////////////////////////////////////
//         ///////////////////////////////////////////////////////////////////////////////////
//         ///////////////////////////////////////////////////////////////////////////////////
// 
// //        fprintf(stdout,"AREA DIFFERENCE 2: %lf \n \n",xb[6]-xb[9]);
// 
// 
// 
//     // Solutions is applied, and right boundary is updated.
//     PV->Anew[N] = xb[ 6];
//     PV->Qnew[N] = xb[ 0];
//     B1->Anew[0] = xb[ 9];
//     B1->Qnew[0] = xb[ 3];
//     
// //    if (alam<1.0) fprintf(stdout,"Armijo: %lf\n",alam);
//     
//     
// //   fprintf(stdout,"Delta P:%lf\n",(PV->P(PV->N,xb[6],WM)-B1->P(0,xb[9],WM))*rho*g*Lr/conv);
// //   fprintf(stdout,"Delta A:%lf\n",(xb[9]-xb[6])/xb[6]);
// //    fprintf(stdout,"Delta Q:%lf\n",(xb[3]-xb[0])/xb[0]);
//     if (j >=ntrial){
//         error ("arteries.C","Root not found in the screen");
//         fprintf(stdout,"Armijo: %lf\n",alam);
//         fprintf(stdout,"Pressure drop: %5.5lf\n",(PV->P(N,xb[6],WM)-B1->P(0,xb[9],WM))/conv);
//         for (int i=0; i<12; i++) fprintf(stdout,"%5.5lf  %5.5lf  ",xinit[i], xb[i]);
//         fprintf(stdout,"\n \n \n");
//         for (int i=0; i<12; i++) fprintf(stdout,"  %5.5lf  ",fvec[i]);
//         fprintf(stdout,"\n \n \n");
//     }
// //    else if (j>1000)
// //    {
// //        fprintf(stdout,"j %d: fvec[10]: %5.5lf fvec[11]: %5.5lf\n",j,fvec[10],fvec[11]);
// //    }
// }


///////////////////////////////////////////////////////////////
/* ADDED BY MJC: Create a solver for a stenosis that uses nearly the same
   details as bound bif, but only solves for a connection between a parent
   and one daughter. For a stenosis, we employ this and add a loss term at
   the inlet of the daughter branch
   11/8/2019
 */
void bound_sten(double theta, double gamma,Tube *Arteries[], int parent, int D1, int WM)
{
    // MJC: Try to define the daughters here
    Tube* PV = Arteries[ parent]; // Parent vessel
    Tube* B1 = Arteries[ D1];      // Branch 1
    int N = PV->N;
    int j = 1;
    int ok = false;
    const int ntrial = 100;
    
    double k1[2];
    double k2[2];
    double k3[2];
    double k4[2];
        
        // Define the stenosis parameters based on derivations
    double Kv,Kt,Ku,Ls,As,Ap;//,sten_factor;
        // For now, just set values based on literature
        // MJC: for flow through a web-like stenoses, set Kt=0
        Kt = 0.0;//1.52; // Dissapation due to turbulent forces
        Ku = 1.2; // Dissapation due to inertial forces
        Ls = (PV->L)/4.0;
    //    fprintf(stdout,"ERROR NOT HERE\n");
//        sten_factor = 0.10; // Describes what percent stenosis
        Ap = M_PI*sq(PV->Aold[N]);//M_PI*sq(PV->rtop);
        As = Ap*(1.0-PV->sten_factor);
        Kv = 16.0*(Ls/As)*sq(Ap/As);
        
        
    //    fprintf(stdout,"Stenosis Parameters: Kt:%lf Ku:%lf  Kv:%lf \n  Ls:%lf  Sten:%lf  Ap:%lf As:%lf\n",Kt,Ku,Kv,Ls,sten_factor,Ap,As);

        double dp1,dp2,dp3;
    if (PV->sten_factor==1.0) {
        dp1=0.0; dp2=0.0; dp3=0.0;
    }
    else
    {
        dp1 = mu*Kv/(2.0*M_PI);
        dp2 = rho*Kt/2.0;
        dp3 = rho*Ku*Ls;
    }
    double delta_PA,delta_PB,term1A,term2A,term3A,term1B,term2B,term3B;
    
    
    k1[0] = PV->Qold[N]  + theta*(PV->R2h[N-1]) + gamma*(PV->S2h[N-1]);
    k1[1] = B1->Qold[0] - theta*(B1->R2h[0]) + gamma*(B1->S2h[0]);
    
    k2[0] = PV->Aold[N] + theta*(PV->R1h[N-1]);
    k2[1] = B1->Aold[0] - theta*(B1->R1h[0]);
    
    k3[0]   = PV->Qh[N-1]/2.0;
    k3[1]   = B1->Qh[0]/2.0;
    
    k4[0]   = PV->Ah[N-1]/2.0;
    k4[1]   = B1->Ah[0]/2.0;
    
    double xb[12];
    
    // The approximative initial guesses are applied.
    xb[ 0] =  PV->Qh[N-1];                      //Initial guess for Q1_xb n+1
    xb[ 1] = (PV->Qold[N-1] + PV->Qold[N])/2.0;       //Initial guess for Q1_xb^n+0.5
    xb[ 2] =  PV->Qold[N];                      //Initial guess for Q1_xb+0.5 n+0.5
    xb[ 3] =  B1->Qh[0];                    //Initial guess for Q2_xb n+1
    xb[ 4] = (B1->Qold[0] + B1->Qold[1])/2.0; //Initial guess for Q2_xb n+0.5
    xb[ 5] =  B1->Qold[0];                  //Initial guess for Q2_xb+0.5 n+0.5
    xb[ 6] =  PV->Ah[N-1];                      //Initial guess for A1_xb n+1
    xb[ 7] = (PV->Aold[N-1] + PV->Aold[N])/2.0;       //Initial guess for A1_xb^n+0.5
    xb[ 8] =  PV->Aold[N];                      //Initial guess for A1_xb+0.5 n+0.5
    xb[ 9] =  B1->Ah[0];                    //Initial guess for A2_xb n+1
    xb[10] = (B1->Aold[0] + B1->Aold[1])/2.0; //Initial guess for A2_xb n+0.5
    xb[11] =  B1->Aold[0];                  //Initial guess for A2_xb+0.5 n+0.5

    double k7nh  = 0.0;//LD->K_loss/2.0; //32*mu/(2*LD->rtop*rho*q);
    double k7n   = 0.0;//LD->K_loss/2.0; //32*mu/(2*RD->rtop*rho*q);

    
    // The residuals (fvec), and the Jacobian is determined, and if possible
    // the system of equations is solved.
    while (j <= ntrial && ok==false) // Find the zero
    {
        double fvec[12];
        
        // Characteristic Q residual at n+1
        fvec[0]  = k1[0]  - xb[0] -
        theta*(sq(xb[2])/xb[8] + PV->Bh(N,xb[8],WM)) +
        gamma*(F(xb[2],xb[8])+PV->dBdx1h(N,xb[8],WM));
        
        fvec[1]  = k1[1]  - xb[3] +
        theta*(sq(xb[5])/xb[11] + B1->Bh(-1,xb[11],WM)) +
        gamma*(F(xb[5],xb[11])  + B1->dBdx1h(-1,xb[11],WM));
        
        // Characteristic A residual at n+1
        fvec[2]  = - theta*xb[2] - xb[6]  + k2[0];
        fvec[3]  =   theta*xb[5] - xb[9]  + k2[1];
        
        // Flow residuals at n+1/2 (ghost points)
        fvec[4]  = - xb[ 1] + xb[ 2]/2.0  + k3[0];
        fvec[5]  = - xb[ 4] + xb[ 5]/2.0  + k3[1];
        
        // Area residuals at n+1/2 (ghost points)
        fvec[6]  = - xb[ 7] + xb[ 8]/2.0  + k4[0];
        fvec[7]  = - xb[10] + xb[11]/2.0  + k4[1];
        
        // Flow conservation residuals (n+1/2 and n+1)
        fvec[8]  = - xb[ 1] + xb[ 4];
        fvec[9]  = - xb[ 0] + xb[ 3];
        
        // Use these terms if you want to have Benoulli loss: u^2=(q/A)^2
        double sq211 = sq(xb[1]/xb[7]);
        double sq110 = sq(xb[0]/xb[6]);
        
        // Define the pressure loss terms. Order is viscous, turbulent, and inertial
        term1A = dp1 * xb[1]/sqrt(cu(xb[7]));
        term2A = (dp2 / sq(xb[7]))*sq(xb[7]/As - 1.0) * ab(xb[1])*xb[1]; // Using absolute value
        term3A = dp3 * (xb[1]-PV->Qold[N])/(xb[7]*2.0*gamma);
        delta_PA =  term1A + term2A + term3A;
        
    
        term1B = dp1 * xb[0]/sqrt(cu(xb[6]));
        term2B = (dp2 / sq(xb[6]))*sq(xb[6]/As - 1.0) * ab(xb[0])*xb[0];// Using absolute value
        term3B = dp3 *(xb[0]-PV->Qold[N])/(xb[6]* 2.0*gamma);
        delta_PB =  term1B + term2B + term3B;
        

            fvec[10] = - PV->P(N,xb[7],WM) + B1->P(0,xb[10],WM) + ab(k7nh)*sq211 + delta_PA;
            fvec[11] = - PV->P(N,xb[6],WM) + B1->P(0,xb[9],WM)  + ab(k7n)*sq110  + delta_PB;

//        fprintf(stdout,"FVEC IN ARTERIES: %lf %lf\n",fvec[10],fvec[11]);
//        fprintf(stdout,"AREA TERMS:%lf %lf\n",(xb[7]/As - 1.0),(xb[6]/As - 1.0));
//        fprintf(stdout,"dp:  %lf %lf %lf\n",dp1,dp2,dp3);
//        fprintf(stdout,"TERMS:  %lf %lf %lf\n",term1A,term2A,term3A);
//        fprintf(stdout,"TERMS:  %lf %lf %lf\n",term1B,term2B,term3B);

        
        
        double chi[10];
        
        // Here are the residuals for the characteristic matching for flow
        chi[0] = -2.0*theta*xb[ 2]/xb[8] + gamma*dFdQ(xb[8]);
        chi[2] =  2.0*theta*xb[ 5]/xb[11] + gamma*dFdQ(xb[11]);
        
        // Here are the residuals for the area characteristic matching
        chi[1] = theta*(sq(xb[2]/xb[8]) - PV->dBdAh(N,xb[8],WM)) +
                   gamma*(dFdA(xb[2],xb[8]) + PV->d2BdAdxh(N,xb[8],WM));
        
        chi[3] = theta*( -sq(xb[5]/xb[11]) + B1->dBdAh(-1,xb[11],WM)) +
                  gamma*(dFdA(xb[5],xb[11]) + B1->d2BdAdxh(-1,xb[11],WM));
        
        // Here is pressure conservation (n+1/2)
        chi[4]  = -PV->dPdA(N,xb[7],WM) - dp1*(-3.0*xb[1])/(2.0*sqrt(cu(xb[7]))) - dp3*(xb[1]-PV->Qold[N])/(sq(xb[7])*2.0*gamma) + dp2*ab(xb[1])*xb[1]*(-2.0/cu(xb[7]) * sq(xb[7]/As - 1.0) + (2.0/As/sq(xb[7]))*(xb[7]/As - 1.0)); // Using abs
        
        chi[5]  =  B1->dPdA(0,xb[10],WM);
        
        // Here is pressure conservation (n+1)
        chi[6] = -PV->dPdA(N,xb[6],WM) - dp1*(-3.0*xb[0])/(2.0*sqrt(cu(xb[6]))) - dp3*(xb[0]-PV->Qold[N])/(sq(xb[6])*2.0*gamma) + dp2*fabs(xb[0])*xb[0]*(-2.0/cu(xb[6])) * sq(xb[6]/As - 1.0) + (2/As/(sq(xb[6]))*(xb[6]/As - 1.0)); //Using abs
        chi[7] =  B1->dPdA(0,xb[9],WM);
        
        // Additional terms that arrive because of the pressure loss (depends on parent flow)
        chi[8] = dp1/sqrt(cu(xb[7]))+(dp2/sq(xb[7]))*sq(xb[7]/As - 1.0) * 2.0*sq(xb[1])/ab(xb[1]) + dp3/(xb[7]*2.0*gamma);
        
        chi[9] = dp1/sqrt(cu(xb[6]))+(dp2/sq(xb[6]))*sq(xb[6]/As - 1.0) * 2.0*sq(xb[0])/ab(xb[0]) + dp3/(xb[6]*2.0*gamma); // Using abs
        
                                                   
        
        for (int row = 0; row < 12; row++)
            for (int col = 0; col < 12; col++)
                fjac[row][col] = 0.0;
        
        // The Jacobian.
        // MJC: Ordered by column
        // Column 0:
        fjac[ 0][ 0] = -1.0;
        fjac[ 0][ 2] = chi[0];
        fjac[ 0][ 8] = chi[1];
        
        
        fjac[ 1][ 3] = -1.0;
        fjac[ 1][ 5] = chi[2];
        fjac[ 1][11] = chi[3];
        
        
        fjac[ 2][2] = -theta;
        fjac[ 2][6] = -1.0;
        
        fjac[ 3][5] =  theta;
        fjac[ 3][9] = -1.0;
        
        
        fjac[ 4][1] = -1.0;
        fjac[ 4][2] =  0.5;
        
        
        fjac[ 5][4] = -1.0;
        fjac[ 5][5] =  0.5;
        
        
        fjac[ 6][ 7] = -1.0;
        fjac[ 6][ 8] = 0.5;
        
        
        fjac[7][10] = -1.0;
        fjac[7][11] = 0.5;
        
        
        fjac[8][1] = -1.0;
        fjac[8][4] =  1.0;
        
        
        fjac[9][0]  = -1.0;
        fjac[ 9][3] =  1.0;
        
        fjac[10][ 1] = chi[8]; //Loss
        fjac[10][ 7] = chi[4];
        fjac[10][10] = chi[5];
        
        fjac[11][ 0] = chi[9]; //Loss
        fjac[11][ 6] = chi[6]; // Might need an if statement for negative flow here
        fjac[11][ 9] = chi[7]; // Might need an if statement for negative flow here
        
        
//        if(j==1){
//            printf("X AT t=%lf\n\n",t);
//            for (int i=0; i<12; i++)
//            {
//             printf("%lf ",xb[i]);
//            }
//            printf("\n\n");
//            printf("JACOBIAN\n");
//            for (int i=0; i<12; i++)
//            {
//                for (int j=0; j<12; j++) {
//                    printf("%lf ",fjac[i][j]);
//                }
//                printf("\n");
//            }
//            printf("\n\n");
//        }

        
//        fprintf(stdout,"AREA DIFFERENCE 1: %lf\n",xb[6]-xb[9]);
        // Check whether solution is close enough. If not run the loop again.
        // int ch = zero (xb, 18, 1.0e-4, 1.0e-4, fvec, fjac);
        int ch = zero (xb, 12, 1.0e-8, 1.0e-8, fvec, fjac);
        if (ch == 1) ok = true;
//        fprintf(stdout,"AREA DIFFERENCE 2: %lf AT t=%lf  \n",xb[6]-xb[9],t);

        j = j+1;
    }
    
    // Solutions is applied, and right boundary is updated.
    PV->Anew[N] = xb[ 6];
    PV->Qnew[N] = xb[ 0];
    B1->Anew[0] = xb[ 9];
    B1->Qnew[0] = xb[ 3];

      if (j >=ntrial) error ("arteries.C","Root not found in the bifurcation");
}



void bound_bif (double theta, double gamma, Tube *Arteries[], int parent, int D1, int D2, int WM)
{
    // MJC: Try to define the daughters here
    Tube* PV  = Arteries[ parent];
    Tube* B1  = Arteries[ D1];
    Tube* B2  = Arteries[ D2];
    
    
  int N = PV->N;
  double PN;
  int j = 1;
  int ok = false;
  const int ntrial = 500;
    

    double k1[3];
    double k2[3];
    double k3[3];
    double k4[3];
    
  k1[0] = PV->Qold[N]  + theta*(PV->R2h[N-1]) + gamma*(PV->S2h[N-1]);
  k1[1] = B1->Qold[0] - theta*(B1->R2h[0]) + gamma*(B1->S2h[0]);
  k1[2]= B2->Qold[0] - theta*(B2->R2h[0]) + gamma*(B2->S2h[0]);
    

  k2[0] = PV->Aold[N] + theta*(PV->R1h[N-1]);
  k2[1] = B1->Aold[0] - theta*(B1->R1h[0]);
  k2[2] = B2->Aold[0] - theta*(B2->R1h[0]);

  k3[0] = PV->Qh[N-1]/2.0;
  k3[1] = B1->Qh[0]/2.0;
  k3[2] = B2->Qh[0]/2.0;

  k4[0]  = PV->Ah[N-1]/2.0;
  k4[1]  = B1->Ah[0]/2.0;
  k4[2]  = B2->Ah[0]/2.0;

  double xb[6*3];

  // The approximative initial guesses are applied.
  xb[ 0] =  PV->Qh[N-1];                        //Initial guess for Q1_xb n+1
  xb[ 1] = (PV->Qold[N-1] + PV->Qold[N])/2.0;       //Initial guess for Q1_xb^n+0.5
  xb[ 2] =  PV->Qold[N];                        //Initial guess for Q1_xb+0.5 n+0.5
  xb[ 3] =  B1->Qh[0];                      //Initial guess for Q2_xb n+1
  xb[ 4] = (B1->Qold[0] + B1->Qold[1])/2.0; //Initial guess for Q2_xb n+0.5
  xb[ 5] =  B1->Qold[0];                    //Initial guess for Q2_xb+0.5 n+0.5
  xb[ 6] =  B2->Qh[0];                      //Initial guess for Q3_xb n+1
  xb[ 7] = (B2->Qold[0] + B2->Qold[1])/2.0; //Initial guess for Q3_xb n+0.5
  xb[ 8] =  B2->Qold[0];                    //Initial guess for Q3_xb+0.5 n+0.5
  xb[ 9] =  PV->Ah[N-1];                        //Initial guess for A1_xb n+1
  xb[10] = (PV->Aold[N-1] + PV->Aold[N])/2.0;       //Initial guess for A1_xb^n+0.5
  xb[11] =  PV->Aold[N];                        //Initial guess for A1_xb+0.5 n+0.5
  xb[12] =  B1->Ah[0];                      //Initial guess for A2_xb n+1
  xb[13] = (B1->Aold[0] + B1->Aold[1])/2.0; //Initial guess for A2_xb n+0.5
  xb[14] =  B1->Aold[0];                    //Initial guess for A2_xb+0.5 n+0.5
  xb[15] =  B2->Ah[0];                      //Initial guess for A3_xb n+1
  xb[16] = (B2->Aold[0] + B2->Aold[1])/2.0; //Initial guess for A3_xb n+0.5
  xb[17] =  B2->Aold[0];                    //Initial guess for A3_xb+0.5 n+0.5

  double k7nh  = 0.0;
  double k7n   = 0.0;
  double k7anh = 0.0;
  double k7an  = 0.0;

  // The residuals (fvec), and the Jacobian is determined, and if possible
  // the system of equations is solved.
  while (j <= ntrial && ok==false) // Find the zero
  {
    double fvec[18];
    // The residuals.
      
      // Characteristic Q residual at n+1
    fvec[0]  = k1[0]  - xb[0] -
    theta*(sq(xb[2])/xb[11] + PV->Bh(N,xb[11],WM)) +
    gamma*(F(xb[2],xb[11])+PV->dBdx1h(N,xb[11],WM));

    fvec[1]  = k1[1]  - xb[3] +
                 theta*(sq(xb[5])/xb[14] + B1->Bh(-1,xb[14],WM)) +
                 gamma*(F(xb[5],xb[14])  + B1->dBdx1h(-1,xb[14],WM));

    fvec[2]  = k1[2] - xb[6] +
    theta*(sq(xb[8])/xb[17] + B2->Bh(-1,xb[17],WM)) +
    gamma*(F(xb[8],xb[17])  + B2->dBdx1h(-1,xb[17],WM));

    // Characteristic A residual at n+1
    fvec[3]  = - theta*xb[2] - xb[ 9] + k2[0];
    fvec[4]  =   theta*xb[5] - xb[12] + k2[1];
    fvec[5]  =   theta*xb[8] - xb[15] + k2[2];
      
    // Flow residuals at n+1/2 (ghost points)
    fvec[6]  = - xb[ 1] + xb[ 2]/2.0 + k3[0];
    fvec[7]  = - xb[ 4] + xb[ 5]/2.0 + k3[1];
    fvec[8]  = - xb[ 7] + xb[ 8]/2.0 + k3[2];
      
    // Area residuals at n+1/2 (ghost points)
    fvec[9]  = - xb[10] + xb[11]/2.0 + k4[0];
    fvec[10] = - xb[13] + xb[14]/2.0 + k4[1];
    fvec[11] = - xb[16] + xb[17]/2.0 + k4[2];
      
    // Flow conservation residuals (n+1/2 and n+1)
    fvec[12] = - xb[ 1] + xb[ 4] + xb[7];
    fvec[13] = - xb[ 0] + xb[ 3] + xb[6];

    // Pressure continuity at the n+1/2 time step
    PN    = PV->P(N,xb[10],WM);
    double u_n_half = sq(xb[1]/xb[10]);

    // The if statements here only matter if a minor loss
    // is included
      fvec[14] =  - PN + B1->P(0,xb[13],WM) + ab(k7nh)*u_n_half;
      fvec[15] =  - PN + B2->P(0,xb[16],WM) + ab(k7anh)*u_n_half;

    // Pressure continuity at the n+1 time step
    PN    = PV->P(N,xb[9],WM);
    double u_n_1 = sq(xb[0]/xb[9]);
      fvec[16] = - PN + B1->P(0,xb[12],WM) + ab(k7n)*u_n_1;
      fvec[17] = - PN + B2->P(0,xb[15],WM) + ab(k7an)*u_n_1;

    for (int row = 0; row < 18; row++)
      for (int col = 0; col < 18; col++)
          fjac[row][col] = 0.0;
//        fjac[row][col] = 1e-12;

      
      
    // The Jacobian.
      
      double chi[12];
      
       // Here are the residuals for the characteristic matching for flow
      chi[0] = -2.0*theta*xb[ 2]/xb[11] + gamma*dFdQ(xb[11]);
      chi[2] =  2.0*theta*xb[ 5]/xb[14] + gamma*dFdQ(xb[14]);
      chi[4] =  2.0*theta*xb[ 8]/xb[17] + gamma*dFdQ(xb[17]);
      
      // Here are the residuals for the area characteristic matching
      chi[1] = theta*(sq(xb[2]/xb[11]) - PV->dBdAh(N,xb[11],WM)) +
                 gamma*(dFdA(xb[2],xb[11]) + PV->d2BdAdxh(N,xb[11],WM));
      
      chi[3] = theta*( -sq(xb[5]/xb[14]) + B1->dBdAh(-1,xb[14],WM)) +
                gamma*(dFdA(xb[5],xb[14]) + B1->d2BdAdxh(-1,xb[14],WM));
      
      chi[5] = theta*( -sq(xb[8]/xb[17]) + B2->dBdAh(-1,xb[17],WM)) +
                gamma*(dFdA(xb[8],xb[17]) + B2->d2BdAdxh(-1,xb[17],WM));
      
      // Here is pressure conservation
      chi[6]  = -PV->dPdA(N,xb[10],WM) + sq(xb[1])/cu(xb[10])*(-2.0*ab(k7nh)); //Loss term
      chi[7]  = B1->dPdA(0,xb[13],WM);
      chi[8] = B2->dPdA(0,xb[16],WM);
      
      chi[9] = -PV->dPdA(N,xb[9],WM) + sq(xb[0])/cu(xb[9])*(-2.0*ab(k7nh)); //Loss term
      chi[10] = B1->dPdA(0,xb[12],WM);
      chi[11] = B2->dPdA(0,xb[15],WM);

      
      //NEW JACOBIAN
            // Order is [row][column]
            fjac[ 0][ 0]  = -1.0;
            fjac[ 0][ 2] = chi[0];
            fjac[ 0][11] = chi[1];
            
            fjac[ 1][ 3] = -1.0;
            fjac[ 1][ 5] = chi[2];
            fjac[ 1][14] = chi[3];

            fjac[ 2][ 6] = -1.0;
            fjac[ 2][ 8] = chi[4];
            fjac[ 2][17] = chi[5];

            fjac[ 3][ 2] = -theta;
            fjac[ 3][ 9] = -1.0;
            
            fjac[ 4][ 5] = theta;
            fjac[ 4][12] = -1.0;
            
            fjac[ 5][ 8] = theta;
            fjac[ 5][15] = -1.0;
            
            fjac[ 6][ 1] = -1.0;
            fjac[ 6][ 2] = 0.5;
            
            fjac[ 7][ 4] = -1.0;
            fjac[ 7][ 5] = 0.5;
            
            fjac[ 8][ 7] = -1.0;
            fjac[ 8][ 8] = 0.5;
            
            fjac[ 9][ 10] = -1.0;
            fjac[ 9][ 11] = 0.5;
            
            fjac[10][13] = -1.0;
            fjac[10][14] = 0.5;
            
            fjac[11][16] = -1.0;
            fjac[11][17] = 0.5;
            
            fjac[12][ 1] = -1.0;
            fjac[12][ 4] = 1.0;
            fjac[12][ 7] = 1.0;
            
            fjac[13][ 0] = -1.0;
            fjac[13][ 3] = 1.0;
            fjac[13][ 6] = 1.0;
            
      //      fjac[14][ 1] = J_lossP1;
            fjac[14][10] = chi[6];
            fjac[14][13] = chi[7];
                                    
                                    
      //      fjac[15][ 1] = J_lossP2;
            fjac[15][10] = chi[6];
            fjac[15][16] = chi[8];
            
      //      fjac[16][ 0] = J_lossP3;
            fjac[16][ 9] = chi[9];
            fjac[16][12] = chi[10];
            
            
      //      fjac[17][ 0] = J_lossP4;
            fjac[17][ 9] = chi[9];
            fjac[17][15] = chi[11];
    // Check whether solution is close enough. If not run the loop again.
    int ch = zero (xb, 18, 1.0e-8, 1.0e-8, fvec, fjac);
    if (ch == 1) ok = true;
//      fprintf(stdout,"parent %d new flow:\n",parent);
//      fprintf(stdout,"%lf %lf %lf\n",xb[0],xb[3],xb[6]);
//      fprintf(stdout,"new area:\n");
//      fprintf(stdout,"%lf %lf %lf\n",xb[9],xb[12],xb[15]);
    j = j+1;
  }

  // Solutions is applied, and right boundary is updated.
  PV->Anew[N] = xb[ 9];
  PV->Qnew[N] = xb[ 0];
  B1->Anew[0] = xb[12];
  B1->Qnew[0] = xb[ 3];
  B2->Anew[0] = xb[15];
  B2->Qnew[0] = xb[ 6];
    
//    fprintf(stdout,"IN BIF: Parent Qold: %lf Daughter 1 Qold: %lf Daughter 2 Qold: %lf\n",PV->Qold[N],B1->Qold[0],B2->Qold[0]);
//    fprintf(stdout,"IN BIF: Parent Qold: %lf Daughter 1 Qold: %lf Daughter 2 Qold: %lf\n",PV->Qold[N],B1->Qold[0],B2->Qold[0]);
    
    if (j >=ntrial) {error ("arteries.C","Root not found in the bifurcation");
        exit(1);}
}


// void bound_sten_bif (double theta, double gamma, Tube *Arteries[], int parent, int D1, int D2, int WM)
// {
//     // MJC: Try to define the daughters here
//     Tube* PV   = Arteries[ parent];
//     Tube* B1  = Arteries[ D1];
//     Tube* B2  = Arteries[ D2];
// 
//     
//   int N = PV->N;
//   int j = 1;
//   int ok = false;
//   const int ntrial = 20000;
// 
//     // Define the stenosis parameters based on derivations
//     // last edited by MJC 5/27/2020
//         double C,Kt,Ku,Ls;//,Ap,Rs,La,As,Rp,sten_factor;
//             // For now, just set values based on literature
//             Kt = 1.52; // Dissapation due to turbulent forces
//             Ku = 1.2; // Dissapation due to inertial forces
//             Ls = PV->sten_length; // Define this outside in future versions
//             C = (1.0-PV->sten_factor); // Define degree of stenosis as As/Ap, or 1-sten_factor
// //            Ap = PV->Aold[N];//M_PI*sq(PV->rtop); // USING DYNAMIC AREA
// //            As = Ap*C;
// //            Rp = sqrt(Ap/M_PI);
// //            Rs = sqrt(C*Ap/M_PI);
// //            La = 0.83*Ls + 3.28*Rs;
// //            fprintf(stdout,"Stenosis Parameters: Kt:%lf Ku:%lf Ls:%lf  Sten:%lf\n",Kt,Ku,Ls,C);
// 
//         double dp1,dp2,dp3;
//         if (C==1.0) {
//             dp1=0.0; dp2=0.0; dp3=0.0;
//         }
//         else
//         {
//             dp1 = (8.0*Ls*mu*M_PI)/sq(C); // From Karniadakis group 2019
//             dp2 = rho*Kt*sq(1.0/C - 1.0)/2.0;
//             dp3 = rho*Ku*Ls;
//         }
// //        fprintf(stdout,"dp terms: dp1:%lf dp2:%lf dp3:%lf\n",dp1,dp2,dp3);
// 
// 
//     double k1[3];
//     double k2[3];
//     double k3[3];
//     double k4[3];
//     
//   k1[0] = PV->Qold[N]  + theta*(PV->R2h[N-1]) + gamma*(PV->S2h[N-1]);
//   k1[1] = B1->Qold[0] - theta*(B1->R2h[0]) + gamma*(B1->S2h[0]);
//   k1[2] = B2->Qold[0] - theta*(B2->R2h[0]) + gamma*(B2->S2h[0]);
//     
// 
//   k2[0] = PV->Aold[N] + theta*(PV->R1h[N-1]);
//   k2[1] = B1->Aold[0] - theta*(B1->R1h[0]);
//   k2[2] = B2->Aold[0] - theta*(B2->R1h[0]);
// 
//   k3[0] = PV->Qh[N-1]/2.0;
//   k3[1] = B1->Qh[0]/2.0;
//   k3[2] = B2->Qh[0]/2.0;
// 
//   k4[0]  = PV->Ah[N-1]/2.0;
//   k4[1]  = B1->Ah[0]/2.0;
//   k4[2]  = B2->Ah[0]/2.0;
// 
//   double xb[6*3];
// 
//   // The approximative initial guesses are applied.
//   xb[ 0] =  PV->Qh[N-1];                        //Initial guess for Q1_xb n+1
//   xb[ 1] = (PV->Qold[N-1] + PV->Qold[N])/2.0;       //Initial guess for Q1_xb^n+0.5
//   xb[ 2] =  PV->Qold[N];                        //Initial guess for Q1_xb+0.5 n+0.5
//   xb[ 3] =  B1->Qh[0];                      //Initial guess for Q2_xb n+1
//   xb[ 4] = (B1->Qold[0] + B1->Qold[1])/2.0; //Initial guess for Q2_xb n+0.5
//   xb[ 5] =  B1->Qold[0];                    //Initial guess for Q2_xb+0.5 n+0.5
//   xb[ 6] =  B2->Qh[0];                      //Initial guess for Q3_xb n+1
//   xb[ 7] = (B2->Qold[0] + B2->Qold[1])/2.0; //Initial guess for Q3_xb n+0.5
//   xb[ 8] =  B2->Qold[0];                    //Initial guess for Q3_xb+0.5 n+0.5
//   xb[ 9] =  PV->Ah[N-1];                        //Initial guess for A1_xb n+1
//   xb[10] = (PV->Aold[N-1] + PV->Aold[N])/2.0;       //Initial guess for A1_xb^n+0.5
//   xb[11] =  PV->Aold[N];                        //Initial guess for A1_xb+0.5 n+0.5
//   xb[12] =  B1->Ah[0];                      //Initial guess for A2_xb n+1
//   xb[13] = (B1->Aold[0] + B1->Aold[1])/2.0; //Initial guess for A2_xb n+0.5
//   xb[14] =  B1->Aold[0];                    //Initial guess for A2_xb+0.5 n+0.5
//   xb[15] =  B2->Ah[0];                      //Initial guess for A3_xb n+1
//   xb[16] = (B2->Aold[0] + B2->Aold[1])/2.0; //Initial guess for A3_xb n+0.5
//   xb[17] =  B2->Aold[0];                    //Initial guess for A3_xb+0.5 n+0.5
// 
//    
// 
// 
// //        sum = 0.0;
// //        stpmax = 0.0;
// //        // Calulcate stpmax for linesearch
// //        for (int i = 0; i<18; i++)
// //            sum+=xb[i]*xb[i];
// //        stpmax = stpmaxMAX*fmax(sqrt(sum),double(18));
// 
//     // MJC: test initial point
//     // The residuals.
// //         double fvec[18];
// 
// 
// //       // Characteristic Q residual at n+1
// //     fvec[0]  = k1[0]  - xb[0] -
// //     theta*(sq(xb[2])/xb[11] + PV->Bh(N,xb[11],WM)) +
// //     gamma*(F(xb[2],xb[11])+PV->dBdx1h(N,xb[11],WM));
// // 
// //     fvec[1]  = k1[1]  - xb[3] +
// //                  theta*(sq(xb[5])/xb[14] + B1->Bh(-1,xb[14],WM)) +
// //                  gamma*(F(xb[5],xb[14])  + B1->dBdx1h(-1,xb[14],WM));
// // 
// //     fvec[2]  = k1[2] - xb[6] +
// //     theta*(sq(xb[8])/xb[17] + B2->Bh(-1,xb[17],WM)) +
// //     gamma*(F(xb[8],xb[17])  + B2->dBdx1h(-1,xb[17],WM));
// // 
// //     // Characteristic A residual at n+1
// //     fvec[3]  = - theta*xb[2] - xb[ 9] + k2[0];
// //     fvec[4]  =   theta*xb[5] - xb[12] + k2[1];
// //     fvec[5]  =   theta*xb[8] - xb[15] + k2[2];
// // 
// //     // Flow residuals at n+1/2 (ghost points)
// //     fvec[6]  = - xb[ 1] + xb[ 2]/2.0 + k3[0];
// //     fvec[7]  = - xb[ 4] + xb[ 5]/2.0 + k3[1];
// //     fvec[8]  = - xb[ 7] + xb[ 8]/2.0 + k3[2];
// // 
// //     // Area residuals at n+1/2 (ghost points)
// //     fvec[9]  = - xb[10] + xb[11]/2.0 + k4[0];
// //     fvec[10] = - xb[13] + xb[14]/2.0 + k4[1];
// //     fvec[11] = - xb[16] + xb[17]/2.0 + k4[2];
// // 
// //     // Flow conservation residuals (n+1/2 and n+1)
// //     fvec[12] = - xb[ 1] + xb[ 4] + xb[7];
// //     fvec[13] = - xb[ 0] + xb[ 3] + xb[6];
// // 
// //       double delta_PA,delta_PB,term1A,term2A,term3A,term1B,term2B,term3B;
// //     // Define the pressure loss terms. Order is viscous, turbulent, and inertial
// //         term1A = dp1 * xb[1]/sq(xb[10]);
// // 
// //         term2A = (dp2 / sq(xb[10])) * ab(xb[1])*xb[1]; // Using absolute value
// // 
// //         term3A = (dp3/xb[10]) * (xb[1]-(PV->Qold[N]))/(2.0*gamma);
// // 
// //         delta_PA =  term1A + term2A + term3A;
// // 
// // 
// //         term1B = dp1 * xb[0]/sq(xb[9]);
// // 
// //         term2B = (dp2 / sq(xb[9])) * ab(xb[0])*xb[0]; // Using absolute value
// // 
// //         term3B = (dp3/xb[9]) * (xb[0]-(PV->Qold[N]))/(2.0*gamma);
// // 
// //         delta_PB =  term1B + term2B + term3B;
// // 
// //     // Pressure loss at the n+1/2 time step
// //       fvec[14] = - PV->P(N,xb[ 10],WM) + B1->P(0,xb[13],WM) + delta_PA;
// //       fvec[15] = - PV->P(N,xb[ 10],WM) + B2->P(0,xb[16],WM) + delta_PA;
// // 
// //     // Pressure continuity at the n+1 time step
// //     fvec[16] = - PV->P(N,xb[ 9],WM) + B1->P(0,xb[12],WM) + delta_PB;
// //     fvec[17] = - PV->P(N,xb[ 9],WM) + B2->P(0,xb[15],WM) + delta_PB;
// // 
// // 
// //       // MJC: test initial point and exit if close enough
// //     // Only use for linesearch
// // //     test=0.0;
// // //     for (int i=0; i<18; i++)
// // //         if (fabs(fvec[i])>test)
// // //             test = fvec[i];
// // //     if (test<tol_res*0.1)
// // //         ok = true;
// 
//    // The residuals (fvec), and the Jacobian is determined, and if possible
//   // the system of equations is solved.
// 
// 
//    
//   
//   while (j <= ntrial && ok==false) // Find the zero
//   {
//    double fvec[18];
//     // The residuals.
// 
//       // Characteristic Q residual at n+1
//     fvec[0]  = k1[0]  - xb[0] -
//     theta*(sq(xb[2])/xb[11] + PV->Bh(N,xb[11],WM)) +
//     gamma*(F(xb[2],xb[11])+PV->dBdx1h(N,xb[11],WM));
// 
//     fvec[1]  = k1[1]  - xb[3] +
//                  theta*(sq(xb[5])/xb[14] + B1->Bh(-1,xb[14],WM)) +
//                  gamma*(F(xb[5],xb[14])  + B1->dBdx1h(-1,xb[14],WM));
// 
//     fvec[2]  = k1[2] - xb[6] +
//     theta*(sq(xb[8])/xb[17] + B2->Bh(-1,xb[17],WM)) +
//     gamma*(F(xb[8],xb[17])  + B2->dBdx1h(-1,xb[17],WM));
// 
//     // Characteristic A residual at n+1
//     fvec[3]  = - theta*xb[2] - xb[ 9] + k2[0];
//     fvec[4]  =   theta*xb[5] - xb[12] + k2[1];
//     fvec[5]  =   theta*xb[8] - xb[15] + k2[2];
// 
//     // Flow residuals at n+1/2 (ghost points)
//     fvec[6]  = - xb[ 1] + xb[ 2]/2.0 + k3[0];
//     fvec[7]  = - xb[ 4] + xb[ 5]/2.0 + k3[1];
//     fvec[8]  = - xb[ 7] + xb[ 8]/2.0 + k3[2];
// 
//     // Area residuals at n+1/2 (ghost points)
//     fvec[9]  = - xb[10] + xb[11]/2.0 + k4[0];
//     fvec[10] = - xb[13] + xb[14]/2.0 + k4[1];
//     fvec[11] = - xb[16] + xb[17]/2.0 + k4[2];
// 
//     // Flow conservation residuals (n+1/2 and n+1)
//     fvec[12] = - xb[ 1] + xb[ 4] + xb[7];
//     fvec[13] = - xb[ 0] + xb[ 3] + xb[6];
// 
//       double delta_PA,delta_PB,term1A,term2A,term3A,term1B,term2B,term3B;
//     // Define the pressure loss terms. Order is viscous, turbulent, and inertial
//         term1A = dp1 * xb[1]/sq(xb[10]);
// 
//         term2A = (dp2 / sq(xb[10])) * ab(xb[1])*xb[1]; // Using absolute value
// 
//         term3A = (dp3/xb[10]) * (xb[1]-(PV->Qold[N]))/(2.0*gamma);
// 
//         delta_PA =  term1A + term2A + term3A;
// 
// 
//         term1B = dp1 * xb[0]/sq(xb[9]);
// 
//         term2B = (dp2 / sq(xb[9])) * ab(xb[0])*xb[0]; // Using absolute value
// 
//         term3B = (dp3/xb[9]) * (xb[0]-(PV->Qold[N]))/(2.0*gamma);
// 
//         delta_PB =  term1B + term2B + term3B;
// 
//     // Pressure loss at the n+1/2 time step
//       fvec[14] = - PV->P(N,xb[ 10],WM) + B1->P(0,xb[13],WM) + delta_PA;
//       fvec[15] = - PV->P(N,xb[ 10],WM) + B2->P(0,xb[16],WM) + delta_PA;
// 
//     // Pressure continuity at the n+1 time step
//     fvec[16] = - PV->P(N,xb[ 9],WM) + B1->P(0,xb[12],WM) + delta_PB;
//     fvec[17] = - PV->P(N,xb[ 9],WM) + B2->P(0,xb[15],WM) + delta_PB;
// 
// 
// 
//     for (int row = 0; row < 18; row++)
//       for (int col = 0; col < 18; col++)
//         fjac[row][col] = 0.0;
// 
// 
// 
//     // The Jacobian.
// 
//       double chi[14];
// 
//        // Here are the residuals for the characteristic matching for flow
//       chi[0] = -2.0*theta*xb[ 2]/xb[11] + gamma*dFdQ(xb[11]);
//       chi[2] =  2.0*theta*xb[ 5]/xb[14] + gamma*dFdQ(xb[14]);
//       chi[4] =  2.0*theta*xb[ 8]/xb[17] + gamma*dFdQ(xb[17]);
// 
//       // Here are the residuals for the area characteristic matching
//       chi[1] = theta*(sq(xb[2]/xb[11]) - PV->dBdAh(N,xb[11],WM)) +
//                  gamma*(dFdA(xb[2],xb[11]) + PV->d2BdAdxh(N,xb[11],WM));
// 
//       chi[3] = theta*( -sq(xb[5]/xb[14]) + B1->dBdAh(-1,xb[14],WM)) +
//                 gamma*(dFdA(xb[5],xb[14]) + B1->d2BdAdxh(-1,xb[14],WM));
// 
//       chi[5] = theta*( -sq(xb[8]/xb[17]) + B2->dBdAh(-1,xb[17],WM)) +
//                 gamma*(dFdA(xb[8],xb[17]) + B2->d2BdAdxh(-1,xb[17],WM));
// 
//       // Here are the pressure loss terms at n+1/2
//       // Derivative with respect to flow
//       chi[6]  = dp1/sq(xb[10]) + (dp2/sq(xb[10]))*2.0*xb[1]/ab(xb[1])
//                            + dp3/(2.0*gamma*xb[10]);
// 
//       // Derivative with respect to area
//       chi[7]  = -PV->dPdA(N,xb[10],WM) - 2.0*dp1*xb[1]/cu(xb[10])
//       - 2.0*dp2*ab(xb[1])*xb[1]/cu(xb[10]) - (dp3/sq(xb[10])) * (xb[1]-(PV->Qold[N]))/(2.0*gamma);
// 
//       chi[8]  = B1->dPdA(0,xb[13],WM);
// 
//       chi[9]  = B2->dPdA(0,xb[16],WM);
// 
// 
//       // Here are the pressure loss terms at n+1
//       // Derivative with respect to flow
//       chi[10]  = dp1/sq(xb[9]) + (dp2/sq(xb[9]))*2.0*xb[0]/ab(xb[0])
//                            + dp3/(2*gamma*xb[9]);
// 
//       // Derivative with respect to area
//       chi[11]  = -PV->dPdA(N,xb[9],WM) - 2.0*dp1*xb[0]/cu(xb[9])
//       - 2.0*dp2*ab(xb[0])*xb[0]/cu(xb[9]) - (dp3/sq(xb[9])) * (xb[0]-(PV->Qold[N]))/(2.0*gamma);
// 
//       chi[12] = B1->dPdA(0,xb[12],WM);
// 
//       chi[13] = B2->dPdA(0,xb[15],WM);
// 
// 
//       //NEW JACOBIAN
//             // Order is [row][column]
//             fjac[ 0][ 0]  = -1.0;
//             fjac[ 0][ 2] = chi[0];
//             fjac[ 0][11] = chi[1];
// 
//             fjac[ 1][ 3] = -1.0;
//             fjac[ 1][ 5] = chi[2];
//             fjac[ 1][14] = chi[3];
// 
//             fjac[ 2][ 6] = -1.0;
//             fjac[ 2][ 8] = chi[4];
//             fjac[ 2][17] = chi[5];
// 
//             fjac[ 3][ 2] = -theta;
//             fjac[ 3][ 9] = -1.0;
// 
//             fjac[ 4][ 5] = theta;
//             fjac[ 4][12] = -1.0;
// 
//             fjac[ 5][ 8] = theta;
//             fjac[ 5][15] = -1.0;
// 
//             fjac[ 6][ 1] = -1.0;
//             fjac[ 6][ 2] = 0.5;
// 
//             fjac[ 7][ 4] = -1.0;
//             fjac[ 7][ 5] = 0.5;
// 
//             fjac[ 8][ 7] = -1.0;
//             fjac[ 8][ 8] = 0.5;
// 
//             fjac[ 9][ 10] = -1.0;
//             fjac[ 9][ 11] = 0.5;
// 
//             fjac[10][13] = -1.0;
//             fjac[10][14] = 0.5;
// 
//             fjac[11][16] = -1.0;
//             fjac[11][17] = 0.5;
// 
//             fjac[12][ 1] = -1.0;
//             fjac[12][ 4] = 1.0;
//             fjac[12][ 7] = 1.0;
// 
//             fjac[13][ 0] = -1.0;
//             fjac[13][ 3] = 1.0;
//             fjac[13][ 6] = 1.0;
// 
//             fjac[14][ 1] = chi[6];
//             fjac[14][10] = chi[7];
//             fjac[14][13] = chi[8];
// 
// 
//             fjac[15][ 1] = chi[6];
//             fjac[15][10] = chi[7];
//             fjac[15][16] = chi[9];
// 
//             fjac[16][ 0] = chi[10];
//             fjac[16][ 9] = chi[11];
//             fjac[16][12] = chi[12];
// 
// 
//             fjac[17][ 0] = chi[10];
//             fjac[17][ 9] = chi[11];
//             fjac[17][15] = chi[13];
// 
// ///////////////////////////////////////////////////////////////////////////////////////////////////
//         // THIS IS FOR NEWTON'S METHOD
// //////////////////////////////////////////////////////////////////////////////////////////////////
// // Check whether solution is close enough. If not run the loop again.
//     int ch = zero (xb, 18, 1.0e-6, 1.0e-6, fvec, fjac);
//     if (ch == 1) ok = true;
// 
// 
// 
// ///////////////////////////////////////////////////////////////////////////////////////////////////
//         // THIS IS FOR THE GRADIENT DESCENT WITH LINESEARCH  
// //////////////////////////////////////////////////////////////////////////////////////////////////
//       /*
// 
//          // ADD NEWTON PARAMETERS w/ LINE SEARCH HERE
//        double tolx = 1e-6;    // Criteria for x
//        double tol_res = 1e-4; // Criteria for residual
//        double tol_min = 1e-6; // Criteria for gradient
//        double stpmaxMAX = 10000.0;
//        double alam, alamin, a, b, disc, f2, rhs1, rhs2, stpmax, slope, sum, temp, test, tmplam;
//        double alam2 = 0.0;
//        double den;
//        double fvec[18], pln[18], grad[18], xold[18], fold[18];
//        int indx[18];
//        double d;
//        double ALF = 1.0e-4;
//        double SSE_old,SSE_new;
// 
//       // Line-search procedures
//           // First, construct a gradient and old x/f values
//           for (int i=0; i<18; i++)
//           {
//               sum = 0.0;
//               for (int k=0; k<18; k++)
//                   sum += fjac[k][i]*fvec[k];
//               grad[i] = sum;
//               pln[i] = -fvec[i]; // Assign Newton direction
//               xold[i] = xb[i];
//               fold[i] = fvec[i];
//   //            fprintf(stdout,"grad: %lf p: %lf xold: %lf fold: %lf\n",grad[i],pln[i], xold[i],fold[i]);
//           }
//               ludcmp (fjac, 18, indx, &d);
//               lubksb (fjac, 18, indx, pln);
//               sum = 0.0;
//               for (int i=0; i<18; i++)
//               {
//                   sum+=pln[i]*pln[i];
//               }
//               // Scale step if too large
//               if (sum>stpmax) {
//       //            fprintf(stdout,"stpmax too large: %lf\n",sum);
//                   for (int i=0; i<18; i++)
//                   {
//                       pln[i]*=stpmax/sum;
// //                      fprintf(stdout,"p %d: %5.5lf",i,pln[i]);
//                   }
// //                  fprintf(stdout,"\n");
//               }
//               slope = 0.0;
//               // Construct the slope for the gradient
//               for (int i=0; i<18 ;i++)
//               {
//                   slope += grad[i]*pln[i];
//       //            fprintf(stdout,"g: %lf p: %lf slope: %lf\n",grad[i],pln[i], slope);
//               }
//               if (slope >= 0.0)
//               {
//                   fprintf(stdout,"Roundoff problem in lnsrch.\n");
//                   exit(1);
//               }
//               test=0.0;
//               for (int i=0;i< 18; i++) {
//                   temp=abs(pln[i])/fmax(abs(xold[i]),1.0);
//                   if (temp > test) test=temp;
//               }
//                 alamin=tolx/test;
//               alam=1.0;
//               alam2 = 0.0;
//               f2 = 0.0;
// 
//               // Construct initial SSE
//               SSE_old = 0.0;
//               for (int i=0; i<18; i++) SSE_old+=fvec[i]*fvec[i];
//               SSE_old*=0.5;
//               // Now iterate until the Armijo condition is satisfied
//               bool armijo_flag = false;
//               while (!armijo_flag) {
//                   // Take an initial full gradient step
//                   for (int i=0; i<18; i++)
//                   {
//                       xb[i] = xold[i] + alam*pln[i];
//                   }
//                   // Now recompute the residual values
//                   // Characteristic Q residual at n+1
//                   fvec[0]  = k1[0]  - xb[0] -
//                   theta*(sq(xb[2])/xb[11] + PV->Bh(N,xb[11],WM)) +
//                   gamma*(F(xb[2],xb[11])+PV->dBdx1h(N,xb[11],WM));
// 
//                   fvec[1]  = k1[1]  - xb[3] +
//                                theta*(sq(xb[5])/xb[14] + B1->Bh(-1,xb[14],WM)) +
//                                gamma*(F(xb[5],xb[14])  + B1->dBdx1h(-1,xb[14],WM));
// 
//                   fvec[2]  = k1[2] - xb[6] +
//                   theta*(sq(xb[8])/xb[17] + B2->Bh(-1,xb[17],WM)) +
//                   gamma*(F(xb[8],xb[17])  + B2->dBdx1h(-1,xb[17],WM));
// 
//                   // Characteristic A residual at n+1
//                   fvec[3]  = - theta*xb[2] - xb[ 9] + k2[0];
//                   fvec[4]  =   theta*xb[5] - xb[12] + k2[1];
//                   fvec[5]  =   theta*xb[8] - xb[15] + k2[2];
// 
//                   // Flow residuals at n+1/2 (ghost points)
//                   fvec[6]  = - xb[ 1] + xb[ 2]/2.0 + k3[0];
//                   fvec[7]  = - xb[ 4] + xb[ 5]/2.0 + k3[1];
//                   fvec[8]  = - xb[ 7] + xb[ 8]/2.0 + k3[2];
// 
//                   // Area residuals at n+1/2 (ghost points)
//                   fvec[9]  = - xb[10] + xb[11]/2.0 + k4[0];
//                   fvec[10] = - xb[13] + xb[14]/2.0 + k4[1];
//                   fvec[11] = - xb[16] + xb[17]/2.0 + k4[2];
// 
//                   // Flow conservation residuals (n+1/2 and n+1)
//                   fvec[12] = - xb[ 1] + xb[ 4] + xb[7];
//                   fvec[13] = - xb[ 0] + xb[ 3] + xb[6];
// 
//                   // Define the pressure loss terms. Order is viscous, turbulent, and inertial
//                       term1A = dp1 * xb[1]/sq(xb[10]);
// 
//                       term2A = (dp2 / sq(xb[10])) * ab(xb[1])*xb[1]; // Using absolute value
// 
//                       term3A = (dp3/xb[10]) * (xb[1]-(PV->Qold[N]))/(2.0*gamma);
// 
//                       delta_PA =  term1A + term2A + term3A;
// 
// 
//                       term1B = dp1 * xb[0]/sq(xb[9]);
// 
//                       term2B = (dp2 / sq(xb[9])) * ab(xb[0])*xb[0]; // Using absolute value
// 
//                       term3B = (dp3/xb[9]) * (xb[0]-(PV->Qold[N]))/(2.0*gamma);
// 
//                       delta_PB =  term1B + term2B + term3B;
// 
//                   // Pressure loss at the n+1/2 time step
//                     fvec[14] = - PV->P(N,xb[ 10],WM) + B1->P(0,xb[13],WM) + delta_PA;
//                     fvec[15] = - PV->P(N,xb[ 10],WM) + B2->P(0,xb[16],WM) + delta_PA;
// 
//                   // Pressure continuity at the n+1 time step
//                   fvec[16] = - PV->P(N,xb[ 9],WM) + B1->P(0,xb[12],WM) + delta_PB;
//                   fvec[17] = - PV->P(N,xb[ 9],WM) + B2->P(0,xb[15],WM) + delta_PB;
// 
//                   // Calculate SSE
//                   SSE_new = 0.0;
//                   for (int i=0; i<18; i++) SSE_new+=fvec[i]*fvec[i];
//                   SSE_new*=0.5;
//                   if (alam<alamin) {
//                       for (int i=0; i<18; i++) xb[i] = xold[i];
//                       armijo_flag = true;
//                   }
//                   else if (SSE_new <= SSE_old +ALF*alam*slope)
//                   {
//                       armijo_flag=true;
// 
//                   }
//                   else
//                   {
//                       if (alam == 1.0)
//                       {
//                           tmplam = -slope/(2.0*(SSE_new-SSE_old-slope));
//                       }
//                       else
//                       {
//                           rhs1 = SSE_new-SSE_old-alam*slope;
//                           rhs2 = f2-SSE_old-alam2*slope;
//                           a    = (rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
//                           b    = (-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
//                           if (a == 0.0) tmplam = -slope/(2.0*b);
//                           else {
//                               disc=b*b-3.0*a*slope;
//                               if (disc < 0.0) tmplam=0.5*alam;
//                               else if (b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
//                               else tmplam=-slope/(b+sqrt(disc));
//                               }
//                               if (tmplam>0.5*alam)
//                               {
// //                                  fprintf(stdout,"tmplam: %12.5lf,  tmplam_new: %12.5lf\n",tmplam,0.5*alam);
//                                   tmplam=0.5*alam;
//                               }
//                       }
//       //                if (isnan(tmplam))
//       //                 {
//       //                     printf(" -  tmplam not a number\n");
//       //                     abort();
//       //                 }
//                       alam2=alam;
//                       f2 = SSE_new;
//                       alam=fmax(tmplam,0.1*alam);
//                   }
// 
//               }
// 
//               // Now that we have run through the linsearch, fvec and xb have been updated.
//               // Now test for convergence and update
// 
// 
//               // Test for residual
//               test = 0.0;
//               for (int i=0; i<18; i++) {
//                   if (fabs(fvec[i]) > test) {
//                       test = fabs(fvec[i]);
//                   }
//               }
//               ok = (test < tol_res);
// //      if (ok==true) fprintf(stdout,"res passed: %5.5lf\n",test);
// 
//               // Test for zero gradient (spurious convergence)
//               if (ok==false)
//               {
//                   test=0.0;
//                   den = fmax(SSE_new,0.5*18.0);
//                   for (int i=0; i<18; i++)
//                   {
//                       temp = grad[i]*grad[i];//fabs(grad[i])*fmax(fabs(xb[i]),1.0)/den;
//                       if (temp > test) test = temp;
//                   }
//                   ok = (test<tol_min);
//                   if (ok==true) fprintf(stdout,"grad passed: %5.5lf\n",test);
//               }
// 
// //              // Test for converged x value
//              if (ok==false)
//                  {
//                  test = 0.0;
//                  for (int i=0; i<18; i++) {
//                      if (fabs(xb[i]-xold[i]) > test) {
//                          test = fabs(xb[i]-xold[i]);
//                      }
//                  }
//                      ok =  (test < 1e-24);//tolx);
// //                  if (ok==true) fprintf(stdout,"x passed: %5.5lf\n",test);
//              }
// 
// 
// //      if (alam<0.001) fprintf(stdout,"Ves: %d  Alam: %5.5lf\n",parent,alam);
// //      if (j==1) fprintf(stdout,"\nj=1 p_parent: %5.5lf   p_d1:  %5.5lf  p_d2:  %5.5lf   dpA: %5.5lf   dpB:   %5.5lf\n",PV->P(N,xb[9],WM)*rho*g*Lr/conv,B1->P(0,xb[12],WM)*rho*g*Lr/conv,B2->P(0,xb[15],WM)*rho*g*Lr/conv,delta_PA*rho*g*Lr/conv,delta_PB*rho*g*Lr/conv);
// */
//     j = j+1;
//   }
// 
//   // Solutions is applied, and right boundary is updated.
//   PV->Anew[N] = xb[ 9];
//   PV->Qnew[N] = xb[ 0];
//   B1->Anew[0] = xb[12];
//   B1->Qnew[0] = xb[ 3];
//   B2->Anew[0] = xb[15];
//   B2->Qnew[0] = xb[ 6];
// 
// 
// 
//     if (j >=ntrial)
//     {
// //         error ("arteries.C","Root not found in the stenotic bifurcation");
// //         exit(1);
//        PV->Anew[N] = PV->Aold[N];//xb[ 9];
//        PV->Qnew[N] = PV->Qold[N];//xb[ 0];
//        B1->Anew[0] = B1->Aold[0];//xb[12];
//        B1->Qnew[0] = B1->Qold[0];//xb[ 3];
//        B2->Anew[0] = B2->Aold[0];//xb[15];
//        B2->Qnew[0] = B2->Qold[0];//xb[ 6];
// 
// //         fprintf(stdout,"ves: %d\n",parent);
// //         for (int i=0; i<18; i++) fprintf(stdout,"   %3.12lf    ",fvec[i]);
// //         fprintf(stdout,"xb:\n");
// //         for (int i=0; i<18; i++) fprintf(stdout,"   %3.12lf    ",xb[i]);
// //         fprintf(stdout,"\n p_parent: %5.5lf   p_d1:  %5.5lf  p_d2:  %5.5lf   t1A: %5.5lf  t2A: %5.5lf   t3A: %5.5lf   t1B: %5.5lf  t2B: %5.5lf   t3B: %5.5lf \n",PV->P(N,xb[9],WM)*rho*g*Lr/conv,B1->P(0,xb[12],WM)*rho*g*Lr/conv,B2->P(0,xb[15],WM)*rho*g*Lr/conv,term1A,term2A,term3A,term1B,term2B,term3B);
// //         fprintf(stdout,"\n\n\n");
// 
// //         PV->Anew[N] = xb[ 9];
// //         PV->Qnew[N] = xb[ 0];
// //         B1->Anew[0] = xb[12];
// //         B1->Qnew[0] = xb[ 3];
// //         B2->Anew[0] = xb[15];
// //         B2->Qnew[0] = xb[ 6];
//     }
// }




void bound_trif (double theta, double gamma, Tube *Arteries[], int parent, int D1,int D2, int D3, int WM)
{
    // MJC: Try to define the daughters here
    Tube* PV  = Arteries[ parent];
    Tube* B1  = Arteries[ D1];
    Tube* B2  = Arteries[ D2];
    Tube* B3  = Arteries[ D3];
//    int ND = 3; // Number daughters: COME BACK AND MAKE XB AND FJAC FUNCTIONS OF THIS NUMBER
      int N = PV->N;
      double PN;
      int j = 1;
      int ok = false;
      const int ntrial = 5000;
    
    double k1[4];
    double k2[4];
    double k3[4];
    double k4[4];
    
      // These are the flows at the half time step.
            k1[0] = PV->Qold[N]  + theta*(PV->R2h[N-1]) + gamma*(PV->S2h[N-1]);
          
            k1[1] = B1->Qold[0] - theta*(B1->R2h[0]) + gamma*(B1->S2h[0]);
          
            k1[2] = B2->Qold[0] - theta*(B2->R2h[0]) + gamma*(B2->S2h[0]);
          
            k1[3] = B3->Qold[0] - theta*(B3->R2h[0]) + gamma*(B3->S2h[0]);
          
            // These are areas at the half time step
            k2[0] = PV->Aold[N] + theta*(PV->R1h[N-1]);
            k2[1] = B1->Aold[0] - theta*(B1->R1h[0]);
            k2[2] = B2->Aold[0] - theta*(B2->R1h[0]);
            k2[3] = B3->Aold[0] - theta*(B3->R1h[0]);
          
            // These are flows at the half time stpe
            k3[0] = PV->Qh[N-1]/2.0;
            k3[1] = B1->Qh[0]/2.0;
            k3[2] = B2->Qh[0]/2.0;
            k3[3] = B3->Qh[0]/2.0;
          
            
            k4[0] = PV->Ah[N-1]/2.0;
            k4[1] = B1->Ah[0]/2.0;
            k4[2] = B2->Ah[0]/2.0;
            k4[3] = B3->Ah[0]/2.0;
          
            double xb[6*4];

            // The approximative initial guesses are applied.
              // Initial flows
          xb[ 0] =  PV->Qh[N-1];                        //Initial guess for Qp_xb n+1
          xb[ 1] = (PV->Qold[N-1] + PV->Qold[N])/2.0;   //Initial guess for Qp_xb^n+0.5
          xb[ 2] =  PV->Qold[N];                        //Initial guess for Qp_xb+0.5 n+0.5
          
          xb[ 3] =  B1->Qh[0];                          //Initial guess for Qd1_xb n+1
          xb[ 4] = (B1->Qold[0] + B1->Qold[1])/2.0;     //Initial guess for Qd1_xb n+0.5
          xb[ 5] =  B1->Qold[0];                        //Initial guess for Qd1_xb+0.5 n+0.5
          
          xb[ 6] =  B2->Qh[0];                          //Initial guess for Qd2_xb n+1
          xb[ 7] = (B2->Qold[0] + B2->Qold[1])/2.0;     //Initial guess for Qd2_xb n+0.5
          xb[ 8] =  B2->Qold[0];                        //Initial guess for Qd2_xb+0.5 n+0.5
          
          xb[ 9] =  B3->Qh[0];                          //Initial guess for Qd3_xb n+1
          xb[10] = (B3->Qold[0] + B3->Qold[1])/2.0;     //Initial guess for Qd3_xb n+0.5
          xb[11] =  B3->Qold[0];                        //Initial guess for Qd3_xb+0.5 n+0.5
          
          xb[12] =  PV->Ah[N-1];                        //Initial guess for Ap_xb n+1
          xb[13] = (PV->Aold[N-1] + PV->Aold[N])/2.0;   //Initial guess for Ap_xb^n+0.5
          xb[14] =  PV->Aold[N];                        //Initial guess for Ap_xb+0.5 n+0.5
          
          xb[15] =  B1->Ah[0];                          //Initial guess for Ad1_xb n+1
          xb[16] = (B1->Aold[0] + B1->Aold[1])/2.0;     //Initial guess for Ad1_xb n+0.5
          xb[17] =  B1->Aold[0];                        //Initial guess for Ad1_xb+0.5 n+0.5
          
          xb[18] =  B2->Ah[0];                          //Initial guess for Ad2_xb n+1
          xb[19] = (B2->Aold[0] + B2->Aold[1])/2.0;     //Initial guess for Ad2_xb n+0.5
          xb[20] =  B2->Aold[0];                        //Initial guess for Ad2_xb+0.5 n+0.5

          xb[21] =  B3->Ah[0];                          //Initial guess for Ad3_xb n+1
          xb[22] = (B3->Aold[0] + B3->Aold[1])/2.0;     //Initial guess for Ad3_xb n+0.5
          xb[23] =  B3->Aold[0];                        //Initial guess for Ad3_xb+0.5 n+0.5

          // This is where a Bernoulli term can be perscribed; set to zero otherwise
          double k7[4];
          double k7h[4];
          k7[ 0] = 0.0;
          k7[ 1] = 0.0;
          k7[ 2] = 0.0;
          k7h[0] = 0.0;
          k7h[1] = 0.0;
          k7h[2] = 0.0;

            // The residuals (fvec), and the Jacobian is determined, and if possible
            // the system of equations is solved.
            while (j <= ntrial && ok==false) // Find the zero
            {
              double fvec[24];
                
                
              // The residuals.
                
                // Characteristic Q residual at n+1
              fvec[0]  = k1[0]  - xb[0] -
                      theta*(sq(xb[2])/xb[14] + PV->Bh(N,xb[14],WM)) +
                      gamma*(F(xb[2],xb[14])+PV->dBdx1h(N,xb[14],WM));

              fvec[1]  = k1[1]  - xb[3] +
                      theta*(sq(xb[5])/xb[17] + B1->Bh(-1,xb[17],WM)) +
                      gamma*(F(xb[5],xb[17])  + B1->dBdx1h(-1,xb[17],WM));

              fvec[2]  = k1[2] - xb[6] +
                      theta*(sq(xb[8])/xb[20] + B2->Bh(-1,xb[20],WM)) +
                      gamma*(F(xb[8],xb[20])  + B2->dBdx1h(-1,xb[20],WM));
                
              fvec[3]  = k1[3] - xb[9] +
                      theta*(sq(xb[11])/xb[23] + B3->Bh(-1,xb[23],WM)) +
                      gamma*(F(xb[11],xb[23])  + B3->dBdx1h(-1,xb[23],WM));
                
                
                // Characteristic A residual at n+1
              fvec[4]  = - theta*xb[ 2] - xb[12]  + k2[0];
              fvec[5]  =   theta*xb[ 5] - xb[15]  + k2[1];
              fvec[6]  =   theta*xb[ 8] - xb[18]  + k2[2];
              fvec[7]  =   theta*xb[11] - xb[21]  + k2[3];
                
              // Flow residuals at n+1/2 (ghost points)
              fvec[ 8]  = - xb[ 1] + xb[ 2]/2.0 + k3[0];
              fvec[ 9]  = - xb[ 4] + xb[ 5]/2.0 + k3[1];
              fvec[10]  = - xb[ 7] + xb[ 8]/2.0 + k3[2];
              fvec[11]  = - xb[10] + xb[11]/2.0 + k3[3];
                
              // Area residuals at n+1/2 (ghost points)
              fvec[12] = - xb[13] + xb[14]/2.0 + k4[0];
              fvec[13] = - xb[16] + xb[17]/2.0 + k4[1];
              fvec[14] = - xb[19] + xb[20]/2.0 + k4[2];
              fvec[15] = - xb[22] + xb[23]/2.0 + k4[3];
                
                
              // Flow conservation residuals (n+1/2 and n+1)
              fvec[16] = - xb[ 1] + xb[ 4] + xb[7] + xb[10];
              fvec[17] = - xb[ 0] + xb[ 3] + xb[6] + xb[9];
                
                
              // Pressure continuity at the n+1/2 time step
              PN    = PV->P(N,xb[13],WM);
              double u_n_half = sq(xb[1]/xb[13]);

              // The if statements here only matter if a minor loss
              // is included.
                fvec[18] =  - PN + B1->P(0,xb[16],WM) + ab(k7h[0])*u_n_half;
                fvec[19] =  - PN + B2->P(0,xb[19],WM) + ab(k7h[1])*u_n_half;
                fvec[20] =  - PN + B3->P(0,xb[22],WM) + ab(k7h[2])*u_n_half;
              // Pressure continuity at the n+1 time step
              PN    = PV->P(N,xb[12],WM);
              double u_n_1 = sq(xb[0]/xb[12]);
                fvec[21] = - PN + B1->P(0,xb[15],WM) + ab(k7[0])*u_n_1;
                fvec[22] = - PN + B2->P(0,xb[18],WM) + ab(k7[1])*u_n_1;
                fvec[23] = - PN + B3->P(0,xb[21],WM) + ab(k7[2])*u_n_1;


              for (int row = 0; row < 4*6; row++)
                for (int col = 0; col < 4*6; col++)
                  fjac[row][col] = 0.0;

                
                
      //        // The Jacobian.
                double chi[16];
                
                // Here are the residuals for the characteristic matching for flow
                chi[0] = -2.0*theta*xb[ 2]/xb[14] + gamma*dFdQ(xb[14]);
                chi[2] =  2.0*theta*xb[ 5]/xb[17] + gamma*dFdQ(xb[17]);
                chi[4] =  2.0*theta*xb[ 8]/xb[20] + gamma*dFdQ(xb[20]);
                chi[6] =  2.0*theta*xb[11]/xb[23] + gamma*dFdQ(xb[23]);
                
                
                // Here are the residuals for the area characteristic matching
                chi[1] = theta*(sq(xb[2]/xb[14]) - PV->dBdAh(N,xb[14],WM)) +
                           gamma*(dFdA(xb[2],xb[14]) + PV->d2BdAdxh(N,xb[14],WM));
                
                chi[3] = theta*( -sq(xb[5]/xb[17]) + B1->dBdAh(-1,xb[17],WM)) +
                          gamma*(dFdA(xb[5],xb[17]) + B1->d2BdAdxh(-1,xb[17],WM));
                
                chi[5] = theta*( -sq(xb[8]/xb[20]) + B2->dBdAh(-1,xb[20],WM)) +
                          gamma*(dFdA(xb[8],xb[20]) + B2->d2BdAdxh(-1,xb[20],WM));
                
                chi[7] = theta*( -sq(xb[11]/xb[23]) + B3->dBdAh(-1,xb[23],WM)) +
                          gamma*(dFdA(xb[11],xb[23]) + B3->d2BdAdxh(-1,xb[23],WM));
                
                // Here is pressure conservation (n+1/2)
                chi[8]  = -PV->dPdA(N,xb[13],WM) + sq(xb[1])/cu(xb[13])*(-2.0*ab(k7h[0])); //Loss term
                chi[9]  = B1->dPdA(0,xb[16],WM);
                chi[10] = B2->dPdA(0,xb[19],WM);
                chi[11] = B3->dPdA(0,xb[22],WM);
                
                // Here is pressure conservation (n+1)
                chi[12] = -PV->dPdA(N,xb[12],WM) + sq(xb[1])/cu(xb[12])*(-2.0*ab(k7h[0])); //Loss term
                chi[13] = B1->dPdA(0,xb[15],WM);
                chi[14] = B2->dPdA(0,xb[18],WM);
                chi[15] = B3->dPdA(0,xb[21],WM);
                
                // The jacobian
                // Order is [row][column]
                
                fjac[0][ 0] = -1.0;
                fjac[0][ 2] = chi[0];
                fjac[0][14] = chi[1];
                
                fjac[1][ 3] = -1.0;
                fjac[1][ 5] = chi[2];
                fjac[1][17] = chi[3];
                
                fjac[2][ 6] = -1.0;
                fjac[2][ 8] = chi[4];
                fjac[2][20] = chi[5];
                
                fjac[3][ 9] = -1.0;
                fjac[3][11] = chi[6];
                fjac[3][23] = chi[7];
                
                fjac[4][ 2] = -theta;
                fjac[4][12] = -1.0;
                
                fjac[5][ 5] = theta;
                fjac[5][15] = -1.0;
                
                fjac[6][ 8] = theta;
                fjac[6][18] = -1.0;
                
                fjac[7][11] = theta;
                fjac[7][21] = -1.0;
                 
                fjac[8][ 1] = -1.0;
                fjac[8][ 2] = 0.5;
                
                fjac[9][ 4] = -1.0;
                fjac[9][ 5] = 0.5;
                
                fjac[10][ 7] = -1.0;
                fjac[10][ 8] = 0.5;
                
                fjac[11][10] = -1.0;
                fjac[11][11] = 0.5;
                
                fjac[12][13] = -1.0;
                fjac[12][14] = 0.5;
                
                fjac[13][16] = -1.0;
                fjac[13][17] = 0.5;
                
                fjac[14][19] = -1.0;
                fjac[14][20] = 0.5;
                
                fjac[15][22] = -1.0;
                fjac[15][23] = 0.5;
                
                fjac[16][ 1] = -1.0;
                fjac[16][ 4] = 1.0;
                fjac[16][ 7] = 1.0;
                fjac[16][10] = 1.0;
                
                fjac[17][ 0] = -1.0;
                fjac[17][ 3] = 1.0;
                fjac[17][ 6] = 1.0;
                fjac[17][ 9] = 1.0;
                
                fjac[18][13] = chi[8];
                fjac[18][16] = chi[9];
                
                fjac[19][13] = chi[8];
                fjac[19][19] = chi[10];
                
                fjac[20][13] = chi[8];
                fjac[20][22] = chi[11];
                
                fjac[21][12] = chi[12];
                fjac[21][15] = chi[13];
                
                fjac[22][12] = chi[12];
                fjac[22][18] = chi[14];
                
                fjac[23][12] = chi[12];
                fjac[23][21] = chi[15];
          
          // For debugging, you can print the jacobian
//                  if(j==1){
//                      printf("X\n");
//                      for (int i=0; i<24; i++)
//                      {
//                       printf("%lf ",xb[i]);
//                      }
//                      printf("\n\n");
//                      printf("JACOBIAN\n");
//                      for (int i=0; i<24; i++)
//                      {
//                          for (int j=0; j<24; j++) {
//                              printf("%lf ",fjac[i][j]);
//                          }
//                          printf("\n");
//                      }
//                      printf("\n\n");
//                  }
        
        // Check whether solution is close enough. If not run the loop again.
        int ch = zero (xb, 24, 1.0e-12, 1.0e-12, fvec, fjac);
        if (ch == 1) ok = true;

        j = j+1;
      }

      // Solutions is applied, and right boundary is updated.
      PV->Anew[N] = xb[ 12];
      PV->Qnew[N] = xb[ 0];
      B1->Anew[0] = xb[15];
      B1->Qnew[0] = xb[ 3];
      B2->Anew[0] = xb[18];
      B2->Qnew[0] = xb[ 6];
      B3->Anew[0] = xb[21];
      B3->Qnew[0] = xb[ 9];
    
//    fprintf(stdout,"Flow -  P: %lf D1: %lf D2: %lf D3: %lf\n",xb[0],xb[3],xb[6],xb[9]);
//
        if (j >=ntrial) {error ("arteries.C","Root not found in the trifurcation");
            exit(1);}
}

// FROM GITLAB

void bound_sten_bif (double theta, double gamma, Tube *Arteries[], int parent, int D1, int D2, int WM)
{
    // MJC: Try to define the daughters here
    Tube* PV   = Arteries[ parent];
    Tube* B1  = Arteries[ D1];
    Tube* B2  = Arteries[ D2];
    
    
  int N = PV->N;
  int j = 1;
  int ok = false;
  const int ntrial = 500000;
    
    // Define the stenosis parameters based on derivations
    // last edited by MJC 5/27/2020
        double C,Kt,Ku,Ls;//,Ap,Rs,La,As,Rp,sten_factor;
            // For now, just set values based on literature
            Kt = 1.52; // Dissapation due to turbulent forces
            Ku = 1.2; // Dissapation due to inertial forces
            Ls = PV->sten_length; // Define this outside in future versions
            C = (1.0-PV->sten_factor); // Define degree of stenosis as As/Ap, or 1-sten_factor
//            Ap = PV->Aold[N];//M_PI*sq(PV->rtop); // USING DYNAMIC AREA
//            As = Ap*C;
//            Rp = sqrt(Ap/M_PI);
//            Rs = sqrt(C*Ap/M_PI);
//            La = 0.83*Ls + 3.28*Rs;
//            fprintf(stdout,"Stenosis Parameters: Kt:%lf Ku:%lf Ls:%lf  Sten:%lf\n",Kt,Ku,Ls,C);

        double dp1,dp2,dp3;
        if (C==1.0) {
            dp1=0.0; dp2=0.0; dp3=0.0;
        }
        else
        {
            dp1 = (8.0*Ls*mu*M_PI)/sq(C); // From Karniadakis group 2019
//            dp1 = (16.0*La*mu)/(sq(C)*Rp); //From Shadden's group 2020
            dp2 = rho*Kt*sq(1.0/C - 1.0)/2.0;
            dp3 = rho*Ku*Ls;
        }
//        fprintf(stdout,"dp terms: dp1:%lf dp2:%lf dp3:%lf\n",dp1,dp2,dp3);
    

    double k1[3];
    double k2[3];
    double k3[3];
    double k4[3];
    
  k1[0] = PV->Qold[N]  + theta*(PV->R2h[N-1]) + gamma*(PV->S2h[N-1]);
  k1[1] = B1->Qold[0] - theta*(B1->R2h[0]) + gamma*(B1->S2h[0]);
  k1[2] = B2->Qold[0] - theta*(B2->R2h[0]) + gamma*(B2->S2h[0]);
    

  k2[0] = PV->Aold[N] + theta*(PV->R1h[N-1]);
  k2[1] = B1->Aold[0] - theta*(B1->R1h[0]);
  k2[2] = B2->Aold[0] - theta*(B2->R1h[0]);

  k3[0] = PV->Qh[N-1]/2.0;
  k3[1] = B1->Qh[0]/2.0;
  k3[2] = B2->Qh[0]/2.0;

  k4[0]  = PV->Ah[N-1]/2.0;
  k4[1]  = B1->Ah[0]/2.0;
  k4[2]  = B2->Ah[0]/2.0;

  double xb[6*3];

  // The approximative initial guesses are applied.
  xb[ 0] =  PV->Qh[N-1];                        //Initial guess for Q1_xb n+1
  xb[ 1] = (PV->Qold[N-1] + PV->Qold[N])/2.0;       //Initial guess for Q1_xb^n+0.5
  xb[ 2] =  PV->Qold[N];                        //Initial guess for Q1_xb+0.5 n+0.5
  xb[ 3] =  B1->Qh[0];                      //Initial guess for Q2_xb n+1
  xb[ 4] = (B1->Qold[0] + B1->Qold[1])/2.0; //Initial guess for Q2_xb n+0.5
  xb[ 5] =  B1->Qold[0];                    //Initial guess for Q2_xb+0.5 n+0.5
  xb[ 6] =  B2->Qh[0];                      //Initial guess for Q3_xb n+1
  xb[ 7] = (B2->Qold[0] + B2->Qold[1])/2.0; //Initial guess for Q3_xb n+0.5
  xb[ 8] =  B2->Qold[0];                    //Initial guess for Q3_xb+0.5 n+0.5
  xb[ 9] =  PV->Ah[N-1];                        //Initial guess for A1_xb n+1
  xb[10] = (PV->Aold[N-1] + PV->Aold[N])/2.0;       //Initial guess for A1_xb^n+0.5
  xb[11] =  PV->Aold[N];                        //Initial guess for A1_xb+0.5 n+0.5
  xb[12] =  B1->Ah[0];                      //Initial guess for A2_xb n+1
  xb[13] = (B1->Aold[0] + B1->Aold[1])/2.0; //Initial guess for A2_xb n+0.5
  xb[14] =  B1->Aold[0];                    //Initial guess for A2_xb+0.5 n+0.5
  xb[15] =  B2->Ah[0];                      //Initial guess for A3_xb n+1
  xb[16] = (B2->Aold[0] + B2->Aold[1])/2.0; //Initial guess for A3_xb n+0.5
  xb[17] =  B2->Aold[0];                    //Initial guess for A3_xb+0.5 n+0.5


    // MJC: reduce area because of pressure loss
//    xb[12] = C*xb[12];
//    xb[13] = C*xb[13];
//    xb[14] = C*xb[14];
//    xb[15] = C*xb[15];
//    xb[16] = C*xb[16];
//    xb[17] = C*xb[17];
  // The residuals (fvec), and the Jacobian is determined, and if possible
  // the system of equations is solved.
  while (j <= ntrial && ok==false) // Find the zero
  {
    double fvec[18];
    // The residuals.
      
      // Characteristic Q residual at n+1
    fvec[0]  = k1[0]  - xb[0] -
    theta*(sq(xb[2])/xb[11] + PV->Bh(N,xb[11],WM)) +
    gamma*(F(xb[2],xb[11])+PV->dBdx1h(N,xb[11],WM));

    fvec[1]  = k1[1]  - xb[3] +
                 theta*(sq(xb[5])/xb[14] + B1->Bh(-1,xb[14],WM)) +
                 gamma*(F(xb[5],xb[14])  + B1->dBdx1h(-1,xb[14],WM));

    fvec[2]  = k1[2] - xb[6] +
    theta*(sq(xb[8])/xb[17] + B2->Bh(-1,xb[17],WM)) +
    gamma*(F(xb[8],xb[17])  + B2->dBdx1h(-1,xb[17],WM));

    // Characteristic A residual at n+1
    fvec[3]  = - theta*xb[2] - xb[ 9] + k2[0];
    fvec[4]  =   theta*xb[5] - xb[12] + k2[1];
    fvec[5]  =   theta*xb[8] - xb[15] + k2[2];
      
    // Flow residuals at n+1/2 (ghost points)
    fvec[6]  = - xb[ 1] + xb[ 2]/2.0 + k3[0];
    fvec[7]  = - xb[ 4] + xb[ 5]/2.0 + k3[1];
    fvec[8]  = - xb[ 7] + xb[ 8]/2.0 + k3[2];
      
    // Area residuals at n+1/2 (ghost points)
    fvec[9]  = - xb[10] + xb[11]/2.0 + k4[0];
    fvec[10] = - xb[13] + xb[14]/2.0 + k4[1];
    fvec[11] = - xb[16] + xb[17]/2.0 + k4[2];
      
    // Flow conservation residuals (n+1/2 and n+1)
    fvec[12] = - xb[ 1] + xb[ 4] + xb[7];
    fvec[13] = - xb[ 0] + xb[ 3] + xb[6];

      double delta_PA,delta_PB,term1A,term2A,term3A,term1B,term2B,term3B;
    // Define the pressure loss terms. Order is viscous, turbulent, and inertial
        term1A = dp1 * xb[1]/sq(xb[10]);
      
        term2A = (dp2 / sq(xb[10])) * ab(xb[1])*xb[1]; // Using absolute value
      
        term3A = (dp3/xb[10]) * (xb[1]-(PV->Qold[N]))/(2.0*gamma);
      
        delta_PA =  term1A + term2A + term3A;
        
    
        term1B = dp1 * xb[0]/sq(xb[9]);
      
        term2B = (dp2 / sq(xb[9])) * ab(xb[0])*xb[0]; // Using absolute value
      
        term3B = (dp3/xb[9]) * (xb[0]-(PV->Qold[N]))/(2.0*gamma);
      
        delta_PB =  term1B + term2B + term3B;

    // Pressure loss at the n+1/2 time step
      fvec[14] = - PV->P(N,xb[ 10],WM) + B1->P(0,xb[13],WM) + delta_PA;
      fvec[15] = - PV->P(N,xb[ 10],WM) + B2->P(0,xb[16],WM) + delta_PA;

    // Pressure continuity at the n+1 time step
    fvec[16] = - PV->P(N,xb[ 9],WM) + B1->P(0,xb[12],WM) + delta_PB;
    fvec[17] = - PV->P(N,xb[ 9],WM) + B2->P(0,xb[15],WM) + delta_PB;

    for (int row = 0; row < 18; row++)
      for (int col = 0; col < 18; col++)
        fjac[row][col] = 0.0;

      
      
    // The Jacobian.
      
      double chi[14];
      
       // Here are the residuals for the characteristic matching for flow
      chi[0] = -2.0*theta*xb[ 2]/xb[11] + gamma*dFdQ(xb[11]);
      chi[2] =  2.0*theta*xb[ 5]/xb[14] + gamma*dFdQ(xb[14]);
      chi[4] =  2.0*theta*xb[ 8]/xb[17] + gamma*dFdQ(xb[17]);
      
      // Here are the residuals for the area characteristic matching
      chi[1] = theta*(sq(xb[2]/xb[11]) - PV->dBdAh(N,xb[11],WM)) +
                 gamma*(dFdA(xb[2],xb[11]) + PV->d2BdAdxh(N,xb[11],WM));
      
      chi[3] = theta*( -sq(xb[5]/xb[14]) + B1->dBdAh(-1,xb[14],WM)) +
                gamma*(dFdA(xb[5],xb[14]) + B1->d2BdAdxh(-1,xb[14],WM));
      
      chi[5] = theta*( -sq(xb[8]/xb[17]) + B2->dBdAh(-1,xb[17],WM)) +
                gamma*(dFdA(xb[8],xb[17]) + B2->d2BdAdxh(-1,xb[17],WM));
      
      // Here are the pressure loss terms at n+1/2
      // Derivative with respect to flow
      chi[6]  = dp1/sq(xb[10]) + (dp2/sq(xb[10]))*2.0*xb[1]/ab(xb[1])
                           + dp3/(2.0*gamma*xb[10]);
      
      // Derivative with respect to area
      chi[7]  = -PV->dPdA(N,xb[10],WM) - 2.0*dp1*xb[1]/cu(xb[10])
      - 2.0*dp2*ab(xb[1])*xb[1]/cu(xb[10]) - (dp3/sq(xb[10])) * (xb[1]-(PV->Qold[N]))/(2.0*gamma);
      
      chi[8]  = B1->dPdA(0,xb[13],WM);
      
      chi[9]  = B2->dPdA(0,xb[16],WM);
      
      
      // Here are the pressure loss terms at n+1
      // Derivative with respect to flow
      chi[10]  = dp1/sq(xb[9]) + (dp2/sq(xb[9]))*2.0*xb[0]/ab(xb[0])
                           + dp3/(2*gamma*xb[9]);
      
      // Derivative with respect to area
      chi[11]  = -PV->dPdA(N,xb[9],WM) - 2.0*dp1*xb[0]/cu(xb[9])
      - 2.0*dp2*ab(xb[0])*xb[0]/cu(xb[9]) - (dp3/sq(xb[9])) * (xb[0]-(PV->Qold[N]))/(2.0*gamma);
      
      chi[12] = B1->dPdA(0,xb[12],WM);
      
      chi[13] = B2->dPdA(0,xb[15],WM);

      
      //NEW JACOBIAN
            // Order is [row][column]
            fjac[ 0][ 0]  = -1.0;
            fjac[ 0][ 2] = chi[0];
            fjac[ 0][11] = chi[1];
            
            fjac[ 1][ 3] = -1.0;
            fjac[ 1][ 5] = chi[2];
            fjac[ 1][14] = chi[3];

            fjac[ 2][ 6] = -1.0;
            fjac[ 2][ 8] = chi[4];
            fjac[ 2][17] = chi[5];

            fjac[ 3][ 2] = -theta;
            fjac[ 3][ 9] = -1.0;
            
            fjac[ 4][ 5] = theta;
            fjac[ 4][12] = -1.0;
            
            fjac[ 5][ 8] = theta;
            fjac[ 5][15] = -1.0;
            
            fjac[ 6][ 1] = -1.0;
            fjac[ 6][ 2] = 0.5;
            
            fjac[ 7][ 4] = -1.0;
            fjac[ 7][ 5] = 0.5;
            
            fjac[ 8][ 7] = -1.0;
            fjac[ 8][ 8] = 0.5;
            
            fjac[ 9][ 10] = -1.0;
            fjac[ 9][ 11] = 0.5;
            
            fjac[10][13] = -1.0;
            fjac[10][14] = 0.5;
            
            fjac[11][16] = -1.0;
            fjac[11][17] = 0.5;
            
            fjac[12][ 1] = -1.0;
            fjac[12][ 4] = 1.0;
            fjac[12][ 7] = 1.0;
            
            fjac[13][ 0] = -1.0;
            fjac[13][ 3] = 1.0;
            fjac[13][ 6] = 1.0;
            
            fjac[14][ 1] = chi[6];
            fjac[14][10] = chi[7];
            fjac[14][13] = chi[8];
                                    
                                    
            fjac[15][ 1] = chi[6];
            fjac[15][10] = chi[7];
            fjac[15][16] = chi[9];
            
            fjac[16][ 0] = chi[10];
            fjac[16][ 9] = chi[11];
            fjac[16][12] = chi[12];
            
            
            fjac[17][ 0] = chi[10];
            fjac[17][ 9] = chi[11];
            fjac[17][15] = chi[13];
      
      //           For debugging, you can print the jacobian
//      if(j==1){
//          printf("PRESSURES IN %d\n",parent);
//          printf("Parent: %lf\n",PV->P(N,xb[ 9],WM));
//          printf("Parent: %lf\n",B1->P(0,xb[ 12],WM));
//          printf("Parent: %lf\n\n\n",B2->P(0,xb[ 15],WM));
//          printf("X\n");
//          for (int i=0; i<18; i++)
//          {
//           printf("%lf ",xb[i]);
//          }
//          printf("\n\n RESIDUAL\n");
//          for (int i=0; i<18; i++)
//          {
//           printf("%lf ",fvec[i]);
//          }
//          printf("\n\n");
//          printf("JACOBIAN\n");
//          for (int i=0; i<18; i++)
//          {
//              for (int j=0; j<18; j++) {
//                  printf("%lf ",fjac[i][j]);
//              }
//              printf("\n");
//          }
//          printf("\n\n");
//      }

    // Check whether solution is close enough. If not run the loop again.
    int ch = zero (xb, 18, 1.0e-6, 1.0e-6, fvec, fjac);
    if (ch == 1) ok = true;

    j = j+1;
//      if (j>ntrial-50 || j<5){
//      fprintf(stdout,"IN BIF:ITERATION j=%d Parent Qold: %lf Daughter 1 Qold: %lf Daughter 2 Qold: %lf\n",j,PV->Qold[N],B1->Qold[0],B2->Qold[0]);
//      fprintf(stdout,"Pressures: P:%lf D1: %lf  D2: %lf\n",PV->P(N,xb[ 9],WM),B1->P(0,xb[12],WM),B2->P(0,xb[15],WM));
//      fprintf(stdout,"Pressure drop: A:%lf B: %lf\n",delta_PA,delta_PB);
//      fprintf(stdout,"A terms: %lf %lf %lf\n",term1A,term2A,term3A);
//      fprintf(stdout,"B terms: %lf %lf %lf\n",term1B,term2B,term3B);
//      }
  }

  // Solutions is applied, and right boundary is updated.
  PV->Anew[N] = xb[ 9];
  PV->Qnew[N] = xb[ 0];
  B1->Anew[0] = xb[12];
  B1->Qnew[0] = xb[ 3];
  B2->Anew[0] = xb[15];
  B2->Qnew[0] = xb[ 6];
    
    

    if (j >=ntrial) 
    {
        error ("arteries.C","Root not found in the stenotic bifurcation");
//         exit(1);
        PV->Anew[N] = PV->Aold[N];//xb[ 9];
        PV->Qnew[N] = PV->Qold[N];//xb[ 0];
        B1->Anew[0] = B1->Aold[0];//xb[12];
        B1->Qnew[0] = B1->Qold[0];//xb[ 3];
        B2->Anew[0] = B2->Aold[0];//xb[15];
        B2->Qnew[0] = B2->Qold[0];//xb[ 6];
        
//         PV->Anew[N] = xb[ 9];
//         PV->Qnew[N] = xb[ 0];
//         B1->Anew[0] = xb[12];
//         B1->Qnew[0] = xb[ 3];
//         B2->Anew[0] = xb[15];
//         B2->Qnew[0] = xb[ 6];
    }
}



//// Original bound_sten without linesearch (From not working version
// void bound_sten_bif (double theta, double gamma, Tube *Arteries[], int parent, int D1, int D2, int WM)
// {
//    // MJC: Try to define the daughters here
//    Tube* PV   = Arteries[ parent];
//    Tube* B1  = Arteries[ D1];
//    Tube* B2  = Arteries[ D2];
// 
// 
//  int N = PV->N;
//  int j = 1;
//  int ok = false;
//  const int ntrial = 500000;
// 
//    // Define the stenosis parameters based on derivations
//    // last edited by MJC 5/27/2020
//        double C,Kt,Ku,Ls;//,Ap,Rs,La,As,Rp,sten_factor;
//            // For now, just set values based on literature
//            Kt = 1.52; // Dissapation due to turbulent forces
//            Ku = 1.2; // Dissapation due to inertial forces
//            Ls = PV->sten_length; // Define this outside in future versions
//            C = (1.0-PV->sten_factor); // Define degree of stenosis as As/Ap, or 1-sten_factor
// //            Ap = PV->Aold[N];//M_PI*sq(PV->rtop); // USING DYNAMIC AREA
// //            As = Ap*C;
// //            Rp = sqrt(Ap/M_PI);
// //            Rs = sqrt(C*Ap/M_PI);
// //            La = 0.83*Ls + 3.28*Rs;
// //            fprintf(stdout,"Stenosis Parameters: Kt:%lf Ku:%lf Ls:%lf  Sten:%lf\n",Kt,Ku,Ls,C);
// 
//        double dp1,dp2,dp3;
//        if (C==1.0) {
//            dp1=0.0; dp2=0.0; dp3=0.0;
//        }
//        else
//        {
//            dp1 = (8.0*Ls*mu*M_PI)/sq(C); // From Karniadakis group 2019
// //            dp1 = (16.0*La*mu)/(sq(C)*Rp); //From Shadden's group 2020
//            dp2 = rho*Kt*sq(1.0/C - 1.0)/2.0;
//            dp3 = rho*Ku*Ls;
//        }
// //        fprintf(stdout,"dp terms: dp1:%lf dp2:%lf dp3:%lf\n",dp1,dp2,dp3);
// 
// 
//    double k1[3];
//    double k2[3];
//    double k3[3];
//    double k4[3];
// 
//  k1[0] = PV->Qold[N]  + theta*(PV->R2h[N-1]) + gamma*(PV->S2h[N-1]);
//  k1[1] = B1->Qold[0] - theta*(B1->R2h[0]) + gamma*(B1->S2h[0]);
//  k1[2] = B2->Qold[0] - theta*(B2->R2h[0]) + gamma*(B2->S2h[0]);
// 
// 
//  k2[0] = PV->Aold[N] + theta*(PV->R1h[N-1]);
//  k2[1] = B1->Aold[0] - theta*(B1->R1h[0]);
//  k2[2] = B2->Aold[0] - theta*(B2->R1h[0]);
// 
//  k3[0] = PV->Qh[N-1]/2.0;
//  k3[1] = B1->Qh[0]/2.0;
//  k3[2] = B2->Qh[0]/2.0;
// 
//  k4[0]  = PV->Ah[N-1]/2.0;
//  k4[1]  = B1->Ah[0]/2.0;
//  k4[2]  = B2->Ah[0]/2.0;
// 
//  double xb[6*3];
// 
//  // The approximative initial guesses are applied.
//  xb[ 0] =  PV->Qh[N-1];                        //Initial guess for Q1_xb n+1
//  xb[ 1] = (PV->Qold[N-1] + PV->Qold[N])/2.0;       //Initial guess for Q1_xb^n+0.5
//  xb[ 2] =  PV->Qold[N];                        //Initial guess for Q1_xb+0.5 n+0.5
//  xb[ 3] =  B1->Qh[0];                      //Initial guess for Q2_xb n+1
//  xb[ 4] = (B1->Qold[0] + B1->Qold[1])/2.0; //Initial guess for Q2_xb n+0.5
//  xb[ 5] =  B1->Qold[0];                    //Initial guess for Q2_xb+0.5 n+0.5
//  xb[ 6] =  B2->Qh[0];                      //Initial guess for Q3_xb n+1
//  xb[ 7] = (B2->Qold[0] + B2->Qold[1])/2.0; //Initial guess for Q3_xb n+0.5
//  xb[ 8] =  B2->Qold[0];                    //Initial guess for Q3_xb+0.5 n+0.5
//  xb[ 9] =  PV->Ah[N-1];                        //Initial guess for A1_xb n+1
//  xb[10] = (PV->Aold[N-1] + PV->Aold[N])/2.0;       //Initial guess for A1_xb^n+0.5
//  xb[11] =  PV->Aold[N];                        //Initial guess for A1_xb+0.5 n+0.5
//  xb[12] =  B1->Ah[0];                      //Initial guess for A2_xb n+1
//  xb[13] = (B1->Aold[0] + B1->Aold[1])/2.0; //Initial guess for A2_xb n+0.5
//  xb[14] =  B1->Aold[0];                    //Initial guess for A2_xb+0.5 n+0.5
//  xb[15] =  B2->Ah[0];                      //Initial guess for A3_xb n+1
//  xb[16] = (B2->Aold[0] + B2->Aold[1])/2.0; //Initial guess for A3_xb n+0.5
//  xb[17] =  B2->Aold[0];                    //Initial guess for A3_xb+0.5 n+0.5
// 
// 
//    // MJC: reduce area because of pressure loss
// //    xb[12] = C*xb[12];
// //    xb[13] = C*xb[13];
// //    xb[14] = C*xb[14];
// //    xb[15] = C*xb[15];
// //    xb[16] = C*xb[16];
// //    xb[17] = C*xb[17];
//  // The residuals (fvec), and the Jacobian is determined, and if possible
//  // the system of equations is solved.
//  while (j <= ntrial && ok==false) // Find the zero
//  {
//    double fvec[18];
//    // The residuals.
// 
//      // Characteristic Q residual at n+1
//    fvec[0]  = k1[0]  - xb[0] -
//    theta*(sq(xb[2])/xb[11] + PV->Bh(N,xb[11],WM)) +
//    gamma*(F(xb[2],xb[11])+PV->dBdx1h(N,xb[11],WM));
// 
//    fvec[1]  = k1[1]  - xb[3] +
//                 theta*(sq(xb[5])/xb[14] + B1->Bh(-1,xb[14],WM)) +
//                 gamma*(F(xb[5],xb[14])  + B1->dBdx1h(-1,xb[14],WM));
// 
//    fvec[2]  = k1[2] - xb[6] +
//    theta*(sq(xb[8])/xb[17] + B2->Bh(-1,xb[17],WM)) +
//    gamma*(F(xb[8],xb[17])  + B2->dBdx1h(-1,xb[17],WM));
// 
//    // Characteristic A residual at n+1
//    fvec[3]  = - theta*xb[2] - xb[ 9] + k2[0];
//    fvec[4]  =   theta*xb[5] - xb[12] + k2[1];
//    fvec[5]  =   theta*xb[8] - xb[15] + k2[2];
// 
//    // Flow residuals at n+1/2 (ghost points)
//    fvec[6]  = - xb[ 1] + xb[ 2]/2.0 + k3[0];
//    fvec[7]  = - xb[ 4] + xb[ 5]/2.0 + k3[1];
//    fvec[8]  = - xb[ 7] + xb[ 8]/2.0 + k3[2];
// 
//    // Area residuals at n+1/2 (ghost points)
//    fvec[9]  = - xb[10] + xb[11]/2.0 + k4[0];
//    fvec[10] = - xb[13] + xb[14]/2.0 + k4[1];
//    fvec[11] = - xb[16] + xb[17]/2.0 + k4[2];
// 
//    // Flow conservation residuals (n+1/2 and n+1)
//    fvec[12] = - xb[ 1] + xb[ 4] + xb[7];
//    fvec[13] = - xb[ 0] + xb[ 3] + xb[6];
// 
//      double delta_PA,delta_PB,term1A,term2A,term3A,term1B,term2B,term3B;
//    // Define the pressure loss terms. Order is viscous, turbulent, and inertial
//        term1A = dp1 * xb[1]/sq(xb[10]);
// 
//        term2A = (dp2 / sq(xb[10])) * ab(xb[1])*xb[1]; // Using absolute value
// 
//        term3A = (dp3/xb[10]) * (xb[1]-(PV->Qold[N]))/(2.0*gamma);
// 
//        delta_PA =  term1A + term2A + term3A;
// 
// 
//        term1B = dp1 * xb[0]/sq(xb[9]);
// 
//        term2B = (dp2 / sq(xb[9])) * ab(xb[0])*xb[0]; // Using absolute value
// 
//        term3B = (dp3/xb[9]) * (xb[0]-(PV->Qold[N]))/(2.0*gamma);
// 
//        delta_PB =  term1B + term2B + term3B;
// 
//    // Pressure loss at the n+1/2 time step
//      fvec[14] = - PV->P(N,xb[ 10],WM) + B1->P(0,xb[13],WM) + delta_PA;
//      fvec[15] = - PV->P(N,xb[ 10],WM) + B2->P(0,xb[16],WM) + delta_PA;
// 
//    // Pressure continuity at the n+1 time step
//    fvec[16] = - PV->P(N,xb[ 9],WM) + B1->P(0,xb[12],WM) + delta_PB;
//    fvec[17] = - PV->P(N,xb[ 9],WM) + B2->P(0,xb[15],WM) + delta_PB;
// 
//    for (int row = 0; row < 18; row++)
//      for (int col = 0; col < 18; col++)
//        fjac[row][col] = 0.0;
// 
// 
// 
//    // The Jacobian.
// 
//      double chi[14];
// 
//       // Here are the residuals for the characteristic matching for flow
//      chi[0] = -2.0*theta*xb[ 2]/xb[11] + gamma*dFdQ(xb[11]);
//      chi[2] =  2.0*theta*xb[ 5]/xb[14] + gamma*dFdQ(xb[14]);
//      chi[4] =  2.0*theta*xb[ 8]/xb[17] + gamma*dFdQ(xb[17]);
// 
//      // Here are the residuals for the area characteristic matching
//      chi[1] = theta*(sq(xb[2]/xb[11]) - PV->dBdAh(N,xb[11],WM)) +
//                 gamma*(dFdA(xb[2],xb[11]) + PV->d2BdAdxh(N,xb[11],WM));
// 
//      chi[3] = theta*( -sq(xb[5]/xb[14]) + B1->dBdAh(-1,xb[14],WM)) +
//                gamma*(dFdA(xb[5],xb[14]) + B1->d2BdAdxh(-1,xb[14],WM));
// 
//      chi[5] = theta*( -sq(xb[8]/xb[17]) + B2->dBdAh(-1,xb[17],WM)) +
//                gamma*(dFdA(xb[8],xb[17]) + B2->d2BdAdxh(-1,xb[17],WM));
// 
//      // Here are the pressure loss terms at n+1/2
//      // Derivative with respect to flow
//      chi[6]  = dp1/sq(xb[10]) + (dp2/sq(xb[10]))*2.0*xb[1]/ab(xb[1])
//                           + dp3/(2.0*gamma*xb[10]);
// 
//      // Derivative with respect to area
//      chi[7]  = -PV->dPdA(N,xb[10],WM) - 2.0*dp1*xb[1]/cu(xb[10])
//      - 2.0*dp2*ab(xb[1])*xb[1]/cu(xb[10]) - (dp3/sq(xb[10])) * (xb[1]-(PV->Qold[N]))/(2.0*gamma);
// 
//      chi[8]  = B1->dPdA(0,xb[13],WM);
// 
//      chi[9]  = B2->dPdA(0,xb[16],WM);
// 
// 
//      // Here are the pressure loss terms at n+1
//      // Derivative with respect to flow
//      chi[10]  = dp1/sq(xb[9]) + (dp2/sq(xb[9]))*2.0*xb[0]/ab(xb[0])
//                           + dp3/(2*gamma*xb[9]);
// 
//      // Derivative with respect to area
//      chi[11]  = -PV->dPdA(N,xb[9],WM) - 2.0*dp1*xb[0]/cu(xb[9])
//      - 2.0*dp2*ab(xb[0])*xb[0]/cu(xb[9]) - (dp3/sq(xb[9])) * (xb[0]-(PV->Qold[N]))/(2.0*gamma);
// 
//      chi[12] = B1->dPdA(0,xb[12],WM);
// 
//      chi[13] = B2->dPdA(0,xb[15],WM);
// 
// 
//      //NEW JACOBIAN
//            // Order is [row][column]
//            fjac[ 0][ 0]  = -1.0;
//            fjac[ 0][ 2] = chi[0];
//            fjac[ 0][11] = chi[1];
// 
//            fjac[ 1][ 3] = -1.0;
//            fjac[ 1][ 5] = chi[2];
//            fjac[ 1][14] = chi[3];
// 
//            fjac[ 2][ 6] = -1.0;
//            fjac[ 2][ 8] = chi[4];
//            fjac[ 2][17] = chi[5];
// 
//            fjac[ 3][ 2] = -theta;
//            fjac[ 3][ 9] = -1.0;
// 
//            fjac[ 4][ 5] = theta;
//            fjac[ 4][12] = -1.0;
// 
//            fjac[ 5][ 8] = theta;
//            fjac[ 5][15] = -1.0;
// 
//            fjac[ 6][ 1] = -1.0;
//            fjac[ 6][ 2] = 0.5;
// 
//            fjac[ 7][ 4] = -1.0;
//            fjac[ 7][ 5] = 0.5;
// 
//            fjac[ 8][ 7] = -1.0;
//            fjac[ 8][ 8] = 0.5;
// 
//            fjac[ 9][ 10] = -1.0;
//            fjac[ 9][ 11] = 0.5;
// 
//            fjac[10][13] = -1.0;
//            fjac[10][14] = 0.5;
// 
//            fjac[11][16] = -1.0;
//            fjac[11][17] = 0.5;
// 
//            fjac[12][ 1] = -1.0;
//            fjac[12][ 4] = 1.0;
//            fjac[12][ 7] = 1.0;
// 
//            fjac[13][ 0] = -1.0;
//            fjac[13][ 3] = 1.0;
//            fjac[13][ 6] = 1.0;
// 
//            fjac[14][ 1] = chi[6];
//            fjac[14][10] = chi[7];
//            fjac[14][13] = chi[8];
// 
// 
//            fjac[15][ 1] = chi[6];
//            fjac[15][10] = chi[7];
//            fjac[15][16] = chi[9];
// 
//            fjac[16][ 0] = chi[10];
//            fjac[16][ 9] = chi[11];
//            fjac[16][12] = chi[12];
// 
// 
//            fjac[17][ 0] = chi[10];
//            fjac[17][ 9] = chi[11];
//            fjac[17][15] = chi[13];
// 
// 
//    // Check whether solution is close enough. If not run the loop again.
//    int ch = zero (xb, 18, 1.0e-6, 1.0e-6, fvec, fjac);
//    if (ch == 1) ok = true;
// 
//    j = j+1;
// //      if (j>ntrial-50 || j<5){
// //      fprintf(stdout,"IN BIF:ITERATION j=%d Parent Qold: %lf Daughter 1 Qold: %lf Daughter 2 Qold: %lf\n",j,PV->Qold[N],B1->Qold[0],B2->Qold[0]);
// //      fprintf(stdout,"Pressures: P:%lf D1: %lf  D2: %lf\n",PV->P(N,xb[ 9],WM),B1->P(0,xb[12],WM),B2->P(0,xb[15],WM));
// //      fprintf(stdout,"Pressure drop: A:%lf B: %lf\n",delta_PA,delta_PB);
// //      fprintf(stdout,"A terms: %lf %lf %lf\n",term1A,term2A,term3A);
// //      fprintf(stdout,"B terms: %lf %lf %lf\n",term1B,term2B,term3B);
// //      }
//  }
// 
//  // Solutions is applied, and right boundary is updated.
//  PV->Anew[N] = xb[ 9];
//  PV->Qnew[N] = xb[ 0];
//  B1->Anew[0] = xb[12];
//  B1->Qnew[0] = xb[ 3];
//  B2->Anew[0] = xb[15];
//  B2->Qnew[0] = xb[ 6];
// 
// 
// 
//    if (j >=ntrial)
//    {
//        error ("arteries.C","Root not found in the stenotic bifurcation");
// //         exit(1);
// //        PV->Anew[N] = PV->Aold[N];//xb[ 9];
// //        PV->Qnew[N] = PV->Qold[N];//xb[ 0];
// //        B1->Anew[0] = B1->Aold[0];//xb[12];
// //        B1->Qnew[0] = B1->Qold[0];//xb[ 3];
// //        B2->Anew[0] = B2->Aold[0];//xb[15];
// //        B2->Qnew[0] = B2->Qold[0];//xb[ 6];
// 
// //         PV->Anew[N] = xb[ 9];
// //         PV->Qnew[N] = xb[ 0];
// //         B1->Anew[0] = xb[12];
// //         B1->Qnew[0] = xb[ 3];
// //         B2->Anew[0] = xb[15];
// //         B2->Qnew[0] = xb[ 6];
//    }
// }

