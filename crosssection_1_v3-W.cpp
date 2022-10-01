//
//  crosssection.cpp
//
//
//  Created by ams52admin on 9/25/18.
//

// Program to compute the diffractive cross section from the
// formula in the 0207027 paper by Enberg Motyka Poludniowski
// using CT14NLO pdf set
//
//

#include "/home/mike/thesis/LHAPDF-6.5.3/include/LHAPDF/LHAPDF.h"
#include "/home/mike/thesis/specialfunctions.h"
#include "/home/mike/thesis/integration.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include "/home/mike/thesis/TGreenFun.hpp"
#include "/home/mike/thesis/tbfk1.hpp"
#include "/home/mike/thesis/stdafx.h"

#include "/home/mike/thesis/statistics.h"
#include "/home/mike/thesis/solvers.h"
#include "/home/mike/thesis/optimization.h"
#include "/home/mike/thesis/linalg.h"
#include "/home/mike/thesis/integration.h"
#include "/home/mike/thesis/fasttransforms.h"
#include "/home/mike/thesis/diffequations.h"
#include "/home/mike/thesis/dataanalysis.h"
#include "/home/mike/thesis/alglibmisc.h"
#include "/home/mike/thesis/alglibinternal.h"
#include "/home/mike/thesis/stdafx.h"
#include "/home/mike/thesis/ap.h"

double ampintegrand(double lk);
double phi0gammaV(double k2,double q2);
double k2tolk(double k2);
double ltok2(double lk);
double qgauss96(double(*func)(double),  double a,  double b);
double PhotonFlux(double W, double sc);
double qcont, gcont;
double quarkcont, gluoncont;

//using namespace alglib;
using namespace LHAPDF;
using namespace std;

double q2sam;
double rap_int;
double tau_int;
const int NTG = 50;
const int NYG = 100;
const int NKG = 150;
//const double YMIN    = 0.0;
//const double YMAX    = 20.0;
const double KMIN2   = 1.0e-4;
const double KMAX2   = 1.0e+10;
const double TAUMIN  = 0.01;
const double TAUMAX  = 5.0; // this is the size in tau for the grid of non-forward BFKL
const double TMIN = 1.;
const double TMAX = 2.;
double L = 10000000.;
double events;
double aux;

TGreenFun nfgluon(NTG, NYG, NKG, TAUMIN,TAUMAX, YMIN, YMAX, KMIN2, KMAX2);
alglib::spline3dinterpolant sglob;


int main(int argc, char* argv[]) {

    if (argc < 3) {
        cerr << "You must specify a PDF set and member number" << endl;
        return 1;
    }
    const string setname = argv[1];
    const string smem = argv[2];
    const int imem = lexical_cast<int>(smem);
    const PDF* pdf = mkPDF(setname, imem);
    vector<int> pids = pdf->flavors();
    const double MINLOGX = -10;
    const double MAXLOGX =   0;
    const double DX = 0.01;
    const int NX = (int) floor((MAXLOGX - MINLOGX)/DX) + 1;
    const double MINLOGQ2 = 1;
    const double MAXLOGQ2 = 8;
    const double DQ2 = 0.01;
    const int NQ2 = (int) floor((MAXLOGQ2 - MINLOGQ2)/DQ2) + 1;
    double xa = 0.0001; double qa2 = 10.0;



    const int precision = 12;  // set precision for the file output
    const int fieldwidth = 22; // set field width for the file output


    // initialize the matrix
    ifstream file_in;

    // Here need to input the correct path to read in the grid
    //file_in.open("/Users/ams52admin/Dropbox/Work/jpsi_loops/current/grids/gf_t_5_km2x5_y20_asb0.1_s0.5_sm_LG.dat",ios::in);
    file_in.open("/Users/douglashowell/Diffraction/grids_nfbfkl/gf_t_5_km2x5_y20_asb0.1_s0.5_sm_LG.dat",ios::in);
    cout << "Starting to read the file... \n" ;


    double *tautemp,*ytemp,*ktemp;
    tautemp = new double[NTG+1];
    ytemp = new double[NYG+1];
    ktemp = new double[NKG+1];
    double *gluholder;
    gluholder = new double[(NTG+1)*(NYG+1)*(NKG+1)];

    if (file_in.is_open())
    {
        double tautemp1,ytemp1,ktemp1;

        for (int i1=0; i1<NTG+1; i1++) {
            for (int i2=0; i2<NYG+1; i2++) {
                for (int i3=0; i3<NKG+1; i3++) {

                    file_in >> tautemp1 ;
                    if (i2==0&&i3==0) {
                        tautemp[i1] = tautemp1;
                    }
                    file_in >> ytemp1 ;
                    if (i1==0&&i3==0) {
                        ytemp[i2] = ytemp1;
                    }
                    file_in >> ktemp1 ;
                    if (i1==0&&i2==0) {
                        ktemp[i3] = ktemp1;
                    }
                    file_in >> nfgluon.GF[i1][i2][i3];
                    gluholder[(NTG+1)*((NYG+1)*i3+i2)+i1] = nfgluon.GF[i1][i2][i3];

                }
            }
        }
    }
    else{cout << "Error opening file " <<endl;}

    cout << "Read file sucessfully... \n";
    file_in.close();

    for (int i1=0; i1<NTG+1; i1++) {
        cout << tautemp[i1] << endl;
    }
    for (int i1=0; i1<NYG+1; i1++) {
        cout << ytemp[i1] << endl;
    }
    for (int i1=0; i1<NKG+1; i1++) {
        cout << ktemp[i1] << endl;
    }


    cout << "***********************************\n";
    // lglib interpolation

    alglib::real_1d_array xx;
    xx.setcontent(NTG+1,tautemp);

    alglib::real_1d_array yy;
    yy.setcontent(NYG+1,ytemp);

    alglib::real_1d_array zz;
    zz.setcontent(NKG+1,ktemp);

    alglib::real_1d_array ff;
    ff.setcontent((NTG+1)*(NYG+1)*(NKG+1),gluholder);

    alglib::spline3dbuildtrilinearv(xx, NTG+1, yy, NYG+1, zz, NKG+1, ff, 1, sglob);

    delete [] gluholder;
    delete [] tautemp;
    delete [] ktemp;
    delete [] ytemp;

    // sample values of tau, rapidity and k2 for which we interpolate
    tau_int = 4.0;
    rap_int = 3.0;
    double k2_int  = 10000.0;

    cout << "\ntau: " << tau_int << "   Y: " << rap_int << "   k2: " << k2_int << " GeV^2 ";
    cout << "\nInterpolated value of G: " << nfgluon.interpolate(tau_int,rap_int,k2_int) << endl;


    //****************************************this is irrelevant
/*    tau_int   = 1.0;
    q2sam = tau_int * MassJPsi * MassJPsi;
    double z = 0.8;
    rap_int   = 2.0 * z / alphas_fixed_grid ;
     double result = qgauss96(ampintegrand,LKMIN,k2tolk(q2sam/4.0)) + qgauss96(ampintegrand,k2tolk(q2sam/4.0),LKMAX);
     result *= 9.0 * pow(q2sam,2.0)/(16.0 * pi); // normalization in paper by LM et al
     result = result / (Ccal * alphas_fixed_grid * alphas_fixed ) * 8.0 / 9.0 * pi ; // 8/9 pi is needed to fix norm.

    cout << " z: " << z << " tau: " << tau_int << " F : " << result << endl;*/
    //  exit(0);
// ****************************************

    ofstream fout;
    fout//.open("/home/michal/jpsi_nonforward_eic/grids_nfbfkl/current/xsec_dt_cteq14nlo_xfix0.01_W100_0.17_asrun_01312020.dat");
            .open("wtest.txt");


    int NY = 200;
    int NTAU = 200;
    double XMIN = 0.01; // the x_min cutoff
    double XMAX = 0.05;
    const double fac=0.5;

/*    double XMID=sqrt(XMIN*XMAX);
    XMIN=pow(XMIN/XMID,fac)*XMID;
    XMAX=pow(XMAX/XMID,fac)*XMID;*/

//for(double XMIN=0.01; XMIN<0.21; XMIN+=0.01)
//  {
    double dely = (log(XMAX/XMIN)) / NY;

    double dsigmahatdt = 0.0;
    double dsigmadt = 0.0;
    const double MJPsi2=MassJPsi * MassJPsi;

    //double deltau = (log(TAUMAX)-log(TAUMIN)) / NTAU;
    double delt = (TMAX - TMIN) / NTAU;
//    double W2; // the energy squared
//    double W;

//double rap_int;

    double W;
    double sc=140.;
    double Ep=275.;
//double rapy;

//for(double dpy=2.; dpy<=6.; dpy+=0.01)
//{

    double dpy=2.;

// W loop
    for(double lW=log(50.); lW<log(sc); lW+=log(sc/50.)/50.)
    {
// Logarithmic grid for W
        double W=exp(lW);
        double W2=W*W;

        dsigmadt = 0;

        // loop over t-momentum transfer
        for (int j=0; j < NTAU; j++)
        {

            //dsigmadt = 0.0;
            aux = 0.0;
            //        double tau = TAUMIN * exp(deltau * j);
            //        tau_int = tau;
            //        q2sam = tau * MassJPsi * MassJPsi;

// Here we get the t value
            double tt = TMIN + j * delt;
            double tau = tt / ( MJPsi2 ) ;
            tau_int = tau;
            q2sam = tt;

            double result = 0.0;

            for (int k=0; k<NY; k++) {

// Rapidity of the last emission in DGLAP
                double yrap =log(1./XMAX)+dely * k;

//	    Calculation of x
                double xxx=exp(-yrap);
                if(xxx>1.) xxx=1.;
                //double aux;

                double rapy;

// Rapidity gap from W, t and x
                rap_int = log(W2/sqrt(tt*(tt+MJPsi2))) - yrap; // the external variable

// Storing the original rapidity gsp size
                rapy=rap_int;

                rap_int = (alphas_fixed  * NC / pi )/ alphas_fixed_grid * rap_int ;
                // shift the rapidity to account
                //for the shift in alphas in the evolution

// Stuff from the original code:
                result = qgauss96(ampintegrand,LKMIN,k2tolk(q2sam/4.0)) + qgauss96(ampintegrand,k2tolk(q2sam/4.0),LKMAX);

                result *= alphas_fixed / alphas_fixed_grid; // rescaling the initial condition by ratio of alphas

                //running coupling : coupling of BFKL to impact factor
                result *= pow(pdf->alphasQ2(tt + MJPsi2) / alphas_fixed,2.0);

// Photon flux
                double PhiGam=PhotonFlux(W,sc);

                dsigmahatdt = 1.0/ (16.0 * pi) * pow(result * pi * (NC * NC - one) / (NC * NC),2.0) ;

// Cross section calculation and adding contribution
                aux += dely * dsigmahatdt * PhiGam * CONV_GeV2mb * CONV_mb2nb *
                       (4.0 * pow(NC,4.0)/pow(NC * NC-1.0, 2.0) * pdf->xfxQ2(21, xxx, tt + MJPsi2)
                        +pdf->xfxQ2(1, xxx, tt + MJPsi2)+pdf->xfxQ2(-1, xxx, tt + MJPsi2)
                        + pdf->xfxQ2(2, xxx, tt + MJPsi2)+pdf->xfxQ2(-2, xxx, tt + MJPsi2)
                        + pdf->xfxQ2(3, xxx, tt + MJPsi2)+pdf->xfxQ2(-3, xxx, tt + MJPsi2)
                        + pdf->xfxQ2(4, xxx, tt + MJPsi2)+pdf->xfxQ2(-4, xxx, tt + MJPsi2)
                        + pdf->xfxQ2(5, xxx, tt + MJPsi2)+pdf->xfxQ2(-5, xxx, tt + MJPsi2)
                       );

                qcont += dely * dsigmahatdt * PhiGam * CONV_GeV2mb * CONV_mb2nb *
                         (pdf->xfxQ2(1, xxx, tt + MJPsi2)+pdf->xfxQ2(-1, xxx, tt + MJPsi2)
                          + pdf->xfxQ2(2, xxx, tt + MJPsi2)+pdf->xfxQ2(-2, xxx, tt + MJPsi2)
                          + pdf->xfxQ2(3, xxx, tt + MJPsi2)+pdf->xfxQ2(-3, xxx, tt + MJPsi2)
                          + pdf->xfxQ2(4, xxx, tt + MJPsi2)+pdf->xfxQ2(-4, xxx, tt + MJPsi2)
                          + pdf->xfxQ2(5, xxx, tt + MJPsi2)+pdf->xfxQ2(-5, xxx, tt + MJPsi2));

                gcont += dely * dsigmahatdt * PhiGam * CONV_GeV2mb * CONV_mb2nb *
                         (4.0 * pow(NC,4.0)/pow(NC * NC-1.0, 2.0) * pdf->xfxQ2(21, xxx, tt + MJPsi2));


// Rapidity cut (minimal rapidity)
                if(rapy<dpy) aux=0., qcont=0., gcont=0.;


// Angular cut
//if(atan(sqrt(tt)/Ep/xxx)<4./180.*pi) aux=0;

// Energy cuts
                if(W>140.) aux=0., qcont=0., gcont=0.;
                if(W<50.) aux=0., qcont=0., gcont=0.;



// Print out
                //cout << fixed;
                //cout << setw(fieldwidth) << xxx << setw(fieldwidth) << tt << setw(fieldwidth) << W << setw(fieldwidth) << dsigmadt << endl;
                //fout << setw(fieldwidth) << xxx << setw(fieldwidth) << tt << setw(fieldwidth) << W << setw(fieldwidth) << dsigmadt << endl;
            }
            dsigmadt += aux * delt;
            quarkcont += qcont * delt;
            gluoncont += gcont * delt;
        }
        events = 2 * W * L * dsigmadt / sc / sc;
        cout << W << setw(fieldwidth) << events << setw(fieldwidth) << dsigmadt << endl;
        fout << W << setw(fieldwidth) << events << endl;
    }

    fout.close();
    // **************************************

    delete pdf;
    return 0;
}
// **************************************************************************************
// Photon Flux
double PhotonFlux(double W, double sc)
{
// Constants photon flux
    double Qmax=2.;
    double Qmax2=Qmax*Qmax;
    double me=0.001;
    double me2=me*me;
    double y=W*W/sc/sc;
    double y2=y*y;
    double alph=1./137.;

// Photon flux
    double PhiGam=(alph*(2*(Qmax2*(-1 + y) + me2*y2) +
                         Qmax2*(2 + (-2 + y)*y)*(log(Qmax2) -
                                                 log(-((me2*y2)/(-1 + y))))))/(2.*pi*Qmax2*y);

    return PhiGam;
}
// ***************************
// *************************************************************************************
double ampintegrand(double lk)
{
    double k2prim = ltok2(lk);
    double rm1 ;
    rm1 = k2prim + q2sam / four + s0const;

    //old interpolation
/*     double result =  k2prim / (rm1 * sqrt(rm1 * rm1 - k2prim * q2sam ) )  * phi0gammaV(k2prim,q2sam)
    * nfgluon.interpolate(tau_int,rap_int,k2prim);
 */

    // new interpolation
    double result =  k2prim / (rm1 * sqrt(rm1 * rm1 - k2prim * q2sam ) )  * phi0gammaV(k2prim,q2sam)
                     * alglib::spline3dcalc(sglob,tau_int,rap_int,k2prim);


    //Born calculation
    /*  double result=  k2prim / (rm1 * sqrt(rm1 * rm1 - k2prim * q2sam ) )  * phi0gammaV(k2prim,q2sam)
      * alphas_fixed_grid;
   */
    return result * two; // factor two comes from the Jacobian; changing dk^2/k^2 into dk/k
}
// **********************************
// *************************************************************************************
// gamma - vector meson impact factor
double phi0gammaV(double k2,double q2)
{
    double qpar2,qperp2;

    qpar2 = MassJPsi * MassJPsi / four;
    qperp2 = qpar2 + q2 / four;

    return Ccal * alphas_fixed / two * (one / qperp2 - one / (k2 + qpar2));
}
// *************************************************************************************
// *************************************************************************************
double k2tolk(double k2)
{
    return 0.5 * log(k2/K2MIN);
}
// *********************************************************************
double ltok2(double lk)
{

    return K2MIN * exp(lk * two);

}
// *************************************************************************************
double qgauss96(double(*func)(double),  double a,  double b)
{

    double xr,xm,dx,s;
    static   double x[]={0.016276744849603,0.04881298513605,0.081297495464426,
                         0.11369585011067,0.1459737146549,0.17809688236762,0.21003131046057,
                         0.24174315616384,0.27319881259105,0.3043649443545,0.33520852289263,
                         0.36569686147231,0.39579764982891,0.4254789884073,0.45470942216774,
                         0.4834579739206,0.51169417715467,0.53938810832436,0.5665104185614,
                         0.59303236477757,0.61892584012547,0.64416340378497,0.66871831004392,
                         0.69256453664217,0.71567681234897,0.7380306437444,0.75960234117665,
                         0.78036904386743,0.80030874413914,0.81940031073793,0.83762351122819,
                         0.8549590334346,0.8713885059093,0.88689451740242,0.90146063531585,
                         0.9150714231209,0.92771245672231,0.93937033975276,0.95003271778444,
                         0.95968829144874,0.96832682846326,0.97593917458514,0.98251726356301,
                         0.98805412632962,0.99254390032376,0.99598184298721,0.99836437586318,
                         0.99968950388323};
    static   double w[]={0.032550614492363,0.032516118713869,0.032447163714063,
                         0.032343822568572,0.03220620479402,0.032034456231969,0.031828758894364,
                         0.031589330770643,0.031316425596721,0.031010332586091,0.03067137612333,
                         0.030299915420329,0.029896344135616,0.029461089957174,0.028994614149195,
                         0.028497411063254,0.027970007616848,0.027412962726029,0.026826866725592,
                         0.026212340735672,0.025570036005349,0.024900633222484,0.024204841792365,
                         0.023483399085926,0.022737069658329,0.021966644438744,0.021172939892191,
                         0.020356797154333,0.019519081140145,0.018660679627411,0.017782502316045,
                         0.016885479864245,0.015970562902562,0.015038721026995,0.014090941772315,
                         0.013128229566962,0.012151604671088,0.011162102099838,0.010160770535008,
                         0.0091486712307834,0.0081268769256987,0.0070964707911537,0.0060585455042356,
                         0.0050142027429258,0.0039645543384356,0.0029107318178558,0.0018539607871533,
                         0.00079679206555206};

    xm=0.5*(a+b);
    xr=0.5*(b-a);
    s=0.0;
    for(int j=0;j<48;j++){
        dx=xr*x[j];
        s+=w[j]*((*func)(xm+dx)+(*func)(xm-dx));
    }
    return s*=xr;
}
// *************************************************************************************

