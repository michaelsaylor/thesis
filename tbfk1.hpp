//
//  tbfkl.hpp
//
//
//  Created by Anna on 8/24/17.
//
//

#ifndef tbfkl_hpp
#define tbfkl_hpp
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <cmath>
#include <time.h>
using namespace std;

typedef unsigned int counter; // define the type for the counter in for loops
const double s0const = 0.5;  // regulator cutoff in the denominators; present in BFKL solution and in the cross section
const double K2MAX = 1.0e+10;  // max value for K2
const double K2MIN = 1.0e-4;  // min value for K2
const double LOGK2MAX = log(K2MAX); // log cutoffs
const double LOGK2MIN = log(K2MIN);
const double LKMIN = 0.0;
const double LKMAX = 0.5 * log(K2MAX/K2MIN);
const double YMAX = 20.0;     // max range in rapidity
const double YMIN =  0.0;     // min range in rapidity
const int    NCHEB = 100;      // number of polynomials
const double DNCHEB = static_cast<double>(NCHEB);   // the same but double precision
const double two_over_cheb = 2.0/DNCHEB;            // conveninent factor
const double four_over_2cheb = 4.0/(DNCHEB*DNCHEB); // -------||--------
const int    NQT = 1; // number of points in q, momentum transfer
const int    NT    = 400;     // number of steps in Y
const double one   = 1.0;
const double two   = 2.0;
const double half  = 0.5;
const double four  = 4.0;
const double pi = atan(1.0) * 4.0;                  //guess what it is?
const double two_pi = 2.0 * pi;                     //
const int    NC = 3; // number of colors
const double alphas_fixed = 0.17; // fixed strong coupling, it is alphasbar
const double alphas_fixed_grid = 0.1; // grid alphasbar
const double MassJPsi = 3.096 ;  // mass of J/Psi
const double MV2 = MassJPsi * MassJPsi;
const double keV_to_GeV = 1.0e-6;
const double Gammaee = 5.55 * keV_to_GeV;
const double Gammall = 11.08 * keV_to_GeV;
const double alphaem = 1.0/137.0356; // alpha electromagnetic
const double Ccal  = pow(3.0 * Gammaee * pow(MassJPsi,3.0) / alphaem,0.5);
const double CONV_GeV2mb = 0.389 ; // GeV2mb
const double CONV_mb2nb = 1.0e+6; //mb to nb
const double CONV_mb2pb = 1.0e+9; //mb to pb
const double FACT = 1.0 / pow(pi,3.0);
const double Qc = 2.0 / 3.0;
const double gJPsi = pow(3.0 * MassJPsi * Gammaee/(16.0 * pi),0.5) / (alphaem * Qc);
const double gJPsi_ll = pow(3.0 * MassJPsi * Gammall/(16.0 * pi),0.5) / (alphaem * Qc);
//const double B_toleptons = 0.05971 + 0.05961;
const double B_toleptons = 0.05961;
const double Cgq = (NC * NC - 4.0) * (NC * NC - 2.0) / (16.0 * pow(NC,7.0)) ;
const double Cgg = 3.0 / 8.0 *  (NC * NC - 4.0) / ((NC * NC - 1.0) * (NC * NC * NC));

#endif /* tbfkl_hpp */