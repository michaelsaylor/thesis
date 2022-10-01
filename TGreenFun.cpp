//
//  TGreenFun.cpp
//  interpolator
//
//  Created by Anna Stasto on 3/30/18.
//
#include <math.h>
#include "/home/mike/thesis/TGreenFun.hpp"
#include <iostream>

using namespace std;

//*************************************************************************************
inline double TGreenFun::ltok(double lk2)
{
    return k2min * exp(lk2);
}
//*************************************************************************************
inline double TGreenFun::k2tol(double k2)
{
    return log(k2 / k2min);
}
//*************************************************************************************
double TGreenFun::del(int nd,double am,double ax)
{
    if (nd > 0) {
        return (ax - am) / static_cast<double>(nd);
    }
    else exit(1);
}
//*************************************************************************************
// routine for trilinear interpolation
// arguments are : tau , rapidity and k2
// interpolation is in rapidity, logtau and log (k2)
double TGreenFun::interpolate(double tau,double y,double k2){


    double ltau = log(tau / tmin);
    double lk2 = k2tol(k2);

    double yd = (y - ymin) / (ymax - ymin);
    double lk2d = (lk2 - lk2min) / (lk2max - lk2min);
    double ltaud = (ltau  - ltaumin) / (ltaumax - ltaumin) ;

    // find one corner of the cube
    int i1 = static_cast<int>(ltau / delt);
    int i2 = static_cast<int>(y / dely);
    int i3 = static_cast<int>(lk2 / delk);

    double c00 = GF[i1][i2][i3]     * (1.0 - ltaud) + GF[i1 + 1][i2][i3]          * ltaud;
    double c01 = GF[i1][i2][i3+1]   * (1.0 - ltaud) + GF[i1 + 1][i2][i3 + 1]      * ltaud;
    double c10 = GF[i1][i2+1][i3]   * (1.0 - ltaud) + GF[i1 + 1][i2 + 1][i3]      * ltaud;
    double c11 = GF[i1][i2+1][i3+1] * (1.0 - ltaud) + GF[i1 + 1][i2 + 1 ][i3 + 1] * ltaud;

    double c0 = c00 * (1.0 - yd) + c10 * yd;
    double c1 = c01 * (1.0 - yd) + c11 * yd;

    return c0 * (1.0 - lk2d) + c1 * lk2d;
}
//*************************************************************************************
// book the matrix
double***  TGreenFun::book_matrix( int size1, int size2, int size3){

    double*** arr = new double**[size1];

    for (int i = 0; i < size1; i++) {
        arr[i] = new double*[size2];
        for (int j = 0; j < size2; j++) {
            arr[i][j] = new double[size3]();
        }
    }

    return arr;
}
//*************************************************************************************
TGreenFun::TGreenFun(int n1,int n2,int n3,double tm,double tx,double y1,double y2,double kmn2,double kmx2){


    NT =  n1;       //  dim. in tau
    NY =  n2;       // dim. in rapidity
    NK =  n3;       // dim. in kt

    ymin = y1;
    ymax = y2;

    k2min = kmn2;
    lk2min = k2tol(k2min);

    k2max  = kmx2;
    lk2max = k2tol(k2max);

    tmin = tm;
    ltaumin = 0.0;

    tmax = tx;
    ltaumax = log(tmax/tmin);

    dely = del(NY, ymin, ymax);
    delk = del(NK,lk2min,lk2max);
    delt = del(NT, ltaumin, ltaumax);

    GF = book_matrix(NT+1, NY+1, NK+1); // note that the sizes of the grid are always Nx+1

}