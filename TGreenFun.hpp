//
//  TGreenFun.hpp
//  interpolator
//
//  Created by Anna Stasto on 3/30/18.
//

#ifndef TGreenFun_hpp
#define TGreenFun_hpp

class TGreenFun{

public:

    double *** GF; // this will be a matrix to hold the solutions
    int NT;       //  dim. in tau
    int NY;       // dim. in rapidity
    int NK;       // dim. in kt

    double k2min; // min. value in k2
    double lk2min; // min. value in logarithm k2

    double k2max; // max. value in k2
    double lk2max; // max value in logarithm k2

    double ymin; // min. value in rapidity
    double ymax; //max value in rapidity

    double tmin;  // min. value in  logtau
    double ltaumin; // max value in logtau

    double tmax; // max value in tau
    double ltaumax; // max value in log tau

    double dely,delt,delk; //spacings in grid in rapidity,log(tau),log(kT2)

    // function to return the interpolated value
    double interpolate(double tau,double y,double k2);

    // book matrix
    double *** book_matrix( int s1, int s2, int s3);

    // conversion from logarithmic variable to k2

    inline double ltok(double lk2);

    // conversion from k2 to logarithmic variable
    inline double k2tol(double k2);

    // returns the delta given number of divisions and the range
    double del(int Ndiv,double argmin,double argmax);

    TGreenFun(int n1,int n2,int n3,double tm,double tx,double y1,double y2,double kmn2,double kmx2);
};



#endif /* TGreenFun_hpp */