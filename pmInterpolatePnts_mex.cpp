#include "mex.h"
#include <stdlib.h>
#include <cmath>

// struct for type double2
struct double2{
  double2() : x(0), y(0) {}
  double2(double xx, double yy) : x(xx),y(yy) {}
  double x;
  double y;
};

inline static double2 operator+(const double2& A, const double2& B) {
    return double2(A.x+B.x, A.y+B.y);
}
inline static double2 operator-(const double2& A, const double2& B) {
    return double2(A.x-B.x, A.y-B.y);
}
inline static double2 operator-(const double2& A) {
    return double2(-A.x, -A.y);
}
inline static double2 operator*(const double2& A, const double& s) {
    return double2(A.x*s, A.y*s);
}
inline static double2 operator*(const double& s, const double2& A) {
    return double2(A.x*s, A.y*s);
}
inline static double2 operator/(const double2& A, const double& s) {
    return double2(A.x/s, A.y/s);
}

void doWork(double2* k0, double2* k1, double2* k2, double t, double* rr, int N, double2* x)
{
    for (int i=0; i<N; ++i)
    {
        double r = rr[i];
        double c = t/r;
        double2 b0,b1,b2,b3;
        if (t <= r)
        {
            c = t/r;
            b0 = k0[i];
            b1=-(k0[i]*(r*r-5*r+4)-k2[i]*r*r-k1[i]*r+2*k1[i])/(6*r-6);
            b2=-(k0[i]*(r*r-2*r+1)-k2[i]*r*r-k1[i]*r+2*k1[i])/(3*r-3);
            b3 = k1[i];
        }
        else
        {
            c = (t-r)/(1-r);
            b0 = k1[i];
            b1=-(k0[i]*(r*r-2*r+1)-k2[i]*r*r-k1[i]*r-k1[i])/(3*r);
            b2=-(k0[i]*(r*r-2*r+1)+k2[i]*(-r*r-3*r)-k1[i]*r-k1[i])/(6*r);
            b3 = k2[i];
        }

        x[i] =     (1-c)*(1-c)*(1-c)*b0
             + 3.0*(1-c)*(1-c)*   c *b1
             + 3.0*(1-c)*   c *   c *b2
             +        c *   c *   c *b3;
    }
}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    double2 *k0, *k1, *k2, *x;
    double *t, *r;
    
	int nOfRows, nOfColumns;
	
    if (nrhs != 5) {
        mexErrMsgIdAndTxt("pmInterpolatePnts_mex:nrhs","Five inputs required.");
        return;
    }
    
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("pmInterpolatePnts_mex:nlhs","One output required.");
        return;
    }
    
	/* Input arguments */
	nOfRows    = (int)mxGetM(prhs[0]);  // 2
	nOfColumns = (int)mxGetN(prhs[0]);  // N
	k0 = (double2*)mxGetPr(prhs[0]);
	k1 = (double2*)mxGetPr(prhs[1]);
    k2 = (double2*)mxGetPr(prhs[2]);
    
    if (nOfRows != 2) {
        mexErrMsgIdAndTxt("pmInterpolatePnts_mex:prhs","Points must be 2D.");
        return;
    }
    
    t = (double*)mxGetPr(prhs[3]);
	r = (double*)mxGetPr(prhs[4]);

    /* Output arguments */
    plhs[0] = mxCreateNumericMatrix(nOfRows, nOfColumns, mxDOUBLE_CLASS, mxREAL);  // 2xN
    x = (double2*) mxGetData(plhs[0]);
    
	/* Call C-function */
	doWork(k0, k1, k2, *t, r, nOfColumns, x);
}