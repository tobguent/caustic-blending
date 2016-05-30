#include "mex.h"

#include <stdlib.h>

void doWork(double* density_in, int width, int height, int* rx, int* ry, double* flux, double* density_out, int num_samples)
{
    // copy old data
    for (int i=0; i<width*height; ++i)
        density_out[i] = density_in[i];
            
    // add new data
    for (int i=0; i<num_samples; ++i)
    {
        int x = rx[i];
        int y = ry[i];
        if (x >= 0 && x < width && y>=0 && y<height)
            density_out[y*width+x] += flux[i];
    }
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	int width, height;
	double *density_in, *density_out;
    
    int num_samples;   
    double *flux;

    int *rx, *ry;
    
    if (nrhs != 4) {
        mexErrMsgIdAndTxt("pmScatter:nrhs","Four inputs required.");
        return;
    }
    
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("pmScatter:nlhs","One output required.");
        return;
    }
    
	/* Input arguments */
	width = (int)mxGetM(prhs[0]);
	height = (int)mxGetN(prhs[0]);
    
	density_in = (double*)mxGetPr(prhs[0]);
    rx = (int*)mxGetPr(prhs[1]);
    ry = (int*)mxGetPr(prhs[2]);
    flux = (double*)mxGetPr(prhs[3]);
    
    num_samples = (int)mxGetN(prhs[3]);
    
    /* Output arguments */
    plhs[0] = mxCreateNumericMatrix(width, height, mxDOUBLE_CLASS, mxREAL);
    density_out = (double*) mxGetData(plhs[0]);
    
	/* Call C-function */
	doWork(density_in, width, height, rx, ry, flux, density_out, num_samples);
}