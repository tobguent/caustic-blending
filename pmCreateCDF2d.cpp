#include "mex.h"
#include <stdlib.h>
#include <cmath>
#include <vector>

void doWork(double *image, int width, int height, double *colCDF, double *rowCDF, double* prob)
{
    // accumulate function values for normalization
	double totalValue = 0;
	for (int i=0; i<width*height; ++i)
		totalValue += image[i];
	
	// create probability distribution
	for (int i=0; i<width*height; ++i)
		prob[i] = image[i] / totalValue;
    
	// compute CDF
	for (int y=0; y<height; ++y)
	{
		// prefix sum of row
		double rowSum = 0;
		for (int x=0; x<width; ++x)
		{
			int linIndex = y*width+x;

			rowSum += prob[linIndex];
			colCDF[linIndex] = rowSum;
		}

		// store the sum in a column vector
		rowCDF[y] = rowSum;

		// normalize prefix sum of row to sum up to 1
        if (rowSum != 0)
        {
            for (int x=0; x<width; ++x)
            {
                int linIndex = y*width+x;
                colCDF[linIndex] /= rowSum;
            }
        }
	}

	// prefix sum of the column to actually turn it into the CDF
	{
		double sum=0;
		for (int y=0; y<height; ++y)
		{
			sum += rowCDF[y];
			rowCDF[y] = sum;
		}
	}
}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    double *image;
    
	int width, height;
	double *colCDF, *rowCDF, *prob;

    if (nrhs != 1) {
        mexErrMsgIdAndTxt("pmInverseCDF2d:nrhs","One input required.");
        return;
    }
    
    if (nlhs != 3) {
        mexErrMsgIdAndTxt("pmInverseCDF2d:nlhs","Three outputs required.");
        return;
    }
    
	/* Input arguments */
	width = (int)mxGetM(prhs[0]);
	height = (int)mxGetN(prhs[0]);
	image = (double*)mxGetPr(prhs[0]);
	
    /* Output arguments */
    plhs[0] = mxCreateNumericMatrix(width, height, mxDOUBLE_CLASS, mxREAL); // width x height
    colCDF = (double*) mxGetData(plhs[0]);
    
    plhs[1] = mxCreateNumericMatrix(1, height, mxDOUBLE_CLASS, mxREAL);  // 1 x height
    rowCDF = (double*) mxGetData(plhs[1]);
    
    plhs[2] = mxCreateNumericMatrix(width, height, mxDOUBLE_CLASS, mxREAL); // width x height
    prob = (double*) mxGetData(plhs[2]);
    
    
	/* Call C-function */
	doWork(image, width, height, colCDF, rowCDF, prob);
}