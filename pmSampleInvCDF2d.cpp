#include "mex.h"

#include "stdlib.h"
#include <cmath>
#include <vector>
#include <random>
#include <time.h>


int binarySearch(int leftIndex, int rightIndex, double value, double* buffer, int offset)
{
	while (rightIndex>leftIndex)
	{
		int centerIndex = (leftIndex+rightIndex)/2;
		double centerValue = buffer[offset+centerIndex]; 

		if (value < centerValue)
			rightIndex = centerIndex;
		else
			leftIndex = centerIndex+1;
	}
	return leftIndex;
}

void doWork(double* colCDF, double* rowCDF, double* prob, double* image, int width, int height, double* samples_in, double* samples_out, int num_samples, int num_channels, double* flux)
{
    const double OUTLIER_COORD = -1;
    const double AreaElement = 2.0 / width * 2.0 / height;
    
    for (int i=0; i<num_samples; ++i)
    {
        // randomly pick a row
        double pRow = samples_in[4*i+0];
        if (pRow < 0 || pRow > 1) {     // if outside, reject (flux=0)
            samples_out[i*2+0] = OUTLIER_COORD;
            samples_out[i*2+1] = OUTLIER_COORD;
            for (int c=0; c<num_channels; ++c)
                flux[i*num_channels+c] = 0;
            continue;
        }
        double tRow = samples_in[4*i+2] - 0.5;  //  [-0.5, 0.5]
        int iRow = binarySearch(0, height-1, pRow, rowCDF, 0);

        // randomly pick a column
        double pCol = samples_in[4*i+1];
        if (pCol < 0 || pCol > 1) {     // if outside, reject (flux=0)
            samples_out[i*2+0] = OUTLIER_COORD;
            samples_out[i*2+1] = OUTLIER_COORD;
            for (int c=0; c<num_channels; ++c)
                flux[i*num_channels+c] = 0;
            continue;
        }
        double tCol = samples_in[4*i+3] - 0.5;  //  [-0.5, 0.5]
        int iCol = binarySearch(0, width-1, pCol, colCDF, iRow*width);

        // sample randomly in the pixel
        samples_out[i*2+1] = (1-(iCol + tCol) / width) * 2 - 1;
        samples_out[i*2+0] = ((iRow + tRow) / height) * 2 - 1;
        
        // read probability of this pixel
        double probability = prob[iRow*width+iCol];
        
        for (int c=0; c<num_channels; ++c)
        {
            double value = image[iCol + iRow * width + c * (width * height)];
            // generate the flux
            flux[i*num_channels+c] = value / probability / num_samples * AreaElement;
        }
    }
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	int width, height;
	double *colCDF, *rowCDF, *prob, *image;
    
    int num_samples, num_channels;
    double *samples_in, *samples_out, *flux;

    if (nrhs != 5) {
        mexErrMsgIdAndTxt("pmSampleInvCDF2d:nrhs","Five inputs required.");
        return;
    }
    
    if (nlhs != 2) {
        mexErrMsgIdAndTxt("pmSampleInvCDF2d:nlhs","Two outputs required.");
        return;
    }
    
	/* Input arguments */
	width = (int)mxGetM(prhs[0]);
	height = (int)mxGetN(prhs[0]);
    
	colCDF = (double*)mxGetPr(prhs[0]);
    rowCDF = (double*)mxGetPr(prhs[1]);
    prob = (double*)mxGetPr(prhs[2]);
    image = (double*)mxGetPr(prhs[3]);
    samples_in = (double*)mxGetPr(prhs[4]);
    
    num_samples = (int)mxGetN(prhs[4]);
    num_channels = 1;
    if (mxGetNumberOfDimensions(prhs[3]) == 3)
        num_channels = (int)mxGetDimensions(prhs[3])[2];
    
    /* Output arguments */
    plhs[0] = mxCreateNumericMatrix(2, num_samples, mxDOUBLE_CLASS, mxREAL);
    samples_out = (double*) mxGetData(plhs[0]);
    plhs[1] = mxCreateNumericMatrix(num_channels, num_samples, mxDOUBLE_CLASS, mxREAL);
    flux = (double*) mxGetData(plhs[1]);
    
	/* Call C-function */
	doWork(colCDF, rowCDF, prob, image, width, height, samples_in, samples_out, num_samples, num_channels, flux);
}