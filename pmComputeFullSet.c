#include "mex.h"
#include <stdlib.h>
#include <math.h>

# if defined(__cplusplus)
#  define INLINE inline
# else
#  define INLINE
# endif

/* struct for type double2 */
typedef struct {
  double x;
  double y;
} double2;

/* computes distance between two points */
INLINE static double double2_dist(double2* A, double2* B) {
    return sqrt((A->x-B->x)*(A->x-B->x) + (A->y-B->y)*(A->y-B->y));
}

/* computes a distance matrix */
INLINE static void compute_dist(double2* A, double2* B, int N, double* dist)
{
    int x,y;
    for (y=0; y<N; ++y)
        for (x=0; x<N; ++x)
            dist[y*N+x] = double2_dist(&A[x], &B[y]);
}

/* computes the frobenius norm between the difference of two matrices (the latter is permutated) */
INLINE static double frobenius_norm(double* distA, double* distB, int N, int* p)
{
    int x,y;
    double temp, result = 0;
    for (y=0; y<N; ++y)
        for (x=0; x<N; ++x)
        {
            temp = (distA[y*N+x] - distB[p[y]*N+p[x]]);
            result += temp*temp;
        }
    return sqrt(result);
}

/* computes the trace of a matrix (colums are permutated) */
INLINE static double trace(double* distAB, int N, int* p)
{
    int x;
    double result = 0;
    for (x=0; x<N; ++x)
        result += distAB[p[x]*N+x];
    return result;
}

/* evaluate error for permutation p */
INLINE static void evalerr(double* distA, double* distB, double* distAB, int N, int* p, double beta, double* e, double* ei, double* ex)
{
    *ei = frobenius_norm(distA, distB, N, p);   /* ei=norm(A-B(p,p),'fro'); */
    *ex = trace(distAB, N, p);                  /* ex=trace(AB(:,p)); */
	*e=(1-beta)*(*ei)/N+beta*(*ex);
}

/* permutes two elements in a vector */
INLINE static void permute(int* p, int j, int k) {
    int temp = p[j];
    p[j] = p[k];
    p[k] = temp;
}

/* computes squared frobenius norm of the difference matrix between A and B, only on columns and rows (j,k) */
INLINE static double frob_sqnorm_sel(double* distA, double* distB, int N, int* p, int j, int k)
{
    int x,y;
    double temp, result = 0;

    /* column cost */
    for (y=0; y<N; ++y)
    {
        temp = (distA[y*N+j] - distB[p[y]*N+p[j]]);
        result += temp*temp;
        temp = (distA[y*N+k] - distB[p[y]*N+p[k]]);
        result += temp*temp;
    }

    /* row cost */
    for (x=0; x<N; ++x)
    {
        temp = (distA[j*N+x] - distB[p[j]*N+p[x]]);
        result += temp*temp;
        temp = (distA[k*N+x] - distB[p[k]*N+p[x]]);
        result += temp*temp;
    }

    return result;
}

INLINE static void doWork(double2* A, double2* B, int N, int* p, double beta)
{
    int i, j, k;
    double *distA, *distB, *distAB;
    double e, ei, ex, f, g, emin, emini, eminx;
    int progress=0;

    /* allocate distance matrices */
    distA=(double*) malloc(N*N*sizeof(double));
    distB=(double*) malloc(N*N*sizeof(double));
    distAB=(double*) malloc(N*N*sizeof(double));

    /* compute the distance matrices */
    compute_dist(A, A, N, distA);
    compute_dist(B, B, N, distB);
    compute_dist(A, B, N, distAB);

    /* initialize permutation p */
    for (i=0; i<N; ++i)
        p[i] = i;

    /* compute the initial error */
    evalerr(distA, distB, distAB, N, p, beta, &emin, &emini, &eminx);

    for (i=0; i<N; ++i)  /* % N iterations is just a "guess" */
    {
        progress=0;
        for (j=0; j<N; ++j)
        {
            for (k=0; k<N; ++k)
            {
                if (j != k)
                {
                    /* compute cost of current state (j,k) */
                    f = emini * emini;
                    g = eminx;
                    f -= frob_sqnorm_sel(distA, distB, N, p, j, k);
                    g -= distAB[p[j]*N+j] + distAB[p[k]*N+k];

                    /* permute */
                    permute(p, j, k);

                    /* compute cost of new state (k,j) */
                    f += frob_sqnorm_sel(distA, distB, N, p, j, k);
                    g += distAB[p[j]*N+j] + distAB[p[k]*N+k];

                    f=sqrt(f);
                    e=(1-beta)*f/N+beta*g;

                    /* accept if error decreased */
                    if (e<emin)
                    {
                        emin=e;
                        emini=f;
                        eminx=g;
                        progress=1;
                    }
                    else
                       permute(p, j, k);    /* undo the permutation */
                }
            }
        }

        if (progress==0) /* more? */
            break;
    }

    /* free the distance matrices */
    free(distA);
    free(distB);
    free(distAB);

    /* add the index shift for matlab... */
    for (i=0; i<N; ++i)
        p[i]++;
}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double2 *setA, *setB;
	int nOfRowsA, nOfColumnsA, nOfRowsB, nOfColumnsB;
	int* p;
    double beta;

    if (nrhs != 3) {
        mexErrMsgIdAndTxt("pmComputeFullSet:nrhs","Three inputs required.");
        return;
    }

    if (nlhs != 1) {
        mexErrMsgIdAndTxt("pmComputeFullSet:nlhs","One output required.");
        return;
    }

	/* Input arguments */
	nOfRowsA    = (int)mxGetM(prhs[0]);  /* 2 */
	nOfColumnsA = (int)mxGetN(prhs[0]);  /* N */
	setA = (double2*)mxGetPr(prhs[0]);

    nOfRowsB    = (int)mxGetM(prhs[1]);  /* 2 */
	nOfColumnsB = (int)mxGetN(prhs[1]);  /* N */
	setB = (double2*)mxGetPr(prhs[1]);

    beta = mxGetScalar(prhs[2]);

    if (nOfRowsA != nOfRowsB || nOfColumnsA != nOfColumnsB) {
        mexErrMsgIdAndTxt("pmComputeFullSet:prhs","Input matrices of different size.");
        return;
    }

    if (nOfRowsA != 2) {
        mexErrMsgIdAndTxt("pmComputeFullSet:prhs","Points must be 2D.");
        return;
    }

    /* Output arguments */
    plhs[0] = mxCreateNumericMatrix(1, nOfColumnsA, mxINT32_CLASS, mxREAL);  /* 1xN */
    p = (int*) mxGetData(plhs[0]);

	/* Call C-function */
	doWork(setA, setB, nOfColumnsA, p, beta);
}
