#include "mex.h"
#define KDTREE_USE_SQUARED_DISTANCE
#include <kdtree++/kdtree.hpp>
#include <stdlib.h>
#include <list>
#include <set>

struct Pnt {
	int id;
	double coord[2];
    bool operator==(const Pnt& other) const {
		return id == other.id && coord[0]==other.coord[0] && coord[1]==other.coord[1];
	}
};
inline double tac2( Pnt t, size_t k) { return t.coord[k]; }
typedef KDTree::KDTree<2, Pnt, std::pointer_to_binary_function<Pnt,size_t,double> > tree_type2;

// struct for type double2
typedef struct {
  double x;
  double y;
} double2;

struct Match {
  double d;
  int A;
  int B;
};

struct Comparator {
  bool operator()(const Match& lhs, const Match& rhs) {
    return lhs.d > rhs.d;
  }
};

/* computes squared distance between two points */
inline static double double2_dist2(double2* A, double2* B) {
	return (A->x-B->x)*(A->x-B->x) + (A->y-B->y)*(A->y-B->y);
}

void doWork(double2* Adata, double2* Bdata, int N, int* cdata)
{
    // construct kd-tree
	tree_type2 Btree(std::ptr_fun(tac2));
	for (int i=0; i<N; ++i)
	{
		Pnt p = {i, Bdata[i].x, Bdata[i].y};
		Btree.insert(p);
	}
	Btree.optimise(); // balance the tree

	// generate a multiset for points in A (sorts them by match distance and stores the matches)
	std::multiset<Match, Comparator> Aset;

	// list of incoming matches at B
	std::vector<std::list<int>> Blist(N);

	// first step: points in A choose closest point in B (fill the multiset)
	for (int i=0; i<N; ++i)
	{
		Pnt query = {-1, Adata[i].x, Adata[i].y};
		std::pair<tree_type2::const_iterator,double> found = Btree.find_nearest(query);
		if (found.first != Btree.end())
		{
			Match m = {found.second,i,found.first->id};
			Aset.insert(m);
			Blist[m.B].push_back(m.A);
		}
		else
		{
            mexErrMsgIdAndTxt("pmRefinement:initMatch","Couldn't find initial match.");
            return;
		}
	}

	// loop until all have a match
	while (!Aset.empty())
	{
		// get the match of the assignment with largest distance (among the candidates that haven't been fixed yet)
		Match m = *Aset.begin();

		// store the match
		cdata[m.A] = m.B;

		// remove the match from the candidate list
		Aset.erase(m);

		// remove B point from the kd tree
		Pnt pntB = {m.B, Bdata[m.B].x, Bdata[m.B].y};
		Btree.erase_exact(pntB);

		// iterate the other matches that tried to connect with B too and update them
		if (Blist[m.B].size() > 1)
		{
			for (auto it=Blist[m.B].begin(); it!=Blist[m.B].end(); ++it)
			{
				int a = *it;
				if (a != m.A) // skip the match that we found in this iteration
				{
					// do a new query (now, with B removed)
					Pnt query = {-1, Adata[a].x, Adata[a].y};
					std::pair<tree_type2::const_iterator,double> found = Btree.find_nearest(query);
					int cc = found.first->id;  // new closest match
					double dd = found.second;  // new squared distance

                    // remove the old link from the set
					double dist = (double2_dist2(&Adata[a], &Bdata[m.B]));
					Match md = {dist, a, m.B};
					Aset.erase(md);

                    // add in the new link into the set
					md.B = cc;
					md.d = dd;
					Aset.insert(md);
					Blist[md.B].push_back(md.A); // record this incoming link at the new B node
				}
			}

		}
		Blist[m.B].clear(); // the processed node does no longer need to store its incoming links.
	}
}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double2 *setA, *setB;
	int nOfRowsA, nOfColumnsA, nOfRowsB, nOfColumnsB;
	int* p;

    if (nrhs != 2) {
        mexErrMsgIdAndTxt("pmRefinement:nrhs","Two inputs required.");
        return;
    }

    if (nlhs != 1) {
        mexErrMsgIdAndTxt("pmRefinement:nlhs","One output required.");
        return;
    }

	/* Input arguments */
	nOfRowsA    = (int)mxGetM(prhs[0]);  // 2
	nOfColumnsA = (int)mxGetN(prhs[0]);  // N
	setA = (double2*)mxGetPr(prhs[0]);

    nOfRowsB    = (int)mxGetM(prhs[1]);  // 2
	nOfColumnsB = (int)mxGetN(prhs[1]);  // N
	setB = (double2*)mxGetPr(prhs[1]);

    if (nOfRowsA != nOfRowsB || nOfColumnsA != nOfColumnsB) {
        mexErrMsgIdAndTxt("pmRefinement:prhs","Input matrices of different size.");
        return;
    }

    if (nOfRowsA != 2) {
        mexErrMsgIdAndTxt("pmRefinement:prhs","Points must be 2D.");
        return;
    }

    /* Output arguments */
    plhs[0] = mxCreateNumericMatrix(1, nOfColumnsA, mxINT32_CLASS, mxREAL);  // 1xN
    p = (int*) mxGetData(plhs[0]);

	/* Call C-function */
	doWork(setA, setB, nOfColumnsA, p);

    /* add the index shift for matlab... */
    for (int i=0; i<nOfColumnsA; ++i)
        p[i] += 1;
}
