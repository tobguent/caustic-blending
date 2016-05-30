mex pmComputeFullSet.c
mex pmCreateCDF2d.cpp
mex CXXFLAGS="\$CXXFLAGS -std=c++11"  pmSampleInvCDF2d.cpp
mex CXXFLAGS="\$CXXFLAGS -std=c++11"  pmRefinement.cpp -I"./";
mex pmInterpolatePnts_mex.cpp;
mex pmScatter.cpp
