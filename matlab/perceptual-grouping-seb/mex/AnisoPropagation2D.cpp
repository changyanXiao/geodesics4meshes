// -------------------------------------------------------------------
// -------------------------------------------------------------------
// 2D Anisotropic Tsitsiklis front propagation from a given point set
// (C) 2008 by Sebastien Bougleux, sebastien.bougleux@unicaen.fr
// -------------------------------------------------------------------
// -------------------------------------------------------------------
#include<mex.h>
#include "AnisoPropagation2D.h"
// -------------------------------------------------------------------
// -------------------------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] ) 
{
  //------------------------------------------------------------------
  // retrieve arguments
  if ((nrhs != 2) && (nrhs != 3) && (nrhs != 4))
      mexErrMsgTxt("[U,V] = AnisoPropagation2D(tensor_field, point_set, stopping_ratio=0, [hx hy]=[1 1])");
  //------------------------------------------------------------------
  // First input argument: tensor field
  if ((mxGetNumberOfDimensions(prhs[0]) != 4) ||
      (mxGetDimensions(prhs[0])[2] != 2) ||
      (mxGetDimensions(prhs[0])[3] != 2))
      mexErrMsgTxt("Library error: tensor field must be a 2Dx2x2 symmetric definite matrix.");	
  int nx = mxGetDimensions(prhs[0])[0], ny = mxGetDimensions(prhs[0])[1];
  double *TF = mxGetPr(prhs[0]);
  //------------------------------------------------------------------
  // Second input argument: seeds array
  double *seeds = mxGetPr(prhs[1]);
  int nb_seeds = mxGetN(prhs[1]);
  if (nb_seeds == 0) mexErrMsgTxt("Library error: point set must be non-empty.");	
  //------------------------------------------------------------------
  // Third input argument: stopping criterion
  double fmax = 0;
  if (nrhs > 2) fmax = (double)mxGetScalar(prhs[2]);
  //------------------------------------------------------------------
  // Fourth input argument: spacing and dimensions
  double hx = 1, hy = 1;
  if (nrhs > 3)
  {
     const int* dim_h = mxGetDimensions(prhs[3]);
     if ((dim_h[0] != 2) || (dim_h[1] != 1)) mexErrMsgTxt("Library error: h must be a 2x1 array list.");
     hx = mxGetPr(prhs[3])[0];
     hy = mxGetPr(prhs[3])[1];
  }
  //------------------------------------------------------------------
  // First output : minimal action map
  int dims[2] = {nx+2, ny+2};
  plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
  double *U = mxGetPr(plhs[0]);
  //------------------------------------------------------------------
  // Second output : Voronoi diagramm
  int *V = NULL;
  if (nlhs > 1)
  {
    plhs[1] = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
    V = (int*)mxGetPr(plhs[1]);
  }
  //==================================================================
  // Propagation
  AnisoPropagation2D<double, int> AP(TF, nx, ny, U, V, hx, hy);
  AP.propagation(seeds, nb_seeds, fmax);
  AP.resize();
  dims[0] = nx;
  dims[1] = ny;
  mxSetDimensions(plhs[0], dims, 2);
  if (nlhs > 1) mxSetDimensions(plhs[1], dims, 2);
  //==================================================================
  return;
}
