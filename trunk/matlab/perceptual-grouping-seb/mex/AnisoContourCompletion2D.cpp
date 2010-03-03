// -------------------------------------------------------------------
// -------------------------------------------------------------------
// 2D Anisotropic Tsitsiklis front propagation from a given point set
// (C) 2008 by Sebastien Bougleux, sebastien.bougleux@unicaen.fr
// -------------------------------------------------------------------
// -------------------------------------------------------------------
#include<mex.h>
#include "AnisoContourCompletion2D.h"
// -------------------------------------------------------------------
// -------------------------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] ) 
{
  //------------------------------------------------------------------
  // retrieve arguments
  if ((nrhs != 2) && (nrhs != 3) && (nrhs != 4))
      mexErrMsgTxt("[action,Voronoi,edges] = AnisoPropagation2D(tensor_field, point_set, stopping_ratio=0, [hx hy]=[1 1])");
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
  if (nlhs > 2)
  {
    plhs[2] = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
    V = (int*)mxGetPr(plhs[2]);
  }
  //==================================================================
  // Propagation
  AnisoContourCompletion2D<double, int> AP(TF, nx, ny, U, V, hx, hy);
  AP.greedyNearestFronts(seeds, nb_seeds);
  //==================================================================
  AP.resize();
  dims[0] = nx;
  dims[1] = ny;
  mxSetDimensions(plhs[0], dims, 2);
  if (nlhs > 1)
  {
    int dimsp[2] = { 2, AP.nbSaddlePoints() };
    plhs[1] = mxCreateNumericArray(2, dimsp, mxDOUBLE_CLASS, mxREAL);
    double *sp_mtrx = mxGetPr(plhs[1]);
    AP.getSaddlePoints(sp_mtrx);
  }
  if (nlhs > 2) mxSetDimensions(plhs[2], dims, 2);
  if (nlhs > 3)
  {
      int dime[2] = { 2, AP.nbEdges() };
      plhs[3] = mxCreateNumericArray(2, dime, mxINT32_CLASS, mxREAL);
      int *edges_mtrx = (int*)mxGetPr(plhs[3]);
      AP.getEdges(edges_mtrx);
  }
  //==================================================================
  return;
}
