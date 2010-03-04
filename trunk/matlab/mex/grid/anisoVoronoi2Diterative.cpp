#include "anisoVoronoi2Diterative.h"

void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray*prhs[] ) 
{
  
  //------------------------------------------------------------------
  /* retrive arguments */
  //------------------------------------------------------------------
  if (nrhs != 2 && nrhs != 4 && nrhs != 5 && nrhs != 6)
    mexErrMsgTxt("USAGE: [U,V,dUx,dUy] = anisoVoronoi2Diterative(TensorField, SeedMatrix, [DistanceMatrix, LabelMatrix, SeedsSubsetVector, tolerance])");
  if (mxGetNumberOfDimensions(prhs[0]) != 4)
    mexErrMsgTxt("HERE T must be a 2D x 2x2 tensor field of symmetric definite matrices.");
  //------------------------------------------------------------------
  //hx = 1; hy = 1;
  //hx2 = hx*hx; hy2 = hy*hy;
  //hxhy = hx*hy;
  nx = mxGetDimensions(prhs[0])[0];
  ny = mxGetDimensions(prhs[0])[1];
  if( (mxGetDimensions(prhs[0])[2] != 2) || (mxGetDimensions(prhs[0])[3] != 2) )
    mexErrMsgTxt("T must be a 2D x 2x2 tensor field of symmetric definite matrices.");
  Nx = nx+2; Ny = ny+2;
  size = Nx*Ny; nxny = nx*ny;
  //------------------------------------------------------------------
  // 1st argument : Anisotropy matrices 
  T = mxGetPr(prhs[0]);
  //------------------------------------------------------------------
  // 2nd argument : Seeds
  Seeds = (int*)mxGetPr(prhs[1]);
  //------------------------------------------------------------------
  // 3rd and 4th arguments : distance map and label map
  int dimr[2] = {nx,ny};
  if (nrhs >= 3)
    {
      UInit = (double*)mxGetPr(prhs[2]);
      VInit = (int*)mxGetPr(prhs[3]);
      // First output : minimal action map
      //mxSetPr(plhs[0], UInit);
      //mxSetPr(plhs[3], (double*)VInit);
      given_u = true;
    }
  else
    {
      // First output : minimal action map
      plhs[0] = mxCreateNumericArray(2, dimr, mxDOUBLE_CLASS, mxREAL);
      UInit = (double*)mxGetPr(plhs[0]);
      // 4th output : V (Voronoi diagram)
      plhs[1] = mxCreateNumericArray(2, dimr, mxINT32_CLASS, mxREAL);
      VInit = (int*)mxGetPr(plhs[1]);
      given_u = false;
    }
  //------------------------------------------------------------------
  // 5th argument : Seeds (to restrict to a subset of SeedsMatrix)
  if (nrhs >= 5)
    {
      start_points = mxGetPr(prhs[4]); 
      nb_start_points = mxGetN(prhs[4]);
    }
  else
    {
      start_points = NULL;
      nb_start_points = 0;
    }
  //------------------------------------------------------------------
  // 6th argument : Tolerance
  if (nrhs == 6)
    {
      tol = (float)mxGetScalar(prhs[5]);
    }
  else tol = 1e-14;
  //==================================================================
  // Outputs
  int dims[2] = {Nx,Ny};
  if (nrhs >= 3)
    {
      // Second output : dUx
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      dUx = (double*) mxGetPr(plhs[0]);
      // Third output : dUy
      plhs[1] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      dUy = (double*) mxGetPr(plhs[1]);
    }
  else
    {
      // Second output : dUx
      plhs[2] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      dUx = (double*) mxGetPr(plhs[2]);
      // Third output : dUy
      plhs[3] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      dUy = (double*) mxGetPr(plhs[3]);
    }
  //==================================================================
  InitializeNeighborhoods();
  //------------------------------------------------------------------
  InitializeArrays();
  //------------------------------------------------------------------
  InitializeQueue();  
  //------------------------------------------------------------------
  GaussSiedelIterateVoronoi();
  //==================================================================
  resize();
  dims[0] = Nx-2; dims[1] = Ny-2;
  if (nrhs >= 3)
    {
      mxSetDimensions(plhs[0], dims, 2);
      mxSetDimensions(plhs[1], dims, 2);
    }
  else
    {
      mxSetDimensions(plhs[1], dims, 2);
      mxSetDimensions(plhs[2], dims, 2);
    }
  //==================================================================
  DELETEARRAY(NeighborhoodLarge);
  DELETEARRAY(Q);
  DELETEARRAY(Simplicies);
  DELETEARRAY(X1); DELETEARRAY(X2); DELETEARRAY(X12);
  DELETEARRAY(TAB);
  DELETEARRAY(U_n_D);
  DELETEARRAY(Metric);
  DELETEARRAY(Utmp);
  DELETEARRAY(Vtmp);
  delete QSeeds;
  return;
} 