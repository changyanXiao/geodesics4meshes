//================================================================
//================================================================
// File: AnisoEikonalSolver2D.cpp
// (C) Fethallah Benmansour 2010 -- CVLab-EPFL --
//================================================================
//================================================================

#include "AnisoEikonalSolver2D.h"

void mexFunction(	int nlhs, mxArray *plhs[], 
					int nrhs, const mxArray*prhs[] ) 
{
	//------------------------------------------------------------------
	/* retrive arguments */
	//------------------------------------------------------------------
	if( (nrhs!=3) && (nrhs!=4) )
		mexErrMsgTxt("3 or 4 arguments are required.");
	if( mxGetNumberOfDimensions(prhs[1])!= 4 )
		mexErrMsgTxt("HERE T must be a 2D x 2x2 tensor field of symmetric definite matrices.");
	//------------------------------------------------------------------
	// First argument : spacing and dimensions
	const int* dim_h = mxGetDimensions(prhs[0]);
    if ( (dim_h[0]!=2) || (dim_h[1]!=1) )
	  mexErrMsgTxt("Library error: h must be a 2x1 array list.");
	hx = mxGetPr(prhs[0])[0]; hy = mxGetPr(prhs[0])[1];
	hx2 = hx*hx; hy2 = hy*hy;
    hxhy = hx*hy;
	nx = mxGetDimensions(prhs[1])[0];
	ny = mxGetDimensions(prhs[1])[1];
    if( (mxGetDimensions(prhs[1])[2] != 2) || (mxGetDimensions(prhs[1])[3] != 2) )
        mexErrMsgTxt("T must be a 2D x 2x2 tensor field of symmetric definite matrices.");
	Nx = nx+2; Ny = ny+2;
	size = Nx*Ny; nxny = nx*ny;
	//------------------------------------------------------------------
	// Second argument : Anisotropy matrices 
	T = mxGetPr(prhs[1]);
	//------------------------------------------------------------------
	// Third argument : start points
	start_points = mxGetPr(prhs[2]);
	nb_start_points = mxGetN(prhs[2]);
	//------------------------------------------------------------------
    // forth argument : Tolerance
    if(nrhs==4)
    	tol = (float) mxGetScalar(prhs[3]);
    else
        tol = 1e-8;
	//==================================================================
	// Outputs
	int dims[2] = {Nx,Ny};
	// First output : minimal action map
	plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL );
	U = (double*) mxGetPr(plhs[0]);
	// Second output : dUx
	plhs[1] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL );
	dUx = (double*) mxGetPr(plhs[1]);
	// Third output : dUy
	plhs[2] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL );
	dUy = (double*) mxGetPr(plhs[2]);
    // Fourth output : V
    plhs[3] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL );
	V = (double*) mxGetPr(plhs[3]);
	//==================================================================
	InitializeNeighborhoods();
	//------------------------------------------------------------------
	InitializeArrays();
	//------------------------------------------------------------------
	InitializeQueue();
	//------------------------------------------------------------------
	GaussSiedelIterate();
    //==================================================================
    resize();
    dims[0] = Nx-2; dims[1] = Ny-2;
    mxSetDimensions(plhs[0], dims, 2);
    mxSetDimensions(plhs[1], dims, 2);
    mxSetDimensions(plhs[2], dims, 2);
    mxSetDimensions(plhs[3], dims, 2);
    //==================================================================
    DELETEARRAY(NeighborhoodLarge);
    DELETEARRAY(S);
    DELETEARRAY(Q);
    DELETEARRAY(ITER);
    DELETEARRAY(Simplicies);
    DELETEARRAY(X1); DELETEARRAY(X2); DELETEARRAY(X12);
    DELETEARRAY(TAB);
    DELETEARRAY(U_n_D);
	DELETEARRAY(Metric);
    return;
}