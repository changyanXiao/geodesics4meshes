//================================================================
//================================================================
// File: fm2dSubGradiend.cpp
// (C) 02/2010 by Fethallah Benmansour
//================================================================
//================================================================

#include "fm2dSubGradient.h"

void mexFunction(	int nlhs, mxArray *plhs[], 
					int nrhs, const mxArray*prhs[] ) 
{
	int k;
    //------------------------------------------------------------------
	/* retrive arguments */
	//------------------------------------------------------------------
	if( nrhs!=5 && nrhs!=6 )
		mexErrMsgTxt("5 or 6 arguments are required.");
	if( mxGetNumberOfDimensions(prhs[1])!= 2 )
		mexErrMsgTxt("W must be a 2D double array.");
	//------------------------------------------------------------------
	// first argument : spacing and dimensions
	const int* dim_h = mxGetDimensions(prhs[0]);
    if ( (dim_h[0]!=2) || (dim_h[1]!=1) )
	  mexErrMsgTxt("Library error: h must be a 2x1 array list.");
	hx = mxGetPr(prhs[0])[0]; hy = mxGetPr(prhs[0])[1];
	hx2 = hx*hx; hy2 = hy*hy;
	hx2hy2 = hx*hx*hy*hy;
	hx2_plus_hy2 = hx*hx + hy*hy;
	nx = mxGetDimensions(prhs[1])[0];
	ny = mxGetDimensions(prhs[1])[1];
	Nx = nx+2; Ny = ny+2;
	size = Nx*Ny;
	//------------------------------------------------------------------
	// Second argument : Weight 
	WW = mxGetPr(prhs[1]);
	//------------------------------------------------------------------
	// Third argument : regularization weight 
	w = (float) mxGetScalar(prhs[2]);
	//------------------------------------------------------------------
	// Fourth argument : start points
	start_points = mxGetPr(prhs[3]);
	nb_start_points = mxGetN(prhs[3]);
    if(nb_start_points > 1)
        mexErrMsgTxt("Chose exactely one source point\n");
    //------------------------------------------------------------------
	// Fourth argument : end points
	end_points = mxGetPr(prhs[4]);
	nb_end_points = mxGetN(prhs[4]);
    if(nrhs==6)
        Obstacle = (bool*) mxGetPr(prhs[5]);
    else{
        Obstacle = new bool[size];
        for(k = 0;  k < size; k++)
            Obstacle[k] = false;
    }
	//==================================================================
	// Outputs
	int dims[2] = {Nx,Ny};
	// First output : minimal action map
	plhs[0] = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL );
	U = (float*) mxGetPr(plhs[0]);
	// Second output : Gradient of U / P on point end_point
    int Dims[3] = {Nx,Ny, nb_end_points};
	plhs[1] = mxCreateNumericArray(3, Dims, mxSINGLE_CLASS, mxREAL );
	OutputGradient = (float*) mxGetPr(plhs[1]);
	//==================================================================
	InitializeNeighborhoods();
	//------------------------------------------------------------------
	InitializeArrays();
	//------------------------------------------------------------------
	InitializeOpenHeap();
	//------------------------------------------------------------------
	RunPropagation();
    //==================================================================
    resize();
    dims[0] = Nx-2; dims[1] = Ny-2;
    Dims[0] = Nx-2; Dims[1] = Ny-2;
    mxSetDimensions(plhs[0], dims, 2);
    mxSetDimensions(plhs[1], Dims, 3);
    //==================================================================
	DELETEARRAY(S);
	DELETEARRAY(W);
    cleanGradients();
    DELETEARRAY(Gradients);
    DELETEARRAY(dUx); DELETEARRAY(dUy);
    DELETEARRAY(Trial);
    DELETEARRAY(Tree);
    if(nrhs==5)
        DELETEARRAY(Obstacle);
    return;
}