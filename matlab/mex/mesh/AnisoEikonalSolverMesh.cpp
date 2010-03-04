//================================================================
//================================================================
// File: AnisoEikonalSolverMesh.cpp
// (C) 02/2010 by Fethallah Benmansour & Gabriel Peyrée
//================================================================
//================================================================

#include "AnisoEikonalSolverMesh.h"

void mexFunction(	int nlhs, mxArray *plhs[], 
				 int nrhs, const mxArray*prhs[] ) 
{
    //==================================================================
	/* retrive arguments */
	if( nrhs<4 ) 
		mexErrMsgTxt("4 or 5 input arguments are required.");
	if( nlhs<1  || nlhs>5 ) 
		mexErrMsgTxt("1 up to 5 output arguments are required.");
    //==================================================================
	// arg1 : vertex
	vertex = mxGetPr(prhs[0]);
	nverts = mxGetN(prhs[0]); 
	if( mxGetM(prhs[0])!=3 )
		mexErrMsgTxt("vertex must be of size 3 x nverts."); 
	//==================================================================
    // arg2 : faces
	faces = mxGetPr(prhs[1]);
	nfaces = mxGetN(prhs[1]);
	if( mxGetM(prhs[1])!=3 )
		mexErrMsgTxt("face must be of size 3 x nfaces."); 
	//==================================================================
    // arg3 : Metric (should be symmetric definite positive on the tangent plan)
    // order of the 6 components for xx, yy, zz, xy, yz, zx
	if( (mxGetDimensions(prhs[2])[0] != 6) || (mxGetDimensions(prhs[2])[1] != nverts) )
        mexErrMsgTxt("T must be of same size as vertex with 6 components : 6xnverts.");
	T = mxGetPr(prhs[2]);
    //==================================================================
	// arg4 : start_points
	start_points = mxGetPr(prhs[3]);
	nstart = mxGetM(prhs[3]);
	//==================================================================	
	// argument 5 and 6: if initial distance values at source points are given
    // In this case, provide the intial voronoi indices as well
	if(nrhs==5)
        mexErrMsgTxt("4 or 6 argumets required");
    if( nrhs==6 )
	{
		U_init_values = mxGetPr(prhs[4]);
        V_init_values = mxGetPr(prhs[5]);
		if( mxGetM(prhs[4])==0 && mxGetN(prhs[4])==0 )
			U_init_values=NULL;
        if( mxGetM(prhs[5])==0 && mxGetN(prhs[5])==0 )
			V_init_values=NULL;
		if( U_init_values!=NULL && (mxGetM(prhs[4])!=nstart || mxGetN(prhs[4])!=1) )
			mexErrMsgTxt("values must be of size nb_start_points x 1."); 
        if( V_init_values!=NULL && (mxGetM(prhs[5])!=nstart || mxGetN(prhs[5])!=1) )
			mexErrMsgTxt("values must be of size nb_start_points x 1.");
	}
	else{
		U_init_values = NULL;
        V_init_values = NULL;
	}
	//==================================================================
	// first ouput : geodesic distance
	plhs[0] = mxCreateNumericArray(1,&nverts, mxDOUBLE_CLASS, mxREAL ); 
	U = (double*) mxGetPr(plhs[0]);
    // second output : voronoi
	plhs[1] = mxCreateNumericArray(1,&nverts, mxDOUBLE_CLASS, mxREAL ); 
	Vor = (double*) mxGetPr(plhs[1]);
    // directions of characteristics
	plhs[2] = mxCreateNumericArray(1,&nverts, mxDOUBLE_CLASS, mxREAL ); 
	dUx = (double*) mxGetPr(plhs[2]);
    plhs[3] = mxCreateNumericArray(1,&nverts, mxDOUBLE_CLASS, mxREAL ); 
	dUy = (double*) mxGetPr(plhs[3]);
    plhs[4] = mxCreateNumericArray(1,&nverts, mxDOUBLE_CLASS, mxREAL ); 
	dUz = (double*) mxGetPr(plhs[4]);
	//==================================================================
    create_the_mesh();
    //------------------------------------------------------------------
    InitializeArrays();
	//------------------------------------------------------------------
	InitializeQueue();
    //------------------------------------------------------------------
	GaussSiedelIterate();
	//==================================================================
    DELETEARRAY(S);
    DELETEARRAY(Q);
    DELETEARRAY(ITER);
    DELETEARRAY(TAB);
    DELETEARRAY(U_n_D);
};