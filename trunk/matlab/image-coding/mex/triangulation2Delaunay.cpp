// tu peux trouver en pièce jointe le code pour les flips. Pour l'utiliser :
//   [LF,E] = triangulation2Delaunay(faces,Pts);
// où LF est la liste des flips (du 1er au dernier flip effectué), il faut donc l'inverser pour retrouver la triangulation initale à partir de celle de Delaunay.
//  E est la liste des arêtes de Delaunay
//   faces est la liste des triangles représentés par les indices (entiers) des points Pts (à partir de 1 obligatoirement).

// File: triangulation2Delaunay.cpp
// author: Sebastien Bougleux (sebastien.bougleux@unicaen.fr),
//         Universite de Caen, GREYC
// date: 10 feb. 2009

#include<mex.h>
#include "triangulationDS.h"

void mexFunction(int nlhs, mxArray *plhs[], 
				 int nrhs, const mxArray*prhs[] ) 
{
	//------------------------------------------------------------------
	/* retrive arguments */
	//------------------------------------------------------------------
	if (nrhs != 2) mexErrMsgTxt("2  arguments are required.");
	//------------------------------------------------------------------
	// First argument : Faces
	int *faces = (int*)mxGetPr(prhs[0]);
	int nbdim = mxGetDimensions(prhs[0])[0];
	int nbf = mxGetDimensions(prhs[0])[1];
    //------------------------------------------------------------------
	// First argument : Points
	double *pts = mxGetPr(prhs[1]);
	int nbdim2 = mxGetDimensions(prhs[1])[0];
	int nbp = mxGetDimensions(prhs[1])[1];
    //==================================================================
	Triangulation<int,double> *T = new Triangulation<int,double>;
	T->loadFacesAndFlip(faces, nbf*nbdim, pts, nbp);
    //=================================================>=================
    // Outputs
    // First output : edges of the triangulation
	int dims[2] = { 2, T->nbFlipEdges() };
	plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
	double *edges_flip = (double*)mxGetPr(plhs[0]);
    T->getFlipList(edges_flip);
	// 2nd output : edges of the Delaunay triangulation
    if (nlhs > 1)
    {
        int dime[2] = { 2, T->nbEdges() };
        plhs[1] = mxCreateNumericArray(2, dime, mxDOUBLE_CLASS, mxREAL);
        double *edges_del = (double*)mxGetPr(plhs[1]);
        T->getEdges(edges_del);
    }
	//=================================================>=================
    delete T;
    //=================================================>=================
    return;
}
