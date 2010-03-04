//================================================================
//================================================================
// File: AnisoEikonalSolver2D.h
// (C) Fethallah Benmansour 2010 -- CVLab-EPFL --
//================================================================
//================================================================
#include "fm.h"
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <functional>
#include <queue>

using namespace std;
// Define a template class vector of int
typedef vector<int, allocator<int> > StrVector ;

//Define an iterator for template class vector of strings
typedef StrVector::iterator StrVectorIt ;

#define isEqual(a, b) (a==b ? 1 : 0);
#define kDead -1 // only for source points
#define kEstimated -2
#define kFar -3
#define kBorder -4
#define kEnqueued 0
#define kUnqueued 1

queue<int> waiting_Points;
/* Global variables */
int nx;			// real size on X
int ny;			// real size on Y
int nxny;
int Nx, Ny; // size for computing
int size;
double hx, hy;// spacing
double hx2, hy2, hxhy;
double*   U  = NULL;// action map
double*  dUx = NULL; // X direction of characteristic 
double*  dUy = NULL; // Y direction of characteristic 
double*  V   = NULL; // Voronoi
short*   S   = NULL; // states
short*   Q   = NULL; // Queue
double* Metric = NULL; // Metric
double* T = NULL;
int* ITER = NULL; // Number of iterations per point
//================================================================
int   connectivity_large;
int*    NeighborhoodLarge = NULL;
int NB_simplicies;
int* Simplicies = NULL;
double *X1 = NULL, *X2 = NULL, *X12 = NULL;
double *TAB = NULL;
double *U_n_D = NULL;
//================================================================
double w; // regularisation parameter
double Dmax;
double* start_points = NULL;
int nb_start_points;
double tol = 0.1;

//================================================================
void InitializeNeighborhoods()
//================================================================
{
	//-----------------------------------------------------------
	connectivity_large = 8;
	NeighborhoodLarge = new int[connectivity_large];
    int i,j,l,pos=0;
    for(i = -1; i <=1; i++){
	for(j = -1; j <=1; j++){
        if( (j!=0) || (i!=0) )
            NeighborhoodLarge[pos++]= j*Nx+i;
    }}
    //-----------------------------------------------------------
    NB_simplicies = 8;
    Simplicies  = new int  [2*NB_simplicies];
    X12         = new double[2*NB_simplicies];
    X1          = new double[2*NB_simplicies];
    X2          = new double[2*NB_simplicies];
    //-----------------------------------------------------------
    const int NB_Dimension = 2;
    StrVector Pattern(NB_Dimension);
    StrVectorIt start, end;
    start = Pattern.begin() ;   // location of first element of Pattern
    end = Pattern.end() ;       // one past the location last element of Pattern
    pos=0;
    //Initialize vector Pattern
    Pattern[0] = 1;
    Pattern[1] = Nx;
    int reff[2]={1,Nx};
    // Generate all possible permutations
    int* sign = new int[2];
    do
    {
        for(i = -1; i <=1; i=i+2){
        for(j = -1; j <=1; j=j+2){
            sign[0] = i; sign[1] = j;
            for(l=0; l<2; l++){
                Simplicies[2*pos+l] =              sign[l]*start[l];
                X1 [2*pos+l]        =             -sign[0]*isEqual(start[0], reff[l]);
                X2 [2*pos+l]        = X1[2*pos+l] -sign[1]*isEqual(start[1], reff[l]);
                X12[2*pos+l]        = X1[2*pos+l] - X2[2*pos+l];
            }
            pos++;
        }}
    }
    while ( next_permutation(start, end) );
    Pattern.clear();
    DELETEARRAY(sign);
}

//================================================================
void InitializeArrays()
//================================================================
{
	int x, y, point;
	//copy the weight list and initialize
	S = new short[size];
    Q = new short[size];
    ITER = new int[size];
    Metric = new double[3*size];
    TAB    = new double[5];
    U_n_D  = new double[4];
    //------------------------------------------------------------
    for(x = 0; x < nx; x++){
		for(y = 0; y < ny; y++){
			point = (x+1) + (y+1)*Nx;
			Metric[3*point    ] = hx2* T[x + y*nx];
            Metric[3*point + 1] = hy2* T[x + y*nx + 3*nxny];
            Metric[3*point + 2] = hxhy*T[x + y*nx + nxny];
			S[point] = kFar;
            Q[point] = kUnqueued;
		}
	}
	for(x = 0; x < size; x++){
		U[x] = INFINITE;
        V[x] = kBorder;
        ITER[x] = 0;
	}
	//------------------------------------------------------------
	// Initialize Borders
	for(x = 0; x < Nx; x++){
		y = 0;
		point = x + y*Nx;
		S[point] = kBorder;
		y = Ny-1;
		point = x + y*Nx;
		S[point] = kBorder;
	}
	for(y = 0; y < Ny; y++){
		x = 0;
		point = x + y*Nx;
		S[point] = kBorder;
		x = Nx-1;
		point = x + y*Nx;
		S[point] = kBorder;
	}
}

//================================================================
void InitializeQueue()
//================================================================
{
	int point, npoint, i, j, k, s;

    for( s=0; s<nb_start_points; ++s ){
		i = round(start_points[2*s]);
		j = round(start_points[1+2*s]);
		point = i + j*Nx;
		//--------------------------------------------------------
		if(point >=size)
			mexErrMsgTxt("start_points should in the domaine.");
		//--------------------------------------------------------
		if( U[point]==0 )
			mexErrMsgTxt("start_points should not contain duplicates.");
		//--------------------------------------------------------
		U[point] = 0.0; S[point] = kDead; V[point] = s;
		// add to heap
        for(k = 0; k < connectivity_large; k++){
            npoint = point + NeighborhoodLarge[k];
            if( (S[npoint] != kBorder) && (S[npoint] != kEstimated) )
                waiting_Points.push(npoint);
                Q[npoint] = kEnqueued;
        }
		//--------------------------------------------------------
	}
}

//================================================================
void DotProductMetric(double *res, double *V1, double *V2, double *M)
/*
COMMENTS : 
 * compute V1^t M V2
 * M is symetric
*/
//================================================================
{
    *res    = M[0]*V1[0]*V2[0] + M[1]*V1[1]*V2[1] 
            + M[2]*( V1[0]*V2[1] + V2[0]*V1[1] );
}

//================================================================
void TsitsiklisNorm(double* res, double* M, double alpha,
                    double* y1, double* y2)
/*
COMMENTS : 
 * compute M^{-1} norm of vector (alpha*y1 + y2)
*/
//================================================================
{
    double DetM;
    double x0, x1;
    x0 = alpha*y1[0] + y2[0];
    x1 = alpha*y1[1] + y2[1];
    DetM = M[0]*M[1] -M[2]*M[2];
    *res = M[1]*x0*x0 + M[0]*x1*x1 - 2.0*M[2]*x0*x1;
    *res /= DetM;
    *res = sqrt(*res);
}


//================================================================
void TsitsiklisOnePoint(double *res, double* M, double u, double* Z)
//================================================================
{
    double NormZ_M;
    DotProductMetric(&NormZ_M, Z, Z, M);
    *res = u + sqrt(NormZ_M);
}

//================================================================
void TsitsiklisTwoPoints(double* res, double* M, double  k, double u, double* z1, double* z2)
/*
COMMENTS : 
 * computes the minimum and the arg min of :
 * \alpha*k + u + \| \alpha* z_1 + z_2 \|_M ; \alpha \in [0, 1]
 * res is a double array of size 2 containing the argmin and the min
*/
//================================================================
{
    double r11, r22, r12, R;
    DotProductMetric(&r11, z1, z1, M);
    DotProductMetric(&r22, z2, z2, M);
    DotProductMetric(&r12, z1, z2, M);
    R = r11*r22-r12*r12;//since M is symetric definite positive, R > 0.
    if( k >= sqrt(r11) ){
        res[0] = 0.0;
        res[1] = u + sqrt(r22);
    }
    else if(k <= -sqrt(r11) ){
        res[0] = 1.0;
        res[1] = k + u + sqrt(r11 + r22 + 2.0*r12);
    }
    else{
        if( r12>= -k*sqrt(R/(r11-k*k)) ){
            res[0] = 0.0;
            res[1] = u + sqrt(r22);
        }
        else if( r12 <= (-r11-k*sqrt(R/(r11-k*k)))  ){
            res[0] = 1.0;
            res[1] = k + u + sqrt(r11 + r22 + 2.0*r12);
        }
        else{
            res[0] = -(r12 + k*sqrt(R/(r11-k*k))) / r11;
            res[1] = res[0]*k + u + sqrt(R/(r11-k*k));
        }
    }
}

//================================================================
void TsitsiklisQuadrant(double* res, int point, int simplex_idx)
//================================================================
{
    int	  npoint1, npoint2, i;
    double k, u;
    double* M;
    double U1, U2, V1, V2;
	double Ur = U[point];
	double *y1, *y2;
    bool b_point1, b_point2;
    //--------------------------------------------------------------
    M = Metric + 3*point;
    //--------------------------------------------------------------
    res[1] = Ur+1.0;// assign bigger value
    //--------------------------------------------------------------
    npoint1 =  point  + Simplicies[2*simplex_idx];
    npoint2 = npoint1 + Simplicies[2*simplex_idx+1];
    //----------------------------------------------------------
	b_point1 = (S[npoint1]==kEstimated) || (S[npoint1]==kDead);
    b_point2 = (S[npoint2]==kEstimated) || (S[npoint2]==kDead);
    if(b_point1){  U1 = U[npoint1]; V1 = V[npoint1];}
    else        {  U1 = INFINITE;   V1 = kBorder;}
    if(b_point2){  U2 = U[npoint2]; V2 = V[npoint2];}
    else        {  U2 = INFINITE;   V2 = kBorder;}
    //----------------------------------------------------------

    if(b_point1 && b_point2){
        y1 = X12+2*simplex_idx; y2 = X2+2*simplex_idx;
        k = U1-U2; u = U2;
        TsitsiklisTwoPoints(res, M, k, u, y1, y2);
        for(i = 2; i < 4; i++){
            res[i] = (res[0]*y1[i-2] + y2[i-2]);
        }
        res[4] = res[0]*V1 +(1-res[0])*V2;
    }
    //----------------------------------------------------------
    else if(b_point1 && !b_point2 ){
        y1 = X1+2*simplex_idx; y2 = X1+2*simplex_idx;
        u = U1;
        TsitsiklisOnePoint(res, M, u, y1);
        res[1] = res[0];
        res[0] = 1.0;
        for(i = 2; i < 4; i++){
            res[i] = y2[i-2];
        }
        res[4] = V1;
    }
    //----------------------------------------------------------
    else if(!b_point1 && b_point2 ){
        y1 = X2+2*simplex_idx; y2 = X2+2*simplex_idx;
        u = U2;
        TsitsiklisOnePoint(res, M, u, y1);
        res[1] = res[0];
        res[0] = 0;
        for(i = 2; i < 4; i++){
            res[i] = y2[i-2];
        }
        res[4] = V2;
    }
};

//================================================================
void TsitsiklisUpdate(int point, double* res)
/*
COMMENTS : 
*/
//================================================================
{
	int i;
    bool is_updated = false;
    double *mat;
    double Ur = U[point], dUxr, dUyr, Vr;
    mat = Metric + 3*point;
    //--------------------------------------------------------------
    for(i=0; i<NB_simplicies; i++){
		//----------------------------------------------------------
        TsitsiklisQuadrant(TAB, point, i);
        //----------------------------------------------------------
        if (TAB[1]<Ur){
            Ur   = TAB[1];
            dUxr = TAB[2]; dUyr = TAB[3]; Vr = TAB[4];
            is_updated = true;
    	}
	}
	//--------------------------------------------------------------
	if (is_updated){
		res[0] = Ur;
        res[1] = dUxr;
        res[2] = dUyr;
        res[3] = Vr;
	}
    else{// copy old values
		res[0] = U[point];
        res[1] = dUx[point];
        res[2] = dUy[point];
        res[3] = V[point];
	}
    //--------------------------------------------------------------
}

//================================================================
void GaussSiedelIterate()
//================================================================
{
    //------------------------------------------------------------
    int point,npoint,k;
    int iter = 0;
	//------------------------------------------------------------
    while( !waiting_Points.empty()  ) {
        point = waiting_Points.front();
        waiting_Points.pop();
        TsitsiklisUpdate(point, U_n_D);
        Q[point] = kUnqueued;        
        S[point] = kEstimated; 
        ITER[point] = ITER[point] + 1;
        iter = MAX(iter, ITER[point]);
        if(fabs(U_n_D[0] - U[point]) > tol){
            U[point]   = U_n_D[0];
            dUx[point] = U_n_D[1];
            dUy[point] = U_n_D[2];
            V[point]   = U_n_D[3];
            for(k = connectivity_large-1; k >=0; k--){
                npoint = point + NeighborhoodLarge[k];
                if( (S[npoint] != kBorder) && (S[npoint] != kDead) && (Q[npoint] != kEnqueued)){//
                    waiting_Points.push(npoint);
                    Q[npoint] = kEnqueued;
                }
            }
       }
    }
    mexPrintf("max iter for a point = %d\n", iter);
}

//================================================================
void resize()
//================================================================
{
    int x, y, point, Point;
    for(y=0;y<ny;y++)
        for(x=0;x<nx;x++){
            point = x+y*nx;
            Point = (x+1)+(y+1)*Nx;
            U[point]   = U[Point];
            dUx[point] = dUx[Point];
            dUy[point] = dUy[Point];
            V[point]   = V[Point];
        }
}
