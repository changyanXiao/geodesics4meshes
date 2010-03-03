#include "fm.h"
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <functional>
#include <queue>
#include <map>
#include<cmath>
using namespace std;
// Define a template class vector of int
typedef vector<int, allocator<int> > StrVector ;

//Define an iterator for template class vector of strings
typedef StrVector::iterator StrVectorIt ;

#define isEqual(a, b) (a==b ? 1 : 0);
#define kBorder -3
#define kUnqueued -2
#define kEnqueued -1
#define kSeed 0
#define kEnqueuedEstimated 1
#define kUnqueuedEstimated 2

std::queue<int> waiting_Points;
/* Global variables */
int nx;			// real size on X
int ny;			// real size on Y
int nxny;
int Nx, Ny; // size for computing
int size;
//double hx, hy;// spacing
//double hx2, hy2, hxhy;
//double* U = NULL;// action map
double* UInit = NULL;// action map (given)
double*  dUx = NULL;
double*  dUy = NULL;
//int* V = NULL;  // Voronoi
int* VInit = NULL;  // Voronoi (given if UInit is given)
int* Seeds = NULL;  // Labels (given if UInit is given)
short* Q = NULL;
double* Metric = NULL; // Metric
double* T = NULL;
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
int nb_start_points = 0;
double tol = 0.1;
bool given_u = false;
std::queue<int> *QSeeds = NULL;

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
    }
  }
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
        }
      }
    }
  while ( next_permutation(start, end) );
  Pattern.clear();
  DELETEARRAY(sign);
}

//================================================================
void InitializeArrays()
//================================================================
{
  int x, y, point, Point, v;
  Q = new short[size];
  Metric = new double[3*size];
  TAB    = new double[5];
  U_n_D  = new double[4];
  QSeeds = new std::queue<int>;
  //------------------------------------------------------------
  if (given_u == false)  // no given distance map (and label map)
    {
      for (x = 0; x < nx; x++)
	{
	  for (y = 0; y < ny; y++)
	    {
	      Point = x + y*nx;           // original matrix
	      point = (x+1) + (y+1)*Nx;   // extended matrix
	      // metric computed on extended matrix (TO MODIFY)
	      Metric[3*point    ] = T[Point];
	      Metric[3*point + 1] = T[Point + 3*nxny];
	      Metric[3*point + 2] = T[Point + nxny];
	      //check domain constraints
	      v = Seeds[Point];
	      if (v < 0)  // not in the domain (mask)
		{
		  UInit[Point] = -1;
		  VInit[Point] = -1;
		  Q[point] = kBorder;
		  continue;
		}
	      if (v > 0)  // a seed
		{
		  UInit[Point] = 0.0;
		  VInit[Point] = v;  // take its label
		  Q[point] = kSeed;
		  QSeeds->push(point);
		  continue;
		}
	      // a point in the domain not yet processed
	      UInit[Point] = INFINITE;
	      VInit[Point] = 0;
	      Q[point] = kUnqueued;
	    }
	}
    }
  else // a given distance map (and Voronoi map)
    {
      for (x = 0; x < nx; x++)
	{
	  for (y = 0; y < ny; y++)
	    {
	      Point = x + y*nx;           // original matrix
	      point = (x+1) + (y+1)*Nx;   // extended matrix
	      // metric computed on extended matrix (TO MODIFY)
	      Metric[3*point    ] = T[Point];
	      Metric[3*point + 1] = T[Point + 3*nxny];
	      Metric[3*point + 2] = T[Point + nxny];
	      //check domain constraints
	      v = Seeds[Point];
	      if (v < 0)  // not in the domain (mask)
		{
		  Q[point] = kBorder;
		  continue;
		}
	      if (v > 0)
		{
		  Q[point] = kSeed;
		  QSeeds->push(point);
		  continue;
		}
	      if (VInit[Point] > 0) // already processed
		{
		  Q[point] = kUnqueuedEstimated;
		}
	      else Q[point] = kUnqueued;
	    }
	}
    }
  //------------------------------------------------------------
  // Initialize Borders
  for (x = 0; x < Nx; x++)
    {
      y = 0;
      point = x + y*Nx;
      Q[point] = kBorder;
      //U[point] = INFINITE;
      //V[point] = -1;
      y = Ny-1;
      point = x + y*Nx;
      Q[point] = kBorder;
      //U[point] = INFINITE;
      //V[point] = -1;
    }
  for (y = 0; y < Ny; y++)
    {
      x = 0;
      point = x + y*Nx;
      Q[point] = kBorder;
      //U[point] = INFINITE;
      //V[point] = -1;
      x = Nx-1;
      point = x + y*Nx;
      Q[point] = kBorder;
      //U[point] = INFINITE;
      //V[point] = -1;
    }
  
}

//================================================================
void pos_to_point(const int &pos, int &x, int &y)
//================================================================
{
  y = pos / Nx;
  x = pos - Nx * y;
}

//================================================================
int comp_to_real(const int &pos)
//================================================================
{
  int y = pos / Nx;
  return (pos-Nx*y-1) + (y - 1)*nx;
}

//================================================================
void InitializeQueue()
//================================================================
{
  int point, npoint, i, j, k, s;
  short q;
  //--------------------------------------------------------
  if (nb_start_points == 0)
    {
      while (!QSeeds->empty())
	{
	  point = QSeeds->front();
	  QSeeds->pop();
	  // add to heap
	  for (k = 0; k < connectivity_large; k++)
	    {
	      npoint = point + NeighborhoodLarge[k];
	      q = Q[npoint];
	      if (q == kUnqueued)
		{
		  waiting_Points.push(npoint);
		  Q[npoint] = kEnqueued;
		  continue;
		}
	      if (q == kUnqueuedEstimated)
		{
		  waiting_Points.push(npoint);
		  Q[npoint] = kEnqueuedEstimated;
		}
	    }
	}
      return;
    }
  //--------------------------------------------------------
  // a subset of Seeds is given
  for (s = 0; s < nb_start_points; ++s) 
    {
      i = (int)round(start_points[2*s]);
      j = (int)round(start_points[1+2*s]);
      point = i + j*Nx;
      //--------------------------------------------------------
      if (point >= size)
	{
	  mexErrMsgTxt("Seeds should be in the domain.");
	}
      //--------------------------------------------------------
      if (Q[point] != kSeed)
	{
	   mexErrMsgTxt("A seed has no label");
	}
      //--------------------------------------------------------
      for (k = 0; k < connectivity_large; k++)
	{
	  npoint = point + NeighborhoodLarge[k];
	  q = Q[npoint];
	  if (q == kUnqueued)
	    {
	      waiting_Points.push(npoint);
	      Q[npoint] = kEnqueued;
	      continue;
	    }
	  if (q == kUnqueuedEstimated)
	    {
	      waiting_Points.push(npoint);
	      Q[npoint] = kEnqueuedEstimated;
	    }
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
  *res = M[0]*V1[0]*V2[0] + M[1]*V1[1]*V2[1] + M[2]*(V1[0]*V2[1] + V2[0]*V1[1]);
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
void TsitsiklisTwoPoints(double* res, double* M, const double &k, 
			 const double &u, double* z1, double* z2,
			 double V1, double V2)
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
  R = r11*r22-r12*r12;     //since M is symetric definite positive, R > 0.
  double sr11 = sqrt(r11);
  if (k >= sr11)
    {
      res[0] = 0.0;
      res[1] = u + sqrt(r22);
      res[4] = V2;
      return;
    }
  if (k <= -sr11)
    {
      res[0] = 1.0;
      res[1] = k + u + sqrt(r11 + r22 + 2.0*r12);
      res[4] = V1;
      return;
    }
  double sR = sqrt(R/(r11-k*k));
  if (r12 >= -k*sR)
    {
      res[0] = 0.0;
      res[1] = u + sqrt(r22);
      res[4] = V2;
      return;
    }
  if (r12 <= (-r11-k*sR))
    {
      res[0] = 1.0;
      res[1] = k + u + sqrt(r11 + r22 + 2.0*r12);
      res[4] = V1;
      return;
    }
  res[0] = -(r12 + k*sR) / r11;
  res[1] = res[0]*k + u + sR;
  // modif : compute new Voronoi label from V1 and V2
  if (V1 == V2)  // case same real label
    {
      res[4] = V1;  // or V2
      return;
    }
  res[4] = (res[0] > 0.5 ? V1 : V2);
}

//================================================================
void TsitsiklisQuadrant(double* res, int point, int simplex_idx)
//================================================================
{
  int npoint1, npoint2, i, simpIdx = 2*simplex_idx, nPoint1, nPoint2;
  double k, u, Norm;
  double* M;
  double U1, U2;
  double V1, V2; // modif
  double Ur = UInit[comp_to_real(point)];
  double *y1, *y2;
  bool b_point1, b_point2;
  //--------------------------------------------------------------
  M = Metric + 3*point;
  //--------------------------------------------------------------
  res[1] = Ur + 1.0;// assign bigger value
  //--------------------------------------------------------------
  npoint1 =  point  + Simplicies[simpIdx];
  npoint2 = npoint1 + Simplicies[simpIdx+1];
  nPoint1 = comp_to_real(npoint1);
  nPoint2 = comp_to_real(npoint2);
  //----------------------------------------------------------
  b_point1 = (Q[npoint1] >= 0);
  b_point2 = (Q[npoint2] >= 0);
  if (b_point1)  { U1 = UInit[nPoint1]; V1 = VInit[nPoint1]; }
  else           { U1 = INFINITE;  V1 = -1; }
  if (b_point2)  { U2 = UInit[nPoint2]; V2 = VInit[nPoint2]; }
  else           { U2 = INFINITE; V2 = -1; }
  //----------------------------------------------------------
  if (b_point1 && b_point2)
    {
      y1 = X12 + simpIdx;
      y2 = X2 + simpIdx;
      k = U1-U2; u = U2;
      //int v = 
      TsitsiklisTwoPoints(res, M, k, u, y1, y2, V1, V2);
      TsitsiklisNorm(&Norm, M, res[0], y1, y2);
      for(i = 2; i < 4; i++) res[i] = (res[0]*y1[i-2] + y2[i-2])/Norm;
      return;
    }
  //----------------------------------------------------------
  if (b_point1 && !b_point2 )
    {
      //std::cerr << "yes1 " << U1 << " " << U2 << "\n";
      y1 = X1 + simpIdx;
      y2 = X1 + simpIdx;
      u = U1;
      TsitsiklisOnePoint(res, M, u, y1);
      res[1] = res[0];
      res[0] = 1.0;
      TsitsiklisNorm(&Norm, M, 0, y1, y2);
      for(i = 2; i < 4; i++) res[i] = y2[i-2]/Norm;
      // compute Voronoi index
      res[4] = V1;
      return;
    }
  //----------------------------------------------------------
  if (!b_point1 && b_point2 )
    {
      //std::cerr << "yes2 " << U1 << " " << U2 << "\n";
      y1 = X2 + simpIdx;
      y2 = X2 + simpIdx;
      u = U2;
      TsitsiklisOnePoint(res, M, u, y1);
      res[1] = res[0];
      res[0] = 0;
      TsitsiklisNorm(&Norm, M, 0, y1, y2);
      for(i = 2; i < 4; i++) res[i] = y2[i-2]/Norm;
      // compute Voronoi index
      res[4] = V2;
    }
  // error !b_point1 && !b_point2 (cannot arrive)
}

//================================================================
bool TsitsiklisUpdate(int point, double* res)
/*
COMMENTS : modified
*/
//================================================================
{
  int i;
  bool is_updated = false;
  double Norm;
  double *mat, det;
  int Point = comp_to_real(point);
  //double Ur = U[point], dUxr, dUyr ;
  double Ur = UInit[Point], dUxr, dUyr ;
  double Vres = -1;
  mat = Metric + 3*point;
  //std::cerr << Ur << std::endl;
  //--------------------------------------------------------------
  for (i = 0; i < NB_simplicies; i++)
    {
      //----------------------------------------------------------
      TsitsiklisQuadrant(TAB, point, i);
      //----------------------------------------------------------
      //std::cerr << TAB[1] << std::endl;
      if (TAB[1] < Ur)
	{
	  Ur   = TAB[1];
	  dUxr = TAB[2];
	  dUyr = TAB[3];
	  Vres = TAB[4];
	  is_updated = true;
	}
    }
  //--------------------------------------------------------------
  if (is_updated)
    {
      res[0] = Ur;
      res[1] = dUxr;
      res[2] = dUyr;
      res[3] = Vres;
    }
  //--------------------------------------------------------------
  //std::cerr << is_updated << std::endl;
  return is_updated;
}

//================================================================
void GaussSiedelIterateVoronoi() // modified
//================================================================
{
  int point, npoint, k, v, Point;
  bool is_updated = false;
  double mx = -1;
  int pmx;
  short q;
  //------------------------------------------------------------
  while (!waiting_Points.empty()) 
    {
      point = waiting_Points.front();
      waiting_Points.pop();
      is_updated = TsitsiklisUpdate(point, U_n_D);
      Q[point] = kUnqueuedEstimated;
      if (is_updated)
	{
	  Point = comp_to_real(point);
	  if (fabs(U_n_D[0] - UInit[Point]) > tol)
	    {                
	      UInit[Point] = U_n_D[0];
	      dUx[point]   = U_n_D[1];
	      dUy[point]   = U_n_D[2];
	      VInit[Point] = (int)U_n_D[3];
	      for (k = 0; k < connectivity_large; k++)
		{
		  npoint = point + NeighborhoodLarge[k];
		  q = Q[npoint];
		  if (q == kUnqueued)
		    {
		      waiting_Points.push(npoint);
		      Q[npoint] = kEnqueued;
		      continue;
		    }
		  if (q == kUnqueuedEstimated)
		    {
		      waiting_Points.push(npoint);
		      Q[npoint] = kEnqueuedEstimated;
		    }
		} 
	    }             
	  else
	    {
	      //if (mx < UInit[point-(Nx+1)]) { mx = UInit[point-(Nx+1)]; pmx = point; } 
	    }
	}
    }
  //return pmx;
}

//================================================================
void resize()
//================================================================
{
  int x, y, point, Point;
  for (y = 0; y < ny; y++)
    for (x = 0; x < nx; x++)
      {
	point = x+y*nx;
	Point = (x+1)+(y+1)*Nx;
	dUx[point] = dUx[Point];
	dUy[point] = dUy[Point];
      }
}

