// ======================================================
// authors: Sebastien Bougleux
// (C) 04 march 2010
// modifications:
// - 
// ======================================================

#include "anisoVoronoi2Diterative.h"
#include <cstring>
//================================================================
int labels = 0;
int* v4n = NULL;

//================================================================
struct BEdge // defines a boundary edge
//================================================================
{
  int label1;
  int label2;
  int pmax;
  int seed1;
  int seed2;
  bool mark;
  BEdge() 
  : label1(-1), label2(-1), pmax(-1), seed1(-1), seed2(-1), mark(false) { }
  BEdge(const int &l, const int &s)
  : label1(l), label2(l), pmax(-1), seed1(s), seed2(-1), mark(false) { }
  BEdge(const int &l1, const int &s1, const int &l2, const int &s2)
  : label1(l1), label2(l2), pmax(-1), seed1(s1), seed2(s2), mark(false) { }
  ~BEdge() {}
};

//================================================================
std::map<int,BEdge*> *BMap = NULL;
std::vector<BEdge*> *BEMap = NULL;

/*//================================================================
int findBoundaryPoint()
// returns the first bounday point encountered
//================================================================
{
  int point, npoint, s = size-Nx-1, k;
  short q;
  for (point = Nx+1; point < s; point++)
    {
      q = Q[point];
      if (q == kBorder || q == kSeed) continue;
      for (k = 0; k < connectivity_large; k++)
	{
	  npoint = point + NeighborhoodLarge[k];
	  if (Q[npoint] == kBorder) // found boundary point
	    {
	      return point;
	    }
	}
    }
  return -1;
}
*/

//================================================================
void GaussSiedelIterateBoundary()
//================================================================
{
  int point, npoint, k, v, Point;
  bool is_updated = false;
  short q;
  //------------------------------------------------------------
  // iterative propagation
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
	      //dUx[point]   = U_n_D[1];
	      //dUy[point]   = U_n_D[2];
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
	}
    }
}

//================================================================
void TsitsiklisTwoPoints(double* res, double* M, const double &k, 
			 const double &u, double* z1, double* z2,
			 int V1, int V2, int &Vr)
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
      Vr = V2;
      return;
    }
  if (k <= -sr11)
    {
      res[0] = 1.0;
      res[1] = k + u + sqrt(r11 + r22 + 2.0*r12);
      Vr = V1;
      return;
    }
  double sR = sqrt(R/(r11-k*k));
  if (r12 >= -k*sR)
    {
      res[0] = 0.0;
      res[1] = u + sqrt(r22);
      Vr = V2;
      return;
    }
  if (r12 <= (-r11-k*sR))
    {
      res[0] = 1.0;
      res[1] = k + u + sqrt(r11 + r22 + 2.0*r12);
      Vr = V1;
      return;
    }
  res[0] = -(r12 + k*sR) / r11;
  res[1] = res[0]*k + u + sR;
  // modif : compute new Voronoi label from V1 and V2
  if (V1 == V2)  // case same real label
    {
      Vr = V1;  // or V2
      return;
    }
  Vr = (res[0] > 0.5 ? V1 : V2);
}

//================================================================
void TsitsiklisQuadrant(double* res, int point, int simplex_idx, 
			const double &Ur, double *UComp, int *VComp,
			int &Vr)
//================================================================
{
  int npoint1, npoint2, i, simpIdx = 2*simplex_idx, nPoint1, nPoint2;
  double k, u, Norm;
  double* M;
  double U1, U2;
  int V1, V2;
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
  if (b_point1)  { U1 = UComp[nPoint1]; V1 = VComp[nPoint1]; }
  else           { U1 = INFINITE;  V1 = -1; }
  if (b_point2)  { U2 = UComp[nPoint2]; V2 = VComp[nPoint2]; }
  else           { U2 = INFINITE; V2 = -1; }
  //----------------------------------------------------------
  if (b_point1 && b_point2)
    {
      y1 = X12 + simpIdx;
      y2 = X2 + simpIdx;
      k = U1-U2; u = U2;
      TsitsiklisTwoPoints(res, M, k, u, y1, y2, V1, V2, Vr);
      //TsitsiklisNorm(&Norm, M, res[0], y1, y2);
      //for(i = 2; i < 4; i++) res[i] = (res[0]*y1[i-2] + y2[i-2])/Norm;
      return;
    }
  //----------------------------------------------------------
  if (b_point1 && !b_point2 )
    {
      y1 = X1 + simpIdx;
      y2 = X1 + simpIdx;
      u = U1;
      TsitsiklisOnePoint(res, M, u, y1);
      res[1] = res[0];
      res[0] = 1.0;
      //TsitsiklisNorm(&Norm, M, 0, y1, y2);
      //for(i = 2; i < 4; i++) res[i] = y2[i-2]/Norm;
      // compute Voronoi index
      Vr = V1;
      return;
    }
  //----------------------------------------------------------
  if (!b_point1 && b_point2 )
    {
      y1 = X2 + simpIdx;
      y2 = X2 + simpIdx;
      u = U2;
      TsitsiklisOnePoint(res, M, u, y1);
      res[1] = res[0];
      res[0] = 0;
      TsitsiklisNorm(&Norm, M, 0, y1, y2);
      for(i = 2; i < 4; i++) res[i] = y2[i-2]/Norm;
      // compute Voronoi index
      Vr = V2;
    }
  // error !b_point1 && !b_point2 (cannot arrive)
}

//================================================================
bool TsitsiklisUpdate(const int &point, double *UComp, int *VComp, // input
		      double &Ur, int &Vr) // input/output 
/*
COMMENTS : Ur must be initialized with UComp[Point]
*/
//================================================================
{
  int i, Vt;
  bool is_updated = false;
  //--------------------------------------------------------------
  for (i = 0; i < NB_simplicies; i++)
    {
      //----------------------------------------------------------
      TsitsiklisQuadrant(TAB, point, i, Ur, UComp, VComp, Vt);
      //----------------------------------------------------------
      if (TAB[1] < Ur)
	{
	  Ur = TAB[1];
	  Vr = Vt;
	  is_updated = true;
	}
    }
  return is_updated;
}

//================================================================
void GaussSiedelIterateInterior(double *UComp, int *VComp)
//================================================================
{
  int point, npoint, k, Point, Vr;
  bool is_updated = false;
  short q;
  double Ur, Uinit;
  while (!waiting_Points.empty()) 
    {
      point = waiting_Points.front();
      waiting_Points.pop();
      Q[point] = kUnqueuedEstimated;
      Point = comp_to_real(point);
      Uinit = Ur = UComp[Point];
      is_updated = TsitsiklisUpdate(point, UComp, VComp, Ur, Vr);
      if (is_updated && fabs(Ur - Uinit) > tol)
	{                
	  UComp[Point] = Ur;
	  VComp[Point] = Vr;
	  for (k = 0; k < connectivity_large; k++)
	    {
	      npoint = point + NeighborhoodLarge[k];
	      q = Q[npoint];
	      if (q == kUnqueuedEstimated)
		{
		  waiting_Points.push(npoint);
		  Q[npoint] = kEnqueuedEstimated;
		}
	    } 
	}
    }
}

/*
//================================================================
void TsitsiklisQuadrantInterior(double* res, int point, int simplex_idx)
//================================================================
{
  int npoint1, npoint2, i, simpIdx = 2*simplex_idx, nPoint1, nPoint2;
  double k, u, Norm;
  double* M;
  double U1, U2;
  double V1, V2; // modif
  double Ur = UComp[comp_to_real(point)];
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
  if (b_point1)  { U1 = UComp[nPoint1]; V1 = VComp[nPoint1]; }
  else           { U1 = INFINITE;  V1 = -1; }
  if (b_point2)  { U2 = UComp[nPoint2]; V2 = VComp[nPoint2]; }
  else           { U2 = INFINITE; V2 = -1; }
  //----------------------------------------------------------
  if (b_point1 && b_point2)
    {
      y1 = X12 + simpIdx;
      y2 = X2 + simpIdx;
      k = U1-U2; u = U2;
      TsitsiklisTwoPoints(res, M, k, u, y1, y2, V1, V2);
      TsitsiklisNorm(&Norm, M, res[0], y1, y2);
      for(i = 2; i < 4; i++) res[i] = (res[0]*y1[i-2] + y2[i-2])/Norm;
      return;
    }
  //----------------------------------------------------------
  if (b_point1 && !b_point2 )
    {
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
  }*/

//================================================================
 //bool TsitsiklisUpdateInterior(int point, double* res)
/*
COMMENTS : 
*/
//================================================================
 /*{
  int i;
  bool is_updated = false;
  double Norm;
  double *mat, det;
  int Point = comp_to_real(point);
  double Ur = UComp[Point], dUxr, dUyr ;
  double Vres = -1;
  mat = Metric + 3*point;
  //--------------------------------------------------------------
  for (i = 0; i < NB_simplicies; i++)
    {
      //----------------------------------------------------------
      TsitsiklisQuadrantInterior(TAB, point, i);
      //----------------------------------------------------------
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
  return is_updated;
}


//================================================================
void GaussSiedelIterateInteriorSeed()
//================================================================
{
  int point, npoint, k, v, Point;
  bool is_updated = false;
  short q;
  //------------------------------------------------------------
  // iterative propagation
  while (!waiting_Points.empty()) 
    {
      point = waiting_Points.front();
      waiting_Points.pop();
      is_updated = TsitsiklisUpdateInterior(point, U_n_D);
      Q[point] = kUnqueuedEstimated;
      if (is_updated)
	{
	  Point = comp_to_real(point);
	  if (fabs(U_n_D[0] - UComp[Point]) > tol)
	    {                
	      UComp[Point] = U_n_D[0];
	      VComp[Point] = (int)U_n_D[3];
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
	}
    }
}

//================================================================
bool addInteriorSeed(const int &point, const int &label, const int &Point)
//================================================================
{
  //------------------------------------------------------------
  // initialize
  VComp[Point] = label;
  UComp[Point] = 0.0;
  Q[point] = kSeed;
  for (int k = 0; k < connectivity_large; k++)
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
  //------------------------------------------------------------
  // compute new distance and Voronoi maps
  GaussSiedelIterateInteriorSeed();
  checkBoundary();
  return res;
}*/

//================================================================
void addSeed(const int &point, const int &label)
//
//================================================================
{
  // add it as a seed
  int Point, k, npoint;
  short q;
  Point = comp_to_real(point);
  //--------------------------------------------------------
  // update maps
  Q[point] = kSeed;
  UInit[Point] = 0.0;
  VInit[Point] = label;
  Seeds[Point] = label;
  //std::cerr << "add seed label " << label << std::endl;
  //--------------------------------------------------------
  // initialize neighbors
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

//================================================================
int findBoundaryMaxPointSecond(std::queue<BEdge*> &q)
//================================================================
{
  std::map<int,BEdge*>::iterator beit = BMap->begin(), bend = BMap->end();
  int point, Point, nb_edges = BEMap->size(), v;
  for (; beit != bend; beit++)
    {
      point = beit->first;
      if (Q[point] == kSeed) continue;
      Point = comp_to_real(point);
      v = VInit[Point];
      if (beit->second->mark == false && 
	  v != beit->second->label1 && v != beit->second->label2)
	{
	  beit->second->mark = true;
	  q.push(beit->second);
	  nb_edges--;
	  if (nb_edges == 0) break;
	}
    }
  return q.size();
}

//================================================================
BEdge* isBounbdaryPoint(const int &p)
//================================================================
{
  std::map<int,BEdge*>::iterator beit = BMap->find(p);
  return (beit != BMap->end() ? beit->second : NULL);
}

//================================================================
bool areBoundaryNeighbors(const int &p, const int &np, const int &k)
//================================================================
{
  switch (k)
    {
    case 0: case 2:
      if (Q[np - Nx] == kBorder) return true;
      if (Q[p - Nx] == kBorder) return true;
      if (Q[np + Nx] == kBorder) return true;
      if (Q[p + Nx] == kBorder) return true;
      break;
    case 1: case 3:
      if (Q[np - 1] == kBorder) return true;
      if (Q[p - 1] == kBorder) return true;
      if (Q[np + 1] == kBorder) return true;
      if (Q[p + 1] == kBorder) return true;
    }
  return false;
}

//================================================================
void findAdjacentBoundaryPoints(const int &p, std::vector<int> &vp)
//================================================================
{
  int np = p+1;
  if (Q[np] != kBorder && Q[np] != kSeed)
    {
      if (areBoundaryNeighbors(p, np, 0)) vp.push_back(np);
    }
  np = p-Nx;
  if (Q[np] != kBorder && Q[np] != kSeed)
    {
      if (areBoundaryNeighbors(p, np, 1)) vp.push_back(np);
    }
  np = p-1;
  if (Q[np] != kBorder && Q[np] != kSeed)
    {
      if (areBoundaryNeighbors(p, np, 2)) vp.push_back(np);
    }
  np = p+Nx;
  if (Q[np] != kBorder && Q[np] != kSeed)
    {
      if (areBoundaryNeighbors(p, np, 3)) vp.push_back(np);
    }
}

//================================================================
bool constructBoundary(const int &p, BEdge* e)
//================================================================
{
  // find two adjacent boundary points of p
  std::vector<int> vp; // adjacent boundary points
  findAdjacentBoundaryPoints(p, vp);
  if (vp.size() == 0) return false;
  if (vp.size() > 2) return false;
  // find boundary from one of these adjacent points
  int point, npoint, k;
  short q;
  point = vp[0];
  Q[point] = kEnqueued;
  BMap->insert(std::make_pair(point, e));
  while (point > 0)
    {
      for (k = 0; k < 4; k++)
	{
	  npoint = point + v4n[k];
	  q = Q[npoint];
	  if (q == kBorder || q == kEnqueued) continue;
	  if (q == kSeed)
	    {
	      point = -1;
	      break;
	    }
	  if (areBoundaryNeighbors(point, npoint, k))
	    {
	      point = npoint;
	      Q[point] = kEnqueued;
	      BMap->insert(std::make_pair(point, e));
	      break;
	    }
	}
    }
  return true;
}

//================================================================
double onePointBoundary(const int &pcomp, const int &pnei,
			const double &u)
//================================================================
{
  double *M = Metric + 3*pcomp;
  int diff = std::abs(pcomp-pnei);
  if (diff == 1) return u+std::sqrt(M[0]);
  if (diff == Nx) return u+std::sqrt(M[1]);
  std::cerr << "error onePointBoundary" << std::endl;
  return -1;
}

//================================================================
class BoundaryUCompare
//================================================================
{
 public:
  BoundaryUCompare() { }
  bool operator() (const int &p1, const int &p2) const
  {
    return (Utmp[comp_to_real(p1)] > Utmp[comp_to_real(p2)]);
  }
};

typedef std::priority_queue<int, std::vector<int>, BoundaryUCompare> BoundaryUQueue;

//================================================================
bool constructBoundaryFM(const int &p, BEdge* e)
// TO CHANGE TO TAKE INTO ACCOUNT ANISOTROPY: ok to check
//================================================================
{
  // find two adjacent boundary points of p
  std::vector<int> vp; // adjacent boundary points
  findAdjacentBoundaryPoints(p, vp);
  if (vp.size() == 0) return false;
  if (vp.size() > 2) return false;
  // find boundary from one of these adjacent points
  int point, npoint, k;
  short q;
  
  // insert them to the min heap and propagate
  //std::queue<int> PQ;
  BoundaryUQueue PQ;
  Q[vp[0]] = kEnqueued;
  Q[vp[1]] = kEnqueued;
  BMap->insert(std::make_pair(vp[0], e));
  BMap->insert(std::make_pair(vp[1], e));
  // update distance
  Utmp[comp_to_real(vp[0])] = onePointBoundary(vp[0], p, 0);
  Utmp[comp_to_real(vp[1])] = onePointBoundary(vp[1], p, 0);
  PQ.push(vp[0]);
  PQ.push(vp[1]);
  while (!PQ.empty())
    {
      point = PQ.top();
      PQ.pop();
      // find next and update distance
      for (k = 0; k < 4; k++)
	{
	  npoint = point + v4n[k];
	  q = Q[npoint];
	  if (q == kBorder || q == kEnqueued || q == kSeed) continue;
	  if (areBoundaryNeighbors(point, npoint, k))
	    {
	      Utmp[comp_to_real(npoint)] = onePointBoundary(npoint,
							    point,
							    Utmp[comp_to_real(point)]);
	      point = npoint;
	      Q[point] = kEnqueued;
	      BMap->insert(std::make_pair(point, e));
	      PQ.push(point);
	      break;
	    }
	}
    }
  // point is the farthest on the boundary from p
  e->pmax = point;
  return true;
}

//================================================================*
bool propagateToSeed(const int &p)
//================================================================*
{
  std::queue<int> PQ;
  PQ.push(p);
  Q[p] = kEnqueued;
  int point, npoint, k;
  short q;
  while (!PQ.empty())
    {
      point = PQ.front();
      PQ.pop();
      // find next and update boundary edge map
      for (k = 0; k < 4; k++)
	{
	  npoint = point + v4n[k];
	  q = Q[npoint];
	  if (q == kBorder || q == kEnqueued) continue;
	  if (q == kSeed)
	    {
	      return true;
	    }
	  if (areBoundaryNeighbors(point, npoint, k))
	    {
	      (*BMap)[npoint] = (*BMap)[point];
	      PQ.push(npoint);
	      Q[npoint] = kEnqueued;
	      break;
	    }
	}
    }
  return false;
}

//================================================================*
bool propagateToSeed(const int &p, std::queue<int> &sauv, int &res)
//================================================================*
{
  std::queue<int> PQ;
  PQ.push(p);
  Q[p] = kEnqueued;
  sauv.push(p);
  int point, npoint, k;
  short q;
  while (!PQ.empty())
    {
      point = PQ.front();
      PQ.pop();
      // find next and update boundary edge map
      for (k = 0; k < 4; k++)
	{
	  npoint = point + v4n[k];
	  q = Q[npoint];
	  if (q == kBorder || q == kEnqueued) continue;
	  if (q == kSeed)
	    {
	      res = npoint;
	      return true;
	    }
	  if (areBoundaryNeighbors(point, npoint, k))
	    {
	      (*BMap)[npoint] = (*BMap)[point];
	      PQ.push(npoint);
	      Q[npoint] = kEnqueued;
	      sauv.push(npoint);
	      break;
	    }
	}
    }
  return false;
}

//================================================================*
bool constructInitialBoundary()
//================================================================
{ 
  kEnqueued = 10;
  // creates a map of boundary points and an initial boundary edge
  BMap = new std::map<int,BEdge*>; // to link points and boundary edges
  BEMap = new std::vector<BEdge*>; // to save boundary edges
  
  int point, s = size-Nx-1, k, npoint, x, y, Point, i;
  short q;
  std::queue<int> PQ_seeds;
      
  // detect boundary curves and add a 1st seed on each curve
  for (x = 0; x < nx; x++)
    {
      for (y = 0; y < ny; y++)
	{
	  Point = x + y*nx;           // original matrix
	  point = (x+1) + (y+1)*Nx;   // extended matrix
	  q = Q[point];
	  if (q == kBorder || q == kSeed || q == kEnqueued) continue;
	  for (k = 0; k < connectivity_large; k++)
	    {
	      npoint = point + NeighborhoodLarge[k];
	      if (Q[npoint] == kBorder) // found boundary point
		{
		  BEdge *be = new BEdge(++labels, point);
		  BEMap->push_back(be);
		  Q[point] = kSeed;
		  Utmp[comp_to_real(point)] = 0.0;
		  PQ_seeds.push(point);
		  Seeds[comp_to_real(point)] = labels;
		  constructBoundaryFM(point, be);
		  break;
		}
	    }
	}
    }

  // Add a second seed on each boundary curve
  // as the max from 1st seed
  // => two sub-curves on each curve
  int eend = BEMap->size();
  std::vector<int> vp;
  double u;
  for (i = 0; i < eend; i++)
    {
      point = (*BEMap)[i]->pmax;  // get the max point from seed
      Utmp[comp_to_real(point)] = 0.0;
      Q[point] = kSeed;
      PQ_seeds.push(point);
      Seeds[comp_to_real(point)] = ++labels;
      // update boundary edge and create a new one
      BEdge *be = new BEdge(labels, point, 
			    (*BEMap)[i]->label1, (*BEMap)[i]->seed1);
      (*BEMap)[i]->seed2 = point;
      (*BEMap)[i]->label2 = labels;
      BEMap->push_back(be);
      // find new max points for each subcurve
      findAdjacentBoundaryPoints(point, vp);
      if (vp.size() != 2)
	{
	  std::cerr << "error: constructInitialBoundary(), adjacent points : " 
		    << vp.size() << std::endl;
	  return false;
	}
      // insert them to the min heap and propagate
      //std::queue<int> PQ;
      BoundaryUQueue PQ;
      kEnqueued++;
      Q[vp[0]] = kEnqueued;
      Q[vp[1]] = kEnqueued;
      (*BMap)[vp[0]] = be; // change boundary edge
      // update distance and propagate
      Utmp[comp_to_real(vp[0])] = onePointBoundary(vp[0], point, 0);
      Utmp[comp_to_real(vp[1])] = onePointBoundary(vp[1], point, 0);
      PQ.push(vp[0]);
      PQ.push(vp[1]);
      vp.clear();
      while (!PQ.empty())
	{
	  point = PQ.top();
	  PQ.pop();
	  // find next and update distance
	  for (k = 0; k < 4; k++)
	    {
	      npoint = point + v4n[k];
	      q = Q[npoint];
	      if (q == kBorder || q == kEnqueued || q == kSeed) continue;
	      if (areBoundaryNeighbors(point, npoint, k))
		{
		  u = onePointBoundary(npoint, point, Utmp[comp_to_real(point)]);
		  if (u  < Utmp[comp_to_real(npoint)])
		    Utmp[comp_to_real(npoint)] = u;
		  else 
		    {
		      (*BMap)[point]->pmax = npoint;
		      (*BMap)[npoint] = (*BMap)[point];
		      propagateToSeed(npoint);
		      k=4;
		      break;
		    }
		  Q[npoint] = kEnqueued;
		  (*BMap)[npoint] = (*BMap)[point];
		  PQ.push(npoint);
		  break;
		}
	    }
	}
    }

  // re-init Q map from boudary data
  kEnqueued = -1;
  std::map<int,BEdge*>::iterator beit = BMap->begin(), bend = BMap->end();
  for (; beit != bend; beit++)
    {
      point  = beit->first;
      Q[point] = kUnqueued;
    }
  // init gauss-siedel
  while (!PQ_seeds.empty())
    {
      point = PQ_seeds.front();
      PQ_seeds.pop();
      addSeed(point, Seeds[comp_to_real(point)]);
    }
  /*std::vector<BEdge*>::iterator eit = BEMap->begin(), beend = BEMap->end();
  for (; eit != beend; eit++)
    {
      point = (*eit)->seed1;
      std::cerr << (*eit)->label1 << " " << (*eit)->label2 << std::endl;
    }
  */
  return true;
}

//================================================================
void addBoundarySeed(BEdge *e)
//================================================================
{
  int point = e->pmax, psauv;
  psauv = point;
  Utmp[comp_to_real(point)] = 0.0;
  Q[point] = kSeed;
  Seeds[comp_to_real(point)] = ++labels;
  // update boundary edge and create a new one
  BEdge *be = new BEdge(labels, point, e->label1, e->seed1);
  e->seed2 = point;
  e->label2 = labels;
  BEMap->push_back(be);
  // find new max points for each subcurve
  std::vector<int> vp;
  findAdjacentBoundaryPoints(point, vp);
  if (vp.size() != 2)
    {
      std::cerr << "error: constructInitialBoundary(), adjacent points : " 
		<< vp.size() << std::endl;
      return;
    }
  // insert them to the min heap and propagate
  std::queue<int> PQ, Reinit;
  PQ.push(vp[0]);
  PQ.push(vp[1]);
      
  kEnqueued = 10;
  int k, npoint, res;
  short q;
  double u;
  Q[vp[0]] = kEnqueued;
  Q[vp[1]] = kEnqueued;
  Reinit.push(vp[0]);
  Reinit.push(vp[1]);
  (*BMap)[vp[0]] = be; // change boundary edge
  // update distance and propagate
  Utmp[comp_to_real(vp[0])] = 1;
  Utmp[comp_to_real(vp[1])] = 1;
  vp.clear();
  while (!PQ.empty())
    {
      point = PQ.front();
      PQ.pop();
      // find next and update distance
      for (k = 0; k < 4; k++)
	{
	  npoint = point + v4n[k];
	  q = Q[npoint];
	  if (q == kBorder || q == kEnqueued || q == kSeed) continue;
	  if (areBoundaryNeighbors(point, npoint, k))
	    {
	      u = Utmp[comp_to_real(point)] + 1.;
	      if (u  < Utmp[comp_to_real(npoint)])
		Utmp[comp_to_real(npoint)] = u;
	      else // finish the edge
		{
		  (*BMap)[point]->pmax = npoint;
		  (*BMap)[npoint] = (*BMap)[point];
		  if (propagateToSeed(npoint, Reinit, res))
		    {
		      (*BMap)[point]->label1 = labels;
		      (*BMap)[point]->seed1 = psauv;
		      (*BMap)[point]->label2 = Seeds[comp_to_real(res)];
		      (*BMap)[point]->seed2 = res;
		    }
		  else
		    {
		      std::cerr << "error : addBoundarySeed" << std::endl;
		    }
		  break;
		 
		}
	      Q[npoint] = kEnqueued;
	      (*BMap)[npoint] = (*BMap)[point];
	      PQ.push(npoint);
	      Reinit.push(npoint);
	      break;
	    }
	}
    }

  // re-init ...
  kEnqueued = -1;
  while (!Reinit.empty())
    {
      point = Reinit.front();
      Reinit.pop();
      Q[point] = kUnqueuedEstimated;
    }
  // add seed
  addSeed(psauv, labels);
}

//================================================================
void initializeFPM()
// get initial boundary seeds and edges with FPS on the domain bounary
//================================================================
{
  v4n = new int[4];
  v4n[0] = 1; v4n[1] = -Nx; v4n[2] = -1; v4n[3] = Nx;
  
  constructInitialBoundary();
  GaussSiedelIterateBoundary();
  BEdge *e;
  std::queue<BEdge*> q;
  findBoundaryMaxPointSecond(q);
  std::vector<BEdge*>::iterator eit, eend;
      
  while (findBoundaryMaxPointSecond(q))
    {
      while (!q.empty())
	{
	  addBoundarySeed(q.front());
	  q.pop();
	}
      GaussSiedelIterateBoundary();
      eit = BEMap->begin();
      eend = BEMap->end();
      for (; eit != eend; eit++) (*eit)->mark = false;
      
    }
}

//================================================================
int findMax(double umax = 0)
//================================================================
{
  int pmax = -1, point, Point;
  double utmp;
  short q;
  for (Point = 0; Point < nxny; Point++)
    {
      point = real_to_comp(Point);
      q = Q[point];
      if (q == kBorder || q == kSeed) continue;
      utmp = UInit[Point];
      if (utmp > umax)
	{
	  umax = utmp;
	  pmax = point;
	}
    }
  return pmax;
}

//================================================================
void addSeedInterior(const int &point, int label)
//
//================================================================
{
  // add it as a seed
  int Point = comp_to_real(point), k, npoint;
  //--------------------------------------------------------
  // update maps
  Q[point] = kSeed;
  UInit[Point] = 0.0;
  VInit[Point] = label;
  Seeds[Point] = label;
//  std::cerr << "add seed label " << label << std::endl;
  //--------------------------------------------------------
  // initialize neighbors
  for (k = 0; k < connectivity_large; k++)
    {
      npoint = point + NeighborhoodLarge[k];
      if (Q[npoint] == kUnqueuedEstimated)
	{
	  waiting_Points.push(npoint);
	  Q[npoint] = kEnqueuedEstimated;
	}
    }
}

//================================================================
template<class T>
void cpy(T *m1, T *m2)
//================================================================
{
  for (int i = 0; i < nxny; i++) m1[i] = m2[i];
}

//================================================================
void gridFPM(int nb_seeds)
//================================================================
{
  initializeFPM();
  // to compute interior seeds
  double *UComp = new double[nxny];
  int *VComp = new int[nxny];
  std::vector<BEdge*>::iterator eit, eend;
  int pmax;
  BEdge* e = NULL;
  while (nb_seeds > 0)
    {
      // find max
      pmax = findMax();
      e = isBounbdaryPoint(pmax);
      if (e)
	{
	  //std::cerr << "boundary seed" << std::endl;
	  addBoundarySeed(e);
	  GaussSiedelIterateBoundary();
	  nb_seeds--;
	  eend = BEMap->end();
	  for (eit = BEMap->begin(); eit != eend; eit++) (*eit)->mark = false;
	  // save data
	  //std::memcpy(UComp, UInit, sizeof(double)*nxny);
	  //std::memcpy(VComp, VInit, sizeof(int)*nxny);
	  continue;
	}
      // save data
      cpy(UComp, UInit);
      cpy(VComp, VInit);
      // insert max as a new seed
      addSeedInterior(pmax, labels+1);
      GaussSiedelIterateVoronoi();//GaussSiedelIterateInterior(UComp, VComp);
      // check boundary
      std::queue<BEdge*> q;
      if (findBoundaryMaxPointSecond(q)) // border intersect new Voronoi cell
	{
	  //std::memcpy(UInit, UComp, sizeof(double)*nxny);
	  //std::memcpy(VInit, VComp, sizeof(int)*nxny);
	  cpy(VInit, VComp);
	  cpy(UInit, UComp);
	  Q[pmax] = kUnqueuedEstimated;
	  Seeds[comp_to_real(pmax)] = 0;
	  //std::cerr << "detected" << std::endl;
	  nb_seeds -= q.size();
	  do
	    {
	      while (!q.empty())
		{
		  addBoundarySeed(q.front());
		  q.pop();
		}
	      GaussSiedelIterateBoundary();
	      eend = BEMap->end();
	      for (eit = BEMap->begin(); eit != eend; eit++) (*eit)->mark = false;
	      //break;
	    } while (findBoundaryMaxPointSecond(q));
	}
      else
	{
	  labels++;
	  nb_seeds--;
	}
    }
  delete[] UComp;
  delete[] VComp;
}

/*//================================================================
struct Triangle
//================================================================
{
  int s1;
  int s2;
  int s3;
  int voro;
};

//================================================================
void extractVoronoi()
//
//================================================================
{
  int Point, point;
  for (point = 0; point < nxny; point++)
    {
      Point = real_to_comp(point);
      if (Q[Point] == kSeed)
	{
	  
	  
	}
      //--------------------------------------------------------
      for (k = 0; k < connectivity_large; k++)
	{
	  npoint = point + NeighborhoodLarge[k];
	  q = Q[npoint];
	  if (q == kBorder)  // boundary point
	    {
	      // add it to the vector
	      break;
	    }
	}
      UInit[point] = ;
      VInit[point] = ;
    }
  
}
*/

//================================================================
void DeleteData()
//================================================================
{
  if (BMap)
    {
      delete BMap;
    }
  if (BEMap)
    {
      for (int i = 0; i < BEMap->size(); i++) delete (*BEMap)[i];
      delete BEMap;
    }
  if (v4n) delete[] v4n;
  labels = 0;
}