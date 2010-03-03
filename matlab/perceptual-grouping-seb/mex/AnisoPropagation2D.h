//================================================================
//================================================================
// File: AnisoPropagation2D.h
// (C) 12/2008 by F. Benmansour and S. Bougleux
//================================================================
//================================================================
#ifndef __ANISOPROPAGATION2D_H__
#define __ANISOPROPAGATION2D_H__
#include<limits>
#include<iostream>
#include "minheap.h"

#define NPOINT(p,k) ((p)+(v8k[(k)]))
#define SIMPLEX(p,k,l) ((p)+(_simplex[(k)][(l)]))
//================================================================
/**
 * @class AnisoPropagation2D
 * @brief Main class of the library
 * @param T type of action values (must be a real type, float by default)
 * @param Integer type of indices, labels (can be a positive integer type, int by default)
 */
//================================================================
template<class T = float, class Integer = int>
class AnisoPropagation2D
//================================================================
{
  //================================================================
 protected:
  //--------------------------------------------------------------
  // given sizes and arrays
  short connectivity;     // point neighborhood
  Integer size;           // number of grid points
  Integer nx, ny, nxny;   // real grid size
  T hx, hy;                     // grid spacing in x and y
  //--------------------------------------------------------------
  T *_U;                          // anisotropic distance map
  Integer *_V;                    // Voronoi map (label map)
  Integer _labels;
  //T *L;                         // Euclidean distance map
  //T *dUx, *dUy;                 // gradient in x and y
  //--------------------------------------------------------------
  // class properties (sizes and arrays)
  Integer Nx, Ny;         // computing grid size
  Integer kDead, kOpen, kFar, kBorder, kSeed;  // state values
  //--------------------------------------------------------------
  T hx2, hy2, hxhy, sqrt_hx2_plus_hy2;
  Integer *S;              // state map
  T *q_gradient;
  //--------------------------------------------------------------
  Integer *v8k;           // 8-adjacency neighborhood
  Integer *simplicies1;
  Integer *simplicies2;
  Integer **_simplex;
  short **_simplex_sign;
  T *signsX, *signsY;
  T *h1;
  T *h2;
  T *h1_h2;
  T *h22;
  T *M1, *M2, *M3;                // metric tensor
  //--------------------------------------------------------------
  typedef MinHeap<Integer, T>  Tree;
  Tree *tree;                     // min heap
  //--------------------------------------------------------------
  // protected functions
  //--------------------------------------------------------------
  /**
   * @brief 
   */
  void InitializeNeighborhoods();
  //--------------------------------------------------------------
  void InitializeArrays();
  void SetTensorField22(T*);
  //--------------------------------------------------------------
  Integer addSeed(Integer, Integer);
  //--------------------------------------------------------------
  /**
   * @brief 
   */
  void clearMaps();
  //--------------------------------------------------------------
  /**
   * @brief Set a set of points to current data
   */
  Integer setSeeds(double*, Integer);
  //--------------------------------------------------------------
  /**
   * @brief Set S with kBorder for each point of the virtual boundary
   */
  void initializeGridBoundary();
  //--------------------------------------------------------------
  //T TsitsiklisQuadrantLength(T&, const T&, const T&, const T&, int);
 //Integer tsitsiklisQuadrantGradient(const T&, const T&, const T&,
 //                                  const T&, const T&,
 //                                  const Integer&, const Integer&, 
 //                                  short, T*);
  //--------------------------------------------------------------
  //virtual void tsitsiklisUpdateNeighbors(const Integer&);
 //bool tsitsiklisUpdate(const Integer&);
 // bool tsitsiklisUpdateGradient(Integer);
 // void propagate(T);
  T solver1D(const Integer&, const Integer&, const short&);
  T solver2D(const Integer&, const Integer&, const Integer&, const short&, Integer&);
  T update(const Integer&, const Integer&, const short&, Integer&);
  //--------------------------------------------------------------
  void pos_to_point(const Integer&, Integer&, Integer&);
  Integer pos_to_point_y(const Integer&);
  short orient2D(Integer, Integer, Integer, Integer, Integer, Integer);
  //--------------------------------------------------------------
  //================================================================
 public:
  T INFINITE;
  Integer INFINITE_INT;
  //--------------------------------------------------------------
  /**
   * @brief Constructor
   * @param tf is a 2D tensor field as a dc*dl*4 matrix
   * @param dc is the number of columns
   * @param dl is the number of lines
   * @param U is the output action map (distance map)
   * @param V is the optional output Voronoi map (label map)
   * @param sx is a scaling factor on lines
   * @param sy is a scaling factor on columns
   */
  AnisoPropagation2D(T* tf, int dx, int dy, T* U, Integer* V, /*T*, T*, T*, */ T sx, T sy);
  //--------------------------------------------------------------
  /**
   * @brief Destructor
   */
  ~AnisoPropagation2D();
  //--------------------------------------------------------------
  /**
   * @brief Add a point to 
   */
  Integer addSeed(int, int, Integer);
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  bool updateNeighbors(const Integer&);
  template<class C> void propagation(C *cl, bool (C::*fct)(const Integer&))
  {
    Integer point;
    while (tree->Size() > 0)
      {
	point = tree->Pop();
	S[point] = kDead;
	(cl->*fct)(point);
      }
  }
  /**
   * @brief Global propagation from the current set of points 
   */
  bool propagation(double*, Integer, T);
  //--------------------------------------------------------------
  Integer GridSizeX();
  Integer GridSizeY();
  //--------------------------------------------------------------
  void resize();
  //--------------------------------------------------------------
  short orient2D(const Integer &p1, const Integer &p2, const Integer &p3);
};
//================================================================
//================================================================
//================================================================
//================================================================
// Functions
//================================================================
//================================================================
template<class T, class Integer>
AnisoPropagation2D<T,Integer>::AnisoPropagation2D(T* TF22, int _nx, int _ny, 
						  T *U, Integer *V = NULL,
						  //T *_L = NULL, T *_dUx = NULL, 
						  //T *_dUy = NULL,
						  T _hx = 1, T _hy = 1)
  : connectivity(8), size((_nx+2)*(_ny+2)), nx(_nx), ny(_ny), nxny(_nx*_ny), hx(_hx), hy(_hy),
     _U(U), _V(V), _labels(0), /*L(_L), dUx(_dUx), dUy(_dUy), */
     Nx(_nx+2), Ny(_ny+2), kDead(3), kOpen(2), kFar(1), kBorder(0), kSeed(1), 
     hx2(_hx*_hx), hy2(_hy*_hy), sqrt_hx2_plus_hy2(std::sqrt(_hx*_hx+_hy*_hy)), 
     S(new Integer[(_nx+2)*(_ny+2)]), q_gradient(new T[3]), 
     v8k(new Integer[9]),
     simplicies1(new Integer[8]), simplicies2(new Integer[8]), _simplex(new Integer*[8]), 
     _simplex_sign(new short*[8]), 
     signsX(new T[8]), signsY(new T[8]), h1(new T[8]), h2(new T[8]),
     h1_h2(new T[8]), h22(new T[8]), M1(NULL), M2(NULL), M3(NULL),
     tree(NULL), INFINITE(std::numeric_limits<T>::max()),
     INFINITE_INT(std::numeric_limits<Integer>::max())
{
  InitializeNeighborhoods();   // initialize neighborhoods
  InitializeArrays();          // initialize Voronoi, Euclidean length and gradient
  SetTensorField22(TF22);     // initialize tensor field
  //tree = new Tree(size, _U);    // initialize min heap
}

//================================================================
template<class T, class Integer>
AnisoPropagation2D<T,Integer>::~AnisoPropagation2D()
{
  if (S) delete[] S;
  if (v8k) delete[] v8k;
  if (simplicies1) delete[] simplicies1;
  if (simplicies2) delete[] simplicies2;
  if (_simplex)
    {
      for (short k = 0; k < connectivity; k++)
	delete[] _simplex[k];
      delete[] _simplex;
    }
  //if (signsX) delete[] signsX;
  //if (signsY) delete[] signsY;
  if (h1) delete[] h1;
  if (h2) delete[] h2;
  if (h1_h2) delete[] h1_h2;
  if (h22) delete[] h22;
  if (M1) delete[] M1;
  if (M2) delete[] M2;
  if (M3) delete[] M3;
  if (tree) delete tree;
  if (q_gradient) delete[] q_gradient;
}

//================================================================
template<class T, class Integer>
void AnisoPropagation2D<T,Integer>::InitializeNeighborhoods()
{
  simplicies2[7] = simplicies2[3] = v8k[0] = v8k[8] = -Nx-1;
  simplicies1[7] = simplicies1[6] = v8k[1] = -Nx; 
  simplicies2[6] = simplicies2[1] = v8k[2] = 1-Nx;
  simplicies1[1] = simplicies1[0] = v8k[3] = 1;
  simplicies2[4] = simplicies2[0] = v8k[4] = Nx+1;
  simplicies1[5] = simplicies1[4] = v8k[5] = Nx;
  simplicies2[5] = simplicies2[2] = v8k[6] = Nx-1;
  simplicies1[3] = simplicies1[2] = v8k[7] = -1;
  //--------------------------------------------------------------
  _simplex[0] = new Integer[2];
  _simplex_sign[0] = new short[2];
  for (short k = 1; k < connectivity; k++)
    {
      _simplex[k] = new Integer[2];
      _simplex_sign[k] = new short[2];
    }
  _simplex[0][0] = 1;
  _simplex[0][1] = Nx;
  _simplex[1][0] = Nx+1;
  _simplex[1][1] = Nx-1;
  _simplex[2][0] = Nx;
  _simplex[2][1] = -1;
  _simplex[3][0] = Nx-1;
  _simplex[3][1] = -Nx-1;
  _simplex[4][0] = -1;
  _simplex[4][1] = -Nx;
  _simplex[5][0] = -Nx-1;
  _simplex[5][1] = -Nx+1;
  _simplex[6][0] = -Nx;
  _simplex[6][1] = 1;
  _simplex[7][0] = -Nx+1;
  _simplex[7][1] = Nx+1;
  //--------------------------------------------------------------
  _simplex_sign[0][0] = 0;
  _simplex_sign[0][1] = 4;
  _simplex_sign[1][0] = 4;
  _simplex_sign[1][1] = 5;
  _simplex_sign[2][0] = 5;
  _simplex_sign[2][1] = 2;
  _simplex_sign[3][0] = 2;
  _simplex_sign[3][1] = 3;
  _simplex_sign[4][0] = 3;
  _simplex_sign[4][1] = 7;
  _simplex_sign[5][0] = 7;
  _simplex_sign[5][1] = 6;
  _simplex_sign[6][0] = 6;
  _simplex_sign[6][1] = 1;
  _simplex_sign[7][0] = 1;
  _simplex_sign[7][1] = 0;
  //--------------------------------------------------------------
  signsX[6] = signsY[5] = signsY[4] = signsX[4] = 
    signsY[2] = signsX[0] = signsX[1] = signsY[0] = -1.;
  signsY[7] = signsX[7] = signsY[6] = signsX[5] = 
    signsY[3] = signsX[3] = signsX[2] = signsY[1] = 1.;
  h2[7] = h2[6] = h2[5] = h2[4] = h1[3] = h1[2] = h1[1] = h1[0] = hy;
  h1[7] = h1[6] = h1[5] = h1[4] = h2[3] = h2[2] = h2[1] = h2[0] = hx;
  h1_h2[3] = h1_h2[2] = h1_h2[1] = h1_h2[0] = hy/hx;
  h22[3] = h22[2] = h22[1] = h22[0] = hx2;
  h1_h2[7] = h1_h2[6] = h1_h2[5] = h1_h2[4] = hx/hy;
  h22[7] = h22[6] = h22[5] = h22[4] = hy2;
}

//================================================================
template<class T, class Integer>
void AnisoPropagation2D<T,Integer>::InitializeArrays()
{
  if (_V == NULL) _V = new Integer[size];
  //if (L == NULL) L = new T[size];
  //if (dUx == NULL) dUx = new T[size];
  //if (dUy == NULL) dUy = new T[size];
}

//================================================================
template<class T, class Integer>
void AnisoPropagation2D<T,Integer>::SetTensorField22(T *TF)
{
  Integer x, y, point;
  if (M1 == NULL) M1 = new T[size];
  if (M2 == NULL) M2 = new T[size];
  if (M3 == NULL) M3 = new T[size];
  //--------------------------------------------------------------
  if ((hx != 1.0) || (hy != 1.0))
    {
      for (x = 0; x < nx; x++)
	for (y = 0; y < ny; y++)
	  {
	    point = (x+1) + (y+1)*Nx;
	    M1[point] = hx2* TF[x + y*nx];
	    M2[point] = hy2* TF[x + y*nx + 3*nxny];
	    M3[point] = hxhy*TF[x + y*nx +   nxny];
	  }
      return;
    }
  //--------------------------------------------------------------
  for (x = 0; x < nx; x++)
    for (y = 0; y < ny; y++)
      {
	point = (x+1) + (y+1)*Nx;
	M1[point] = TF[x + y*nx];
	M2[point] = TF[x + y*nx + 3*nxny];
	M3[point] = TF[x + y*nx + nxny];
      }
}

//================================================================
template<class T, class Integer>
void AnisoPropagation2D<T,Integer>::clearMaps()
{
  for (Integer x = 0; x < size; x++)
    {
      _U[x] = INFINITE;
      //L[x] = INFINITE;
      _V[x] = INFINITE_INT;
      S[x] = kFar;
    }
  if (tree) delete tree;      //->Clear(); ????
  tree = new Tree(size, _U);
  initializeGridBoundary();
}

//================================================================
template<class T, class Integer>
void AnisoPropagation2D<T,Integer>::initializeGridBoundary()
{
  Integer x, y, point;
  //--------------------------------------------------------------
  for (x = 0; x < Nx; x++)
    {
      y = 0;
      point = x + y*Nx;
      S[point] = kBorder;
      y = Ny-1;
      point = x + y*Nx;
      S[point] = kBorder;
    }
  //--------------------------------------------------------------
  for (y = 0; y < Ny; y++)
    {
      x = 0;
      point = x + y*Nx;
      S[point] = kBorder;
      x = Nx-1;
      point = x + y*Nx;
      S[point] = kBorder;
    }
}

//================================================================
/*template<class T, class Integer>
T AnisoPropagation2D<T,Integer>::TsitsiklisQuadrantLength(T &result, const T &Pc,
					      const T &Ua, const T &Ub, 
					      int QuadNb)
{
  T ha, hb, hb2, ha_hb;
  ha = h1[QuadNb];
  hb = h2[QuadNb];
  hb2= h22[QuadNb];
  ha_hb = h1_h2[QuadNb];
  //-----------------------------------------------------------
  if (Ua > Ub)
    {
      if ((Pc*hb) > ((Ua-Ub)*sqrt_hx2_plus_hy2/hb))
	return (Ua + ha_hb*sqrt(Pc*Pc*hb2 - (Ua-Ub)*(Ua-Ub)));
      return (Ub + Pc*sqrt_hx2_plus_hy2);
    }
  //else (Ua <= Ub)
  return (Ua + Pc*ha);
}
*/

//================================================================
/*
template<class T, class Integer>
Integer AnisoPropagation2D<T,Integer>::tsitsiklisQuadrantGradient(const T &m1, const T &m2, const T &m3,
								  const T &Ua, const T &Ub, 
								  const Integer &Va, const Integer &Vb, 
								  short QuadNb, T *result_gradient)
{
  T ha, hb, hb2, ha_hb;
  T k1 = Ua - Ub, k2 = Ub, alpha, r11, r12, r22, R;
  //-----------------------------------------------------------
  if(QuadNb < 4)
    {
      r11 = m2;
      r12 = m2 + m3*(signsX[QuadNb]*signsY[QuadNb]);
    }
  else
    {
      r11 = m1;
      r12 = m1 + m3*(signsX[QuadNb]*signsY[QuadNb]);
    }
  //-----------------------------------------------------------
  r22 = m2 + 2.0*m3*(signsX[QuadNb]*signsY[QuadNb]) + m1;
  R = r11*r22 - r12*r12; // always positive
  //-----------------------------------------------------------
  if (k1 < sqrt(r11))
    {
      if(k1 <= -sqrt(r11))
	{
	  if(QuadNb < 4)
	    {
	      result_gradient[0] = Ua + sqrt(m1);
	      result_gradient[1] = signsX[QuadNb]/sqrt(m1);
	      result_gradient[2] = 0.0;
	    }
	  else
	    {
	      result_gradient[0] = Ua + sqrt(m2);
	      result_gradient[1] = 0.0;
	      result_gradient[2] = signsY[QuadNb]/sqrt(m2);
	    }
	  return Va;
	}
      //-----------------------------------------------------------
      else
	{
	  //-------------------------------------------------------
	  if (r12 <= k1*sqrt(R/(r11-k1*k1)))
	    {
	      result_gradient[0] = Ub + sqrt(r22);
	      result_gradient[1] = signsX[QuadNb]/sqrt(r22);
	      result_gradient[2] = signsY[QuadNb]/sqrt(r22);
	      return Vb;
	    }
	  //-------------------------------------------------------
	  else if (r12 > (r11 + k1*sqrt(R/(r11-k1*k1))))
	    {
	      if (QuadNb < 4)
		{
		  result_gradient[0] = Ua + sqrt(m1);
		  result_gradient[1] = signsX[QuadNb]/sqrt(m1);
		  result_gradient[2] = 0.0;
		}
	      else
		{
		  result_gradient[0] = Ua + sqrt(m2);
		  result_gradient[1] = 0.0;
		  result_gradient[2] = signsY[QuadNb]/sqrt(m2);
		}
	      return Va;
	    }
	  //-------------------------------------------------------
	  else
	    {
	      alpha = (r12 - k1 * sqrt(R/(r11-k1*k1))) / r11;
	      result_gradient[0] = alpha*k1 + k2 + sqrt(R/(r11-k1*k1));
	      if (QuadNb < 4)
		{
		  result_gradient[1] = signsX[QuadNb] / sqrt(R/(r11-k1*k1));
		  result_gradient[2] = (1-alpha) * signsY[QuadNb] / sqrt(R/(r11-k1*k1));
		}
	      else
		{
		  result_gradient[1] = (1-alpha) * signsX[QuadNb] / sqrt(R/(r11-k1*k1));
		  result_gradient[2] = signsY[QuadNb] / sqrt(R/(r11-k1*k1));
		}
	      return ((Ua+sqrt(r11)) < (Ub+sqrt(r22)) ? Va : Vb);
	    }
	}
    }
  //-----------------------------------------------------------
  else // if (k1 >= sqrt(r11))
    {
      result_gradient[0] = Ub + sqrt(r22);
      result_gradient[1] = signsX[QuadNb]/sqrt(r22);
      result_gradient[2] = signsY[QuadNb]/sqrt(r22);
      return Vb;
    }
  return 0;
}
*/
//================================================================
template<class T, class Integer>
T AnisoPropagation2D<T,Integer>::solver2D(const Integer &point,
					  const Integer &npoint1, const Integer &npoint2,
					  const short &k, Integer &V)
{
  T Ua = _U[npoint1], Ub = _U[npoint2];
  T k1 = Ua - Ub, k2 = Ub, alpha, r11, r12, r22, R;
  T m1 = M1[point], m2 = M2[point], m3 = M3[point];
  //-----------------------------------------------------------
  if(k < 4)
    {
      r11 = m2;
      r12 = m2 + m3*(signsX[k]*signsY[k]);
    }
  else
    {
      r11 = m1;
      r12 = m1 + m3*(signsX[k]*signsY[k]);
    }
  //-----------------------------------------------------------
  r22 = m2 + 2.0*m3*(signsX[k]*signsY[k]) + m1;
  R = r11*r22 - r12*r12; // always positive
  //-----------------------------------------------------------
  if (k1 < sqrt(r11))
    {
      if(k1 <= -sqrt(r11))
	{
	  V = _V[npoint1];
	  if(k < 4) return Ua + sqrt(m1);
	  return Ua + sqrt(m2);
	}
      //-----------------------------------------------------------
      else
	{
	  //-------------------------------------------------------
	  if (r12 <= k1*sqrt(R/(r11-k1*k1)))
	    {
	      V = _V[npoint2];
	      return Ub + sqrt(r22);
	    }
	  //-------------------------------------------------------
	  else if (r12 > (r11 + k1*sqrt(R/(r11-k1*k1))))
	    {
	      V = _V[npoint1];
	      if (k < 4) return Ua + sqrt(m1);
	      else return Ua + sqrt(m2);
	    }
	  //-------------------------------------------------------
	  else
	    {
	      V = (Ua+sqrt(r11) < Ub+sqrt(r22) ? _V[npoint1] : _V[npoint2]);
	      return ((r12 - k1 * sqrt(R/(r11-k1*k1))) / r11)*k1 + k2 + sqrt(R/(r11-k1*k1));
	    }
	}
    }
  //-----------------------------------------------------------
  // if (k1 >= sqrt(r11))
  else
    {
      V = _V[npoint2];
      return Ub + sqrt(r22);
    }
  V = INFINITE_INT;
  return INFINITE;
}

//================================================================
template<class T, class Integer>
T AnisoPropagation2D<T,Integer>::solver1D(const Integer &point, const Integer &npoint, const short &k)
{
  switch (k)
    {
    case 0: case 4:
      return _U[npoint] + std::sqrt(M1[point]+M2[point]+2.*M3[point]);
    case 1: case 5:
      return _U[npoint] + std::sqrt(M2[point]);
    case 2: case 6:
      return _U[npoint] + std::sqrt(M1[point]+M2[point]-2.*M3[point]);
    case 3: case 7:
      return _U[npoint] + std::sqrt(M1[point]);
    }
  return -1.;
}

//================================================================
template<class T, class Integer>
T AnisoPropagation2D<T,Integer>::update(const Integer &point, const Integer &npoint, const short &k, 
					Integer &Vr)
{
  Integer npoint1, npoint2, Vtmp;
  T Ur = INFINITE, Utmp;
  bool is_updated = false;
  //--------------------------------------------------------------
  short i = _simplex_sign[k][0];
  npoint1 = point + simplicies1[i];
  npoint2 = point + simplicies2[i];
  if (S[npoint1] == kDead && S[npoint2] == kDead)
    if ((Utmp = solver2D(point, npoint1, npoint2, i, Vtmp))  < Ur)
      {
	Vr = Vtmp;
	Ur = Utmp;
      }
  //--------------------------------------------------------------
  i = _simplex_sign[k][1];
  npoint1 = point + simplicies1[i];
  npoint2 = point + simplicies2[i];
  if (S[npoint1] == kDead && S[npoint2] == kDead)
    if ((Utmp = solver2D(point, npoint1, npoint2, i, Vtmp)) < Ur)
      {
	Vr = Vtmp;
	Ur = Utmp;
      }
  //--------------------------------------------------------------
  if (Ur < INFINITE) return Ur;
  //--------------------------------------------------------------
  Vr = _V[npoint];
  return solver1D(point, npoint, k);
}

//================================================================
/*template<class T, class Integer>
bool AnisoPropagation2D<T,Integer>::tsitsiklisUpdate(const Integer &point)
{
  Integer npoint1, npoint2, Vr, Vtmp;
  T Ua, Ub, Ur = _U[point], Utmp;
  bool is_updated = false;
  //--------------------------------------------------------------
  // Get the U & L values for each neighbor
  for (short i = 0; i < connectivity; i++)
    {
      npoint1 = point + simplicies1[i];
      npoint2 = point + simplicies2[i];
      Ua = (S[npoint1] == kDead ? _U[npoint1] : INFINITE);
      Ub = (S[npoint2] == kDead ? _U[npoint2] : INFINITE);
      Vtmp = tsitsiklisQuadrant(M1[point], M2[point], M3[point],
				Ua, Ub, _V[npoint1], _V[npoint2], i, Utmp);
      if (Utmp < Ur)
        {
          Ur = Utmp;
          Vr = Vtmp;
          is_updated = true;
        }
    }
  //--------------------------------------------------------------
  if (is_updated) { _U[point] = Ur; _V[point] = Vr; }
  //--------------------------------------------------------------
  return is_updated;
}
*/
//================================================================
/*
template<class T, class Integer>
bool AnisoPropagation2D<T,Integer>::tsitsiklisUpdateGradient(Integer point)
{
  Integer npoint1, npoint2, Vr, Vtmp;
  T Ua, Ub, La, Lb, Ur = _U[point], Lr = L[point]; //, dUxr = dUx[point], dUyr = dUy[point];
  bool is_updated = false;
  //--------------------------------------------------------------
  // Get the U & L values for each neighbor
  for (short i = 0; i < connectivity; i++)
    {
      npoint1 = point + simplicies1[i];
      npoint2 = point + simplicies2[i];
      if (S[npoint1] == kDead) { Ua = _U[npoint1]; }//La = L[npoint1]; }
      else { Ua = INFINITE; }//La = INFINITE; }
      if (S[npoint2] == kDead) { Ub = _U[npoint2]; }//Lb = L[npoint2]; }
      else { Ub = INFINITE; }//Lb = INFINITE; }
      Vtmp = TsitsiklisQuadrantGradient(M1[point], M2[point], M3[point],
					Ua, Ub, _V[npoint1], _V[npoint2], i, q_gradient);
      if (q_gradient[0] < Ur)
        {
          Ur   =  q_gradient[0];
          //dUxr =  q_gradient[1];
          //dUyr =  q_gradient[2];
          //TsitsiklisQuadrantLength(Lr, 1, La, Lb,i);
          Vr = Vtmp;
          is_updated = true;
        }
    }
  //--------------------------------------------------------------
  if (is_updated)
    {
      _U[point] = Ur;
      //L[point] = Lr;
      _V[point] = Vr;
      //dUx[point] = dUxr;
      //dUy[point] = dUyr;
    }
  //--------------------------------------------------------------
  return is_updated;
}
*/

//================================================================
template<class T, class Integer>
Integer AnisoPropagation2D<T,Integer>::addSeed(Integer pos, Integer label)
{
  if (pos >= size) return 0;
  _U[pos] = 0.;
  //L[pos] = 0.;
  S[pos] = kSeed;
  _V[pos] = label;
  // add to heap
  tree->Push(pos);
  
  //if (_labels <= label) _labels = label + 1;
  return pos;
}

//================================================================
template<class T, class Integer>
Integer AnisoPropagation2D<T,Integer>::addSeed(int x, int y, Integer label)
{
  x += 1; y += 1;
  return addSeed(x + y*Nx, label);
}

//================================================================
template<class T, class Integer>
Integer AnisoPropagation2D<T,Integer>::setSeeds(double *pts, Integer nb_pts)
{
  Integer i, j, point, nb_inserted = 0;
  for (Integer s = 0; s < nb_pts; s++)
    {
      i = (Integer)round(pts[2*s]);
      j = (Integer)round(pts[1+2*s]);
      point = i + j*Nx;
      if ((point >= size) || (_U[point] == 0.)) continue;
      _U[point] = 0.;
      //L[point] = 0.0;
      S[point] = kOpen;
      _V[point] = _labels;
      _labels++;
      nb_inserted++;
      // add to heap
      tree->Push(point);
    }
  return nb_inserted;
}

//================================================================
/*
template<class T, class Integer>
void AnisoPropagation2D<T,Integer>::tsitsiklisUpdateNeighbors(const Integer &point)
{
  Integer npoint, state_np;
  for (short k = 0; k < connectivity; k++)
    {
      npoint = NPOINT(point, k);
      if ((state_np == S[npoint]) == kFar)
	{
	  tsitsiklisUpdate(npoint);
	  S[npoint] = kOpen;
	  tree->Push(npoint);
	  continue;
	}
      if (state_np = S[npoint])
	{
	  if (tsitsiklisUpdate(npoint) == true)
	    tree->Update((*tree)[npoint]);
	}
    }
}
*/

//================================================================
template<class T, class Integer>
  bool AnisoPropagation2D<T,Integer>::updateNeighbors(const Integer &point)
{
  Integer npoint, state_np, v;
  T u;
  for (short k = 0; k < connectivity; k++)
    {
      npoint = NPOINT(point, k);
      //--------------------------------------------------------------
      if ((state_np = S[npoint]) == kFar)
	{
	  _U[npoint] = update(npoint, point, k, _V[npoint]);
	  S[npoint] = kOpen;
	  tree->Push(npoint);
	  continue;
	}
      //--------------------------------------------------------------
      if (state_np == kOpen)
	{
	  if ((u = update(npoint, point, k, v)) < _U[npoint])
	    {
	      _V[npoint] = v;
	      _U[npoint] = u;
	      tree->Update((*tree)[npoint]);
	    }
	}
    }
  return true;
}

//================================================================
/*
  template<class T, class Integer>
  void AnisoPropagation2D<T,Integer>::propagate(T fmax = 0)
  {
  Integer point;
  //T Lmax = 0.0;
  //bool dmax_test = (dmax > 0);
  //--------------------------------------------------------------
  if (fmax > 0)
  {
     T nb_p = 0.;
     while ((tree->Size() > 0) && (fmax > nb_p/T(size)))
	{
	  point = tree->Pop();
	  S[point] = kDead;
	  nb_p += 1.;
	  tsitsiklisUpdateNeighbors(point);
	}
      return;
    }
  while (tree->Size() > 0)
    {
      point = tree->Pop();
      S[point] = kDead;
      tsitsiklisUpdateNeighbors(point);
      //(Lmax = std::max(Lmax, L[point])) > dmax) break;
    }
    }
*/

//================================================================
template<class T, class Integer>
  bool AnisoPropagation2D<T,Integer>::propagation(double *pts, Integer nb_pts, T fmax = 0)
{
  if (nb_pts < 1 || pts == NULL) return false;
  clearMaps();
  setSeeds(pts, nb_pts);
  propagation(this, &AnisoPropagation2D<T,Integer>::updateNeighbors);
  return true;
}

//================================================================
template<class T, class Integer>
Integer AnisoPropagation2D<T,Integer>::GridSizeX() { return Nx; }

//================================================================
template<class T, class Integer>
Integer AnisoPropagation2D<T,Integer>::GridSizeY() { return Ny; }

//================================================================
template<class T, class Integer>
void AnisoPropagation2D<T,Integer>::resize()
{
  Integer x, y, point, Point;
  for (y = 0; y < ny; y++)
    for (x = 0; x < nx; x++)
      {
	point = x+y*nx;
	Point = (x+1)+(y+1)*Nx;
	_U[point] = _U[Point];
	//L[point] = L[Point];
	//dUx[point] = dUx[Point];
	//dUy[point] = dUy[Point];
	_V[point] = _V[Point];
      }
}
//================================================================
template<class T, class Integer>
void AnisoPropagation2D<T,Integer>::pos_to_point(const Integer &pos,
						 Integer &x, Integer &y)
{
  y = pos / Nx;
  x = pos - Nx * y;
}

template<class T, class Integer>
Integer AnisoPropagation2D<T,Integer>::pos_to_point_y(const Integer &pos)
{
  return pos / Nx;
}

template<class T, class Integer>
short AnisoPropagation2D<T,Integer>::orient2D(Integer ax, Integer ay,
					      Integer bx, Integer by,
					      Integer cx, Integer cy)
{
    Integer d = (ax*(by-cy) - ay*(bx-cx) + (bx*cy - cx*by));
    return (d > 0 ? 1 : (d < 0 ? -1 : 0));
}

template<class T, class Integer>
short AnisoPropagation2D<T,Integer>::orient2D(const Integer &point1,
					      const Integer &point2,
					      const Integer &point3)
{
    Integer y1 = pos_to_point_y(point1),
      y2 = pos_to_point_y(point2),
      y3 = pos_to_point_y(point3);
    return orient2D(point1-Nx*y1, y1, point2-Nx*y2, y2, point3-Nx*y3, y3);
}

#endif
