//================================================================
//================================================================
// File: RecursiveAnisoPropagation2D.h
// (C) 01/2009 by Sebastien Bougleux
//================================================================
//================================================================
#ifndef __RECURSIVEANISOPROPAGATION2D_H__
#define __RECURSIVEANISOPROPAGATION2D_H__

#include "AnisoPropagation2D.h"

//================================================================
/**
 * @class RecursiveAnisoPropagation2D
 * @brief based on A Recursive Anisotropic Fast Marching Approach to Reaction Diffusion Equation: Application to Tumor Growth Modeling by E. Konukoglu, M. Sermesant, O. Clatz, J-M. Peyrat, H. Delingette and N. Ayache
 * @param T type of action values (must be a real type, float by default)
 * @param Integer type of indices, labels (can be a positive integer type, int by default)
 */
//================================================================
template<class T = float, class Integer = int>
class RecursiveAnisoPropagation2D : public AnisoPropagation2D<T,Integer>
//================================================================
{
  //================================================================
 protected:
  //--------------------------------------------------------------
  using AnisoPropagation2D<T,Integer>::size;
  using AnisoPropagation2D<T,Integer>::v8k;
  using AnisoPropagation2D<T,Integer>::_U;
  using AnisoPropagation2D<T,Integer>::_V;
  using AnisoPropagation2D<T,Integer>::S;
  using AnisoPropagation2D<T,Integer>::tree;
  using AnisoPropagation2D<T,Integer>::connectivity;
  using AnisoPropagation2D<T,Integer>::kDead;
  using AnisoPropagation2D<T,Integer>::kOpen;
  using AnisoPropagation2D<T,Integer>::kFar;
  using AnisoPropagation2D<T,Integer>::kBorder;
  using AnisoPropagation2D<T,Integer>::setSeeds;
  typedef typename AnisoPropagation2D<T,Integer>::Tree Tree;
  Tree *_changed;                     // min heap
  T EPSILON_UPDATE;
  //--------------------------------------------------------------
  // protected functions
  //--------------------------------------------------------------
  void clearMaps();
  //--------------------------------------------------------------
  //================================================================
 public:
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
  RecursiveAnisoPropagation2D(T* tf, int dx, int dy, T* U, Integer* V, 
			      T sx, T sy);
  //--------------------------------------------------------------
  /**
   * @brief Destructor
   */
  ~RecursiveAnisoPropagation2D();
  //--------------------------------------------------------------
  //typedef void (RecursiveAnisoPropagation2D<T,Integer>::*UpdateFctPtr)(const Integer&);
  bool updateNeighbors(const Integer&);
  template<class C> void propagation(C *cl, bool (C::*fct)(const Integer&))
  {
    Integer point;
    while (tree->Size() > 0 || _changed->Size() > 0)
      {
	if (_changed->Size() > 0)
	  {
	    point = _changed->Pop();
	  }
	else
	  {
	    point = tree->Pop();
	    S[point] = kDead;
	  }
	(cl->*fct)(point);
      }
  }
  bool propagation(double*, Integer, T, T);
};
//================================================================
//================================================================
//================================================================
//================================================================
// Functions
//================================================================
//================================================================
template<class T, class Integer>
RecursiveAnisoPropagation2D<T,Integer>::RecursiveAnisoPropagation2D(T* TF22, int _nx, int _ny, 
								    T *U, Integer *V = NULL,
								    T _hx = 1, T _hy = 1)
  : AnisoPropagation2D<T,Integer>(TF22, _nx, _ny, U, V, _hx, _hy),
    _changed(NULL), EPSILON_UPDATE(1.e-10)
{

}

//================================================================
template<class T, class Integer>
RecursiveAnisoPropagation2D<T,Integer>::~RecursiveAnisoPropagation2D()
{
  if (_changed) delete _changed;
}

//================================================================
template<class T, class Integer>
void RecursiveAnisoPropagation2D<T,Integer>::clearMaps()
{
  if (_changed) delete _changed;
  _changed = new Tree(size, _U);
  AnisoPropagation2D<T,Integer>::clearMaps();
}

//================================================================
template<class T, class Integer>
bool RecursiveAnisoPropagation2D<T,Integer>::updateNeighbors(const Integer &point)
{
  Integer npoint, state_np, v;
  T u;
  for (short k = 0; k < connectivity; k++)
    {
      npoint = NPOINT(point, k);
      if (S[npoint] == kDead)
	{
	  u = AnisoPropagation2D<T,Integer>::update(npoint, point, k, v);
	  if (_U[npoint] - u > EPSILON_UPDATE)
	    {
	      _U[npoint] = u;
	      _V[npoint] = v;
	      _changed->Push(npoint);
	    }
	}
    }
  AnisoPropagation2D<T,Integer>::updateNeighbors(point);
  return true;
}

//================================================================
template<class T, class Integer>
  bool RecursiveAnisoPropagation2D<T,Integer>::propagation(double *pts, Integer nb_pts, 
							   T eps = 1.e-10, 
							   T fmax = 0)
{
  if (nb_pts < 1 || pts == NULL) return false;
  EPSILON_UPDATE = eps;
  clearMaps();
  setSeeds(pts, nb_pts);
  propagation(this,&RecursiveAnisoPropagation2D<T,Integer>::updateNeighbors);
  return true;
}

#endif
