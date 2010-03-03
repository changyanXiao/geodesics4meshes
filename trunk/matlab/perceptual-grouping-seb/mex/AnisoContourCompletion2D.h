//================================================================
//================================================================
// File: AnisoPropagation2D.h
// (C) 12/2008 by Sebastien Bougleux
//================================================================
//================================================================
#ifndef __ANISOCONTOURCOMPLETION2D_H__
#define __ANISOCONTOURCOMPLETION2D_H__

#include "RecursiveAnisoPropagation2D.h"
#include<map>
//================================================================
/**
 * @class AnisoContourCompletion2D
 * @brief 
 * @param T type of action values (must be a real type, float by default)
 * @param Integer type of indices, labels (can be a positive integer type, int by default)
 */
//================================================================
template<class Integer>
struct DoublePoint
{
  Integer l1;  // first Voronoi label (the smallest one)
  Integer l2;  // second Voronoi label
  //TriplePoint *tpoint_cw;
  //TriplePoint *tpoint_ccw;
  DoublePoint() : l1(0), l2(0) { } //, tpoint_cw(NULL), tpoint_ccw(NULL) { }
  DoublePoint(Integer _l1, Integer _l2) : l1(_l1), l2(_l2)
  {
        if (l2 < l1) std::swap(l1,l2);
  }
  DoublePoint(const DoublePoint &e) : l1(e.l1), l2(e.l2) { } //, tpoint_cw(e.tpoint_cw), tpoint_ccw(e.tpoint_ccw) {}
  bool operator<(const DoublePoint &e) const
  {
    if (l1 == e.l1) return (l2 < e.l2);
    return (l1 < e.l1);
  }
  bool operator==(const DoublePoint &e) const
  {
    return ((l1 == e.l1) && (l2 == e.l2));
  }
  void init(Integer _i1, Integer _i2)
  {
    if (_i1 < _i2) { l1 = _i1; l2 = _i2; }
    else { l1 = _i2; l2 = _i1; }
  }
};
//================================================================
template<class Integer, class T>
struct Edge
{
  Integer idx1;
  Integer idx2;
  Integer point1;
  Integer point2;
  T d1;
  T d2;
  T maxd;
  //Edge() : point1(0), point2(0), d1(-1), d2(-1) { }
  //Edge(const float &_d) : idx1(-1), idx2(-1), point1(-1), point2(-1), d(_d) { }
  Edge(Integer i1, Integer i2, Integer p1, Integer p2, T da, T db) : maxd(std::max(da,db))
  {
    if (i1 < i2) { idx1 = i1; idx2 = i2; point1 = p1; point2 = p2; d1 = da; d2 = db;}
    else { idx1 = i2; idx2 = i1; point1 = p2; point2 = p1; d1 = db; d2 = db; }
  }
  ~Edge() { }
  bool update(const Integer &va, const T &da, const Integer &pointa,
	      const Integer &vb, const T &db, const Integer &pointb)
  {
    T maxab = std::max(da,db);
    if (maxab < maxd)
      {
	maxd = maxab;
	if (va == idx1) { point1 = pointa; d1 = da; point2 = pointb; d2 = db; return true; }
	if (va == idx2) { point2 = pointa; d2 = da; point1 = pointb; d1 = db; return true; }
      }
    return false;
  }
};
//================================================================

template<class T = float, class Integer = int>
class AnisoContourCompletion2D : public RecursiveAnisoPropagation2D<T,Integer>
//================================================================
{
  //================================================================
 protected:
  //--------------------------------------------------------------
  using AnisoPropagation2D<T,Integer>::size;
  using AnisoPropagation2D<T,Integer>::Nx;
  using AnisoPropagation2D<T,Integer>::Ny;
  using AnisoPropagation2D<T,Integer>::_labels;
  using AnisoPropagation2D<T,Integer>::v8k;
  using AnisoPropagation2D<T,Integer>::_U;
  using AnisoPropagation2D<T,Integer>::_V;
  using AnisoPropagation2D<T,Integer>::S;
  using AnisoPropagation2D<T,Integer>::tree;
  using RecursiveAnisoPropagation2D<T,Integer>::_changed;
  using AnisoPropagation2D<T,Integer>::connectivity;
  using AnisoPropagation2D<T,Integer>::kDead;
  using AnisoPropagation2D<T,Integer>::kOpen;
  using AnisoPropagation2D<T,Integer>::kFar;
  using AnisoPropagation2D<T,Integer>::kBorder;
  using AnisoPropagation2D<T,Integer>::setSeeds;
  typedef typename AnisoPropagation2D<T,Integer>::Tree Tree;
  typedef std::map<DoublePoint<Integer>, Edge<Integer,T>* > EdgeMap;
  typedef typename EdgeMap::iterator EdgeMapIterator;
  EdgeMap _edges;
  std::map<Integer,Integer> _seeds;
  Integer *_degree;
  Integer _nb_edges;
  //--------------------------------------------------------------
  // protected functions
  //--------------------------------------------------------------
  Integer setSeeds(double*, Integer);
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
  AnisoContourCompletion2D(T* tf, int dx, int dy, T* U, Integer* V, 
			   /*T*, T*, T*, */ T sx, T sy);
  //--------------------------------------------------------------
  /**
   * @brief Destructor
   */
  ~AnisoContourCompletion2D();
  //--------------------------------------------------------------
  bool updateNeighbors_NF_REC(const Integer&);
  bool updateNeighbors_NF(const Integer&);
  bool greedyNearestFronts(double*, Integer, bool);
  void printEdges();
  void getSeeds(T *tab)
  {
    if (tab == NULL) return;
    typename std::map<Integer,Integer>::iterator it = _seeds.begin(), end = _seeds.end();
    Integer i = 0, x, y;
    for (; it != end; ++it)
      {
	pos_to_point(it->second, x, y);
	tab[i++] = (T)x;
	tab[i++] = (T)y;
      }
  }
  Integer nbEdges() { return _edges.size(); }
  void getEdges(Integer *tab)
  {
    if (tab == NULL) return;
    EdgeMapIterator eit = _edges.begin(), eend = _edges.end();
    Integer i = 0;
    for (; eit != eend; ++eit)
    {
      tab[i++] = eit->second->idx1;
      tab[i++] = eit->second->idx2;
    }
  }
  Integer nbSaddlePoints() { return 2*_edges.size(); }
  void getSaddlePoints(T *tab)
  {
    if (tab == NULL) return;
    EdgeMapIterator eit = _edges.begin(), eend = _edges.end();
    Integer i = 0, x1, y1, x2, y2;
    for (; eit != eend; ++eit)
    {
      pos_to_point(eit->second->point1, x1, y1);
      pos_to_point(eit->second->point2, x2, y2);
      tab[i++] = (T)x1;
      tab[i++] = (T)y1;
      tab[i++] = (T)x2;
      tab[i++] = (T)y2;
    }
  }
};
//================================================================
//================================================================
//================================================================
//================================================================
// Functions
//================================================================
//================================================================
template<class T, class Integer>
AnisoContourCompletion2D<T,Integer>::AnisoContourCompletion2D(T* TF22, int _nx, int _ny, 
							      T *U, Integer *V = NULL,
							      T _hx = 1, T _hy = 1)
  : RecursiveAnisoPropagation2D<T,Integer>(TF22, _nx, _ny, U, V, _hx, _hy),
    _degree(NULL), _nb_edges(0)
{
  
}

//================================================================
template<class T, class Integer>
AnisoContourCompletion2D<T,Integer>::~AnisoContourCompletion2D()
{
  
}

//================================================================
template<class T, class Integer>
void AnisoContourCompletion2D<T,Integer>::clearMaps()
{
  RecursiveAnisoPropagation2D<T,Integer>::clearMaps();
}

//================================================================
template<class T, class Integer>
Integer AnisoContourCompletion2D<T,Integer>::setSeeds(double *pts, Integer nb_pts)
{
  Integer i, j, point, nb_inserted = 0;
  for (Integer s = 0; s < nb_pts; s++)
    {
      i = (Integer)round(pts[2*s]);
      j = (Integer)round(pts[1+2*s]);
      point = i + j*Nx;
      if ((point >= size) || (_U[point] == 0.)) continue;
      _U[point] = 0.;
      S[point] = kOpen;
      _V[point] = _labels;
      _seeds.insert(std::make_pair(_labels, point));
      _labels++;
      nb_inserted++;
      tree->Push(point);       // add to heap
    }
  return nb_inserted;
}

//================================================================
template<class T, class Integer>
bool AnisoContourCompletion2D<T,Integer>::updateNeighbors_NF_REC(const Integer &point)
{
  Integer npoint, state_np, v, vpoint = _V[point];
  T u;
  EdgeMapIterator eit;
  DoublePoint<Integer> dp;
  for (short k = 0; k < connectivity; k++)
    {
      npoint = NPOINT(point, k);
      if (S[npoint] == kDead)
	{
	  //--------------------------------------------------------------
	  // new action is lower
	  u = AnisoPropagation2D<T,Integer>::update(npoint, point, k, v);
	  if (_U[npoint] - u > 1.e-10)
	    {
	      _U[npoint] = u;
	      _V[npoint] = v;
	      _changed->Push(npoint);
	    }
	  //--------------------------------------------------------------
	  // first contact of fronts ?
	  if ((v = _V[npoint]) != vpoint)
	    {
	      dp.init(vpoint, v);
	      //--------------------------------------------------------------
	      if ((eit = _edges.find(dp)) == _edges.end())
		{
		  if ((_nb_edges < _seeds.size()) && (_degree[vpoint] < 2) && (_degree[v] < 2))
		    {
		      _edges.insert(std::make_pair(dp, new Edge<Integer,T>(vpoint, v, point, npoint, 
									   _U[point], _U[npoint])));
		      _nb_edges++;
		      _degree[vpoint]++;
		      _degree[v]++;
		    }
		}
	      //--------------------------------------------------------------
	      else
		{
		  eit->second->update(v, _U[npoint], npoint, vpoint, _U[point], point);
		  //if (eit->second->update(v, _U[npoint], npoint, vpoint, _U[point], point) == true)
		  //{
		  //std::cerr << "update\n";
		  //  }
		}
	    }
	}
    }
  AnisoPropagation2D<T,Integer>::updateNeighbors(point);
  return (_nb_edges < _seeds.size());
}
//================================================================
template<class T, class Integer>
bool AnisoContourCompletion2D<T,Integer>::updateNeighbors_NF(const Integer &point)
{
  Integer npoint, state_np, v, vpoint = _V[point];
  EdgeMapIterator eit;
  T u;
  DoublePoint<Integer> dp;
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
	  continue;
	}
      //--------------------------------------------------------------
      if (S[npoint] == kDead)
	{
	  if ((v = _V[npoint]) != vpoint)
	    {
	      dp.init(vpoint, v);
	      //--------------------------------------------------------------
	      if ((eit = _edges.find(dp)) == _edges.end())
		{
		  if ((_nb_edges < _seeds.size()) && (_degree[vpoint] < 2) && (_degree[v] < 2))
		    {
		      _edges.insert(std::make_pair(dp, new Edge<Integer,T>(vpoint, v, point, npoint, 
									   _U[point], _U[npoint])));
		      _nb_edges++;
		      _degree[vpoint]++;
		      _degree[v]++;
		    }
		}
	      //--------------------------------------------------------------
	      else
		{
		  eit->second->update(v, _U[npoint], npoint, vpoint, _U[point], point);
		  //if (eit->second->update(v, _U[npoint], npoint, vpoint, _U[point], point) == true)
		  //{
		  //std::cerr << "update\n";
		  //}
		}
	    }
	}
    }
  return (_nb_edges < _seeds.size());
}
//================================================================
template<class T, class Integer>
  bool AnisoContourCompletion2D<T,Integer>::greedyNearestFronts(double *pts, Integer nb_pts, bool rec = false)
{
  if (nb_pts < 2 || pts == NULL) return false;
  clearMaps();
  _degree = new Integer[nb_pts+1];
  for (Integer i = 0; i <= nb_pts; i++) _degree[i] = 0;
  setSeeds(pts, nb_pts);
  if (rec == false)
    AnisoPropagation2D<T,Integer>::propagation(this,&AnisoContourCompletion2D<T,Integer>::updateNeighbors_NF);
  else
    RecursiveAnisoPropagation2D<T,Integer>::propagation(this,&AnisoContourCompletion2D<T,Integer>::updateNeighbors_NF_REC);
  return true;
}

//================================================================
template<class T, class Integer>
  void AnisoContourCompletion2D<T,Integer>::printEdges()
{
  EdgeMapIterator eit = _edges.begin(), eend = _edges.end();
  int x1, y1, x2, y2;
  for (; eit != eend; ++eit)
    {
      pos_to_point(eit->second->point1, x1, y1);
      pos_to_point(eit->second->point2, x2, y2);
      std::cerr << "[" << x1 << ";" << y1 << "] [" 
		<< x2 << ";" << y2 << "] (" 
		<< eit->second->idx1 << "," << eit->second->idx2
		<< ")" << std::endl;
    }
  std::cerr << std::endl;
}

#endif
