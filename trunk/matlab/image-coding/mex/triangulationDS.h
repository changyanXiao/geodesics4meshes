//================================================================
// File: triangulationsDS.h
// author: Sebastien Bougleux (sebastien.bougleux@unicaen.fr),
//         IUT Cherbourg-Manche / Universite de Caen
//         GREYC, Equipe Image
// date: 10 feb. 2009
// last modif: 21 march 2009
//================================================================
#include<iostream>
#include<map>
#include<list>
#include<queue>
#include<vector>
#include<cmath>
//================================================================
// returns 1 if ccw orientation, -1 if cw orientation, 
// or 0 if points are aligned
//================================================================
template<class T>
short orient2D(const T &ax, const T &ay,
               const T &bx, const T &by,
               const T &cx, const T &cy)
{
    T d = (ax - cx) * (by - cy) - (ay - cy) * (bx - cx);
    if (d > 0) return 1;
    if (d < 0) return -1;
    return 0;
}

//================================================================
double det4x4(double a00,  double a01,  double a02,  double a03,
              double a10,  double a11,  double a12,  double a13,
              double a20,  double a21,  double a22,  double a23,
              double a30,  double a31,  double a32,  double a33)
{
  // Compute det2x2
  double m01 = a10*a01 - a00*a11;
  double m02 = a20*a01 - a00*a21;
  double m03 = a30*a01 - a00*a31;
  double m12 = a20*a11 - a10*a21;
  double m13 = a30*a11 - a10*a31;
  double m23 = a30*a21 - a20*a31;
  // -------------------------------------------------------
  // Compute minors of rank 3
  double m012 = m12*a02 - m02*a12 + m01*a22;
  double m013 = m13*a02 - m03*a12 + m01*a32;
  double m023 = m23*a02 - m03*a22 + m02*a32;
  double m123 = m23*a12 - m13*a22 + m12*a32;
  // -------------------------------------------------------
  // Compute minors of rank 4
  return m123*a03 - m023*a13 + m013*a23 - m012*a33;
}

//================================================================
// returns true if p4 is in the circle passing by p1, p2 and p3
// suppose points p1, p2 and p3 in ccw order
//================================================================
template<class T>
bool _inCircle2D(const T &x1, const T &y1,
                 const T &x2, const T &y2,
                 const T &x3, const T &y3,
                 const T &x4, const T &y4)
{
   //T ca = (x3 - x2) * (x1 - x2) + (y3 - y2) * (y1 - y2);
   //T cb = (x1 - x4) * (x3 - x4) + (y1 - y4) * (y3 - y4);
   //if ((ca < 0) && (cb < 0)) return true;
   //if ((ca > 0) && (cb > 0)) return false;
   //T sa = (x3 - x2) * (y1 - y2) - (x1 - x2) * (y3 - y2);
   //T sb = (x1 - x4) * (y3 - y4) - (x3 - x4) * (y1 - y4);
   //std::cerr << sa << " " << sb << std::endl;
   //if ((ca * sb + sa * cb) < 0) return true;
   //return false;
   double det = det4x4(x1, y1, (x1*x1+y1*y1), 1,
                       x2, y2, (x2*x2+y2*y2), 1,
                       x3, y3, (x3*x3+y3*y3), 1,
                       x4, y4, (x4*x4+y4*y4), 1);
   
   return (det > 0);
}

//================================================================
// returns true if p4 is in the circle passing by p1, p2 and p3
//================================================================
template<class T>
bool inCircle2D(const T &x1, const T &y1,
                const T &x2, const T &y2,
                const T &x3, const T &y3,
                const T &x4, const T &y4)
{
    if (orient2D(x1, y1, x2, y2, x3, y3) > 0) // ccw case
        return _inCircle2D(x1, y1, x2, y2, x3, y3, x4, y4);
    // cw case
    return _inCircle2D(x3, y3, x2, y2, x1, y1, x4, y4);
}
//================================================================
// class Pair
//================================================================
template<class T>
class Pair
{
  public:
    T idx1;
    T idx2;
    Pair(const T &i1, const T &i2) : idx1(i1), idx2(i2) { if (idx1 > idx2) std::swap(idx1,idx2); }
    ~Pair() { }
    Pair(const Pair<T> &p) : idx1(p.idx1), idx2(p.idx2) { }
    bool operator<(const Pair<T> &p) const
      {
        if (idx1 < p.idx1) return true;
        if (idx1 == p.idx1 && idx2 < p.idx2) return true;
        return false;
      }
    bool operator>(const Pair<T> &p) const
      {
        if (idx1 > p.idx1) return true;
        if (idx1 == p.idx1 && idx2 > p.idx2) return true;
        return false;
      }
    bool operator==(const Pair<T> &p) const
      {
        return (idx1 == p.idx1 && idx2 == p.idx2);
      }
};
//================================================================
// class Edge
//================================================================
template<class T, class FT>
class Edge
{
    public:
        T idx1; // extremity
        T idx2; // extremity
        T idx3; // opposite
        T idx4; // opposite
        T val;  // for edge flip
        //FT length; // of the edge
        Edge() : idx1(0), idx2(0), idx3(0), idx4(0), val(0)//, length(0)
        {
            
        }
        Edge(const T &i1, const T &i2) : idx1(i1), idx2(i2), idx3(0), idx4(0), val(0){ }//, length(0) { }
        Edge(const T &i1, const T &i2, const T &i3) : idx1(i1), idx2(i2), idx3(i3), idx4(0), val(0) { }//, length(0) { }
        //Edge(const T &i1, const T &i2, const T &i3, FT d) : idx1(i1), idx2(i2), idx3(i3), idx4(0), val(0) { }//, length(d) { }
        ~Edge() { }
        /*bool operator<(const Edge<T,FT> &e) const
        {
            return (length < e.length);
        }
        bool operator>(const Edge<T,FT> &e) const
        {
            return (length > e.length);
        }
        bool operator==(const Edge<T,FT> &e) const
        {
            return (length == e.length);
        }*/
            
};

/*template<class T, class FT>
class EdgePtrComparison
{
 public:
    typedef Edge<T,FT>* EdgePtr; 
    EdgePtrComparison() { }
    bool operator() (const EdgePtr &e1, const EdgePtr &e2) const
    {
      return (*e1 < *e2);
    }
};*/

//================================================================
// class Triangulation
//================================================================
template<class T, class FT>
class Triangulation
{
    public:
        typedef Edge<T,FT> Edge;
        typedef typename std::map<Pair<T>,Edge*> EdgeMap;
        typedef typename EdgeMap::iterator EdgeIterator;
        typedef typename std::list<Pair<T> > FlipList;
        typedef typename FlipList::iterator FlipListIterator;
        typedef typename std::queue<Edge*> EdgeQueue;
        //typedef typename std::priority_queue<Edge*,std::vector<Edge*>,EdgePtrComparison<T,FT> > EdgeQueue;
    protected:
        EdgeMap *_emap;
        EdgeQueue *_edgeQ;
        FT *_pointsX;
        FT *_pointsY;
        T _nb_pts;
        FlipList *_flip_list;
        //================================================================
        bool isNotLocallyDelaunay(Edge *e)
        {
            if ((e->idx3 == 0) || (e->idx4 == 0)) return false;       // boundary case
            T x1 = (int)point_x(e->idx1), y1 = (int)point_y(e->idx1);
            T x2 = (int)point_x(e->idx2), y2 = (int)point_y(e->idx2);
            T x3 = (int)point_x(e->idx3), y3 = (int)point_y(e->idx3);
            T x4 = (int)point_x(e->idx4), y4 = (int)point_y(e->idx4);
            // -------------------------------------------------------
            if ((orient2D(x3, y3, x1, y1, x2, y2) != orient2D(x3, y3, x1, y1, x4, y4)) ||
                (orient2D(x3, y3, x2, y2, x1, y1) != orient2D(x3, y3, x2, y2, x4, y4)))
                return false;   // e diagonal of a non-convex quad
            // -------------------------------------------------------
            if (inCircle2D<T>(x1, y1, x3, y3, x2, y2, x4, y4) == true) return true;
            return false;
        }
        //================================================================
        void flipEdge(Edge *e)
        {
            // flip e=(1,2) with (3,4)
            T t1 = e->idx1, t2 = e->idx2;
            e->idx1 = e->idx3;
            e->idx2 = e->idx4;
            e->idx3 = t1;
            e->idx4 = t2;
        }
        //================================================================
        bool insertOrUpdateEdge(const T &i1, const T &i2, const T &i3)
        {
            Pair<T> p(i1, i2);
            EdgeIterator eit = _emap->find(p);
            if (eit == _emap->end())
            {
                _emap->insert(std::make_pair(p, new Edge(p.idx1, p.idx2, i3)));
                return false;
            }
            eit->second->idx4 = i3;                
            return true;
        }
        //================================================================
        bool insertOrUpdateEdgeAndFlipDetect(const T &i1, const T &i2, const T &i3)
        {
            Pair<T> p(i1, i2);
            EdgeIterator eit = _emap->find(p);
            if (eit == _emap->end())
            {
                _emap->insert(std::make_pair(p, new Edge(p.idx1, p.idx2, i3)));
                return false;
            }
            eit->second->idx4 = i3;
            //eit->second->length = distance(eit->second->idx3, i3);
            // -------------------------------------------------------
            if (isNotLocallyDelaunay(eit->second) == true)
            {
                _edgeQ->push(eit->second);
                eit->second->val = 1;
                return true;
            }
            return false;
        }
        //================================================================
        void face2Edges(const T &i1, const T &i2, const T &i3)
        {
            insertOrUpdateEdge(i1, i2, i3);
            insertOrUpdateEdge(i1, i3, i2);
            insertOrUpdateEdge(i2, i3, i1);
        }
        //================================================================
        void face2EdgesAndFlipDetect(const T &i1, const T &i2, const T &i3)
        {
            insertOrUpdateEdgeAndFlipDetect(i1, i2, i3);
            insertOrUpdateEdgeAndFlipDetect(i1, i3, i2);
            insertOrUpdateEdgeAndFlipDetect(i2, i3, i1);
        }
    public:
        //================================================================
        Triangulation() : _emap(0), _edgeQ(0), _pointsX(0), _pointsY(0), _nb_pts(0), _flip_list(0) { }
        //================================================================
        Triangulation(T *vec_faces, T size_vec) : _emap(0), _edgeQ(0), _pointsX(0), _pointsY(0), _nb_pts(0), _flip_list(0)
        {
            loadFaces(vec_faces, size_vec);
        }
        //================================================================
        ~Triangulation()
        {
            if (_emap) delete _emap;
            if (_flip_list) delete _flip_list;
            if (_pointsX) delete[] _pointsX;
            if (_pointsY) delete[] _pointsY;
            if (_edgeQ) delete _edgeQ;
        }
        //================================================================
        FT point_x(T idx)
        {
            return _pointsX[idx-1];
        }
        //================================================================
        FT point_y(T idx)
        {
            return _pointsY[idx-1];
        }
        //================================================================
        void loadFaces(T *vec_faces, T size_vec, FT *points, T nb_pts)
        {
            if (_emap) delete _emap;
            _emap = new EdgeMap;
            T i1, i2, i3;
            for (T f = 0; f < size_vec;)
            {
                i1 = vec_faces[f++];
                i2 = vec_faces[f++];
                i3 = vec_faces[f++];
                face2Edges(i1, i2, i3);
            }
        }
        //================================================================
        FT distance(const T &i1, const T &i2)
        {
            FT dx = point_x(i1) - point_x(i2);
            FT dy = point_y(i1) - point_y(i2);
            return std::sqrt(dx*dx + dy*dy);
            
        }
        //================================================================
        bool updateEdge(const T &i1, const T &i2, const T &i3, const T &i4)
        {
            Pair<T> p(i1, i2);
            EdgeIterator eit;
            if ((eit = _emap->find(p)) != _emap->end())
            {
                if (eit->second->idx3 == i3)
                {
                    eit->second->idx3 = i4;
                }
                else
                {
                    eit->second->idx4 = i4;
                }
                //eit->second->length = distance(eit->second->idx3, eit->second->idx4);
                if (isNotLocallyDelaunay(eit->second) == true)
                {
                    if (eit->second->val == 0)
                    {
                        eit->second->val = 1;
                        _edgeQ->push(eit->second);
                    }
                    
                }
                else eit->second->val = 0;
                return true;
            }
            return false;
        }
        //================================================================
        void loadFacesAndFlip(T *vec_faces, T size_vec, FT *points, T nb_pts)
        {
            if (_emap) delete _emap;
            _emap = new EdgeMap;
            if (_flip_list) delete _flip_list;
            _flip_list = new FlipList;
            if (_pointsX) delete[] _pointsX;
            _pointsX = new FT[nb_pts];
            if (_pointsY) delete[] _pointsY;
            _pointsY = new FT[nb_pts];
            if (_edgeQ) delete _edgeQ;
            _edgeQ = new EdgeQueue;
            // -------------------------------------------------------
            // import faces and prepare flips
            T i1, i2, i3;
            FT px, py;
            T i, j;
            for (i = 0, j = 0; i < nb_pts; i++, j += 2)
            {
                _pointsX[i] = points[j];
                _pointsY[i] = points[j+1];
            }
            for (T f = 0; f < size_vec;)
            {
                i1 = vec_faces[f++];
                i2 = vec_faces[f++];
                i3 = vec_faces[f++];
                face2EdgesAndFlipDetect(i1, i2, i3);
            }
            // -------------------------------------------------------
            // Flip non-Delaunay edges
            Edge *e = 0;
            EdgeIterator eit;
            while (!_edgeQ->empty())
            {
                e = _edgeQ->front();
                //e = _edgeQ->top();
                //std::cerr << e->length << std::endl;
                if (e->val == 0)
                {
                    _edgeQ->pop();
                    continue;
                }
                Pair<T> p(e->idx1, e->idx2);
                eit = _emap->find(p);
                flipEdge(e);
                e->val = 0;
                Pair<T> q(e->idx1, e->idx2);
                _flip_list->push_back(q); // save into the fliped edge list
                _emap->erase(eit);
                //e->length = distance(e->idx3, e->idx4);
                _emap->insert(std::make_pair(q, e));
                // -------------------------------------------------------
                // update adjacent edges
                updateEdge(e->idx1, e->idx3, e->idx4, e->idx2);
                updateEdge(e->idx1, e->idx4, e->idx3, e->idx2);
                updateEdge(e->idx2, e->idx3, e->idx4, e->idx1);
                updateEdge(e->idx2, e->idx4, e->idx3, e->idx1);
                _edgeQ->pop();
            }
            
            
        }
        //================================================================
        void getFlipList(FT *vec)
        {
            FlipListIterator fit = _flip_list->begin(), fend = _flip_list->end();
            T i = 0;
            for (; fit != fend; ++fit)
            {
                vec[i++] = fit->idx1;
                vec[i++] = fit->idx2;
            }
        }
        //================================================================
        T nbFlipEdges() { return _flip_list->size(); }
        //================================================================
        bool getEdges(FT *vec)
        {
            EdgeIterator eit = _emap->begin(), end = _emap->end();
            T i = 0;
            for (; eit != end; ++eit)
            {
                vec[i++] = eit->second->idx1;
                vec[i++] = eit->second->idx2;
            }
        }
        //================================================================
        T nbEdges() { return _emap->size(); }
        
        
};

