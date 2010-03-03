// 

#include "fm.h"
#include<vector>
#include<list>
#include<queue>
#include<iostream>
#include<map>

short kDead = -1;
short kOpen = -2;
short kFar = -3;
short kBorder = -4;
short kSeed = -5;

/* Global variables */
int nx;			// real size on X
int ny;			// real size on Y
int nxny;
int Nx, Ny;     // size for computing
int size;       // number of points
float _labels;
float hx, hy;   // spacing
float hx2, hy2, hxhy;
float hx2hy2;
float hx2_plus_hy2;
float sqrt_hx2_plus_hy2;
float* U = NULL;// action map
float*  dUx = NULL;
float*  dUy = NULL;
//float* L = NULL; // distance map
short* S = NULL; // states
float* V = NULL; // voronoi
float* M1 = NULL; // M1
float* M2 = NULL; // M1
float* M3 = NULL; // M1
double* T = NULL;
//================================================================
int connectivity_large;
int* NeighborhoodLarge = NULL;
int* v8c = NULL;
int* Simplicies1 = NULL;
int* Simplicies2 = NULL;
int* signsX      = NULL;
int* signsY      = NULL;
float* h1 = NULL;
float* h2 = NULL;
float* h1_h2 = NULL;
float* h22 = NULL;
float* q_gradient = NULL;
//================================================================
float w; // regularization parameter
float Dmax;
double* start_points = NULL;
int nb_start_points;
int start_points_size;
// min-heap
int*	  Tree;
int*	  Trial;
int*	  PtrToTrial;
std::map<float,int> *Seeds = NULL;
//================================================================
struct Edge
{
  float label1;
  float label2;
  int midpoint;
  bool mark;
  //std::map<Integer,Integer> split_points;
  std::queue<int> points;
  Edge(const float &p1, const float &p2) : label1(p1), label2(p2), midpoint(0), mark(false)
  {
    //if (label1 > label2) std::swap(label1, label2);
  }
  Edge(const float &p1, const float &p2, const int &mdp)
    : label1(p1), label2(p2), midpoint(mdp), mark(false)
  {
    //if (label1 > label2) std::swap(label1, label2);
  }
  //Edge(const Edge &e)
    //: label1(e.label1), label2(e.label2), midpoint(e.midpoint), mark(e.mark)
  //{ }
  
  ~Edge() { }
  //--------------------------------------------------------------
  //bool operator==(const Edge &e) const
  //{
   // return ((e.label1 == label1) && (e.label2 == label2));
  //}
  /*bool operator<(const Edge &e) const
  {
    return ((label1 < e.label1) || ((label1 == e.label1) && (label2 < e.label2)));
  }
  bool operator>(const Edge &e) const
  {
    return ((label1 > e.label1) || ((label1 == e.label1) && (label2 > e.label2)));
  }
  bool addSplitPoint(const int &point_label, const int &point)
  {
      typename std::map<int, int>::iterator m_it = split_points.find(point_label);
      if (m_it == split_points.end())
      {
          split_points.insert(std::make_pair(point_label, point));
          return true;
      }
      return false;
  }*/
};
//================================================================
std::queue<Edge*> *_edges_to_split = NULL;
//================================================================
struct VoronoiVertex
{
    int point;
    float label;
    float u;
    VoronoiVertex(const int &p, const float &l, const float &d)
    : point(p), label(l), u(d)
    {
        
    }
    VoronoiVertex(const VoronoiVertex &v)
    : point(v.point), label(v.label), u(v.u)
    {
        
    }
    ~VoronoiVertex() { }
    bool operator==(const VoronoiVertex &v) const
      {
        return (v.u == u);
      }
    bool operator<(const VoronoiVertex &v) const
      {
        return (u < v.u);
      }
    bool operator>(const VoronoiVertex &v) const
      {
        return (u > v.u);
      }
};
//================================================================
class VoronoiVertexComparison
{
 public:
  VoronoiVertexComparison() { }
  bool operator() (const VoronoiVertex &v1, const VoronoiVertex &v2) const
    {
      return (v1 < v2);
    }
};
//================================================================
typedef std::priority_queue<VoronoiVertex, std::vector<VoronoiVertex>,
                            VoronoiVertexComparison> VoronoiQueue;
VoronoiQueue *_inverted_triangles = NULL;
//================================================================
struct Triangle
{
    float l1;
    float l2;
    float l3;
    int voronoi_vertex;
    Triangle(const float &label1, const float &label2, const float &label3, const int &vv)
    : l1(label1), l2(label2), l3(label3), voronoi_vertex(vv)
    {
        
    }
    Triangle(const Triangle &t)
    : l1(t.l1), l2(t.l2), l3(t.l3), voronoi_vertex(t.voronoi_vertex)
    {
        
    }
    ~Triangle() { }
    
};
//================================================================
std::list<Triangle> *_triangles = NULL;
//--------------------------------------------------------------
void print_point(int p)
//--------------------------------------------------------------
{
    int y = p / Nx;
    int x = p - Nx * y;
    std::cerr << "[" << x << ";" << y << "]";
}
//================================================================
// MIN-HEAP
//================================================================
//================================================================
int Tree_GetFather(int position)
//================================================================
{
  if(position) {
	if (ceil((float)position/2)==((float)position/2))
	  return (position/2 -1);
	else return((position-1)/2);
  }
  else return -1;
};

//================================================================
// Tree_PushIn
/*
COMMENTS : 
The tree is updated, since the value at the bottom is no longer 
necessarily greater than the value of its father. Starting from 
the bottom of the tree (the node we have just pushed in), the 
value of each node is compared with the value of its father. If 
the value of the father is greater, the two nodes are permuted. 
The process stops when the top of the tree has been reached.
*/
//================================================================
int Tree_PushIn(int NewPosition)
{
  *(++PtrToTrial) = NewPosition;
  int position = (int)(PtrToTrial - Trial);
  Tree[NewPosition] = position;
  int s_idx = Trial[position];
  int f_idx = Trial[Tree_GetFather(position)];
  while((position != 0) && (U[s_idx] < U[f_idx]))
  {
	int buffer = Trial[position];				
	Tree[Trial[position]] = Tree_GetFather(position);
	Trial[position] = Trial[Tree_GetFather(position)];
	Tree[Trial[Tree_GetFather(position)]] = position;
	Trial[Tree_GetFather(position)] = buffer;
	position = Tree_GetFather(position);	
	s_idx = Trial[position];
	f_idx = Trial[Tree_GetFather(position)];
  }
  return (PtrToTrial - Trial +1);
};

//================================================================
bool Tree_IsEmpty()
//================================================================
{ 
  return ((PtrToTrial - Trial + 1) == 0);
};

//================================================================
int Tree_GetRightSon(int position)
//================================================================
{
  if ((2*position+2) < (int)(PtrToTrial - Trial +1))
	return 2*position+2 ;
  else return -1;
};

//================================================================
int Tree_GetLeftSon(int position)
//================================================================
{
  if ((2*position+1) < (int)(PtrToTrial - Trial +1))
	return 2*position+1 ;
  else return -1;
};

//================================================================
// Tree_UpdateDescent
/*
COMMENTS : 
The tree is updated in order to extract the head by marching down
the tree. Starting from the head of the tree, the value of a 
node is compared with the values of its two sons and replaced by
the smallest one. This process is iterated from the son with the
smallest value, until a leaf has been reached.
*/
//================================================================
void Tree_UpdateDescent()
{
  int position = 0;
  bool stop = false;
  while((position >= 0)&&(!stop)) {		
    if((Tree_GetRightSon(position) > 0) && (Tree_GetLeftSon(position) > 0))
	{
	  int ls_idx = Trial[Tree_GetLeftSon(position)];
	  int rs_idx = Trial[Tree_GetRightSon(position)];
	  if(U[ls_idx] <= U[rs_idx]) {
		Trial[position] = Trial[Tree_GetLeftSon(position)];
		Tree[Trial[position]] = position;
		position = Tree_GetLeftSon(position);
	  }
	  else {
		Trial[position] = Trial[Tree_GetRightSon(position)];				
		Tree[Trial[position]] = (position);
		position = Tree_GetRightSon(position);				
	  }
	}
	else
        if (Tree_GetLeftSon(position) > 0) {
            Trial[position] = Trial[Tree_GetLeftSon(position)];				
            Tree[Trial[position]] = (position);
            position = Tree_GetLeftSon(position);
        }
        else stop = true;
    }
    if(position != (PtrToTrial - Trial)) {
        Tree[*PtrToTrial] = position;
        Trial[position]=*PtrToTrial;
        int s_idx = Trial[position];
        int f_idx = Trial[Tree_GetFather(position)];
        while((position!=0)&&(U[s_idx]<U[f_idx])) {
            int buffer = Trial[position];				
            Tree[Trial[position]] = Tree_GetFather(position);
            Trial[position] = Trial[Tree_GetFather(position)];		
            Tree[Trial[Tree_GetFather(position)]] = position;
            Trial[Tree_GetFather(position)] = buffer;
            position = Tree_GetFather(position);			
            s_idx = Trial[position];
            f_idx = Trial[Tree_GetFather(position)];
        }
    }
};

//================================================================
int Tree_PopHead()
//================================================================
{
  if(PtrToTrial - Trial + 1)
  {
	int first = *Trial;
	Tree[first] = -1;
	Tree_UpdateDescent();
	PtrToTrial--;
	return first;
  }
  else return NULL;
};

//================================================================
void Tree_UpdateChange(int position)
//================================================================
{
  int s_idx = Trial[position];
  int f_idx = Trial[Tree_GetFather(position)];
  while ((position!=0) && (U[s_idx]<U[f_idx]))
  {
	int buffer = Trial[position];				
	Tree[Trial[position]] = Tree_GetFather(position);
	Trial[position] = Trial[Tree_GetFather(position)];
	Tree[Trial[Tree_GetFather(position)]] = position;
	Trial[Tree_GetFather(position)] = buffer;
	position = Tree_GetFather(position);			
	s_idx = Trial[position];
	f_idx = Trial[Tree_GetFather(position)];
  }
};

//================================================================
void Tree_PullFromTree(int point)
//================================================================
{
  float Uv = U[point];
  U[point] = 0;
  Tree_UpdateChange(Tree[point]);
  U[Tree_PopHead()] = Uv;
};

//================================================================
int Tree_GetSize()
//================================================================
{
  return PtrToTrial - Trial + 1;
};

//================================================================
//================================================================
//================================================================
Edge **_C = NULL;  // constraint map
//================================================================
//================================================================
void InitializeNeighborhoods()
//================================================================
{
    _labels = 0;
	connectivity_large = 8;
    // Freeman code of the 8-adjacency pixel graph
	NeighborhoodLarge = new int[connectivity_large];
    v8c = new int[connectivity_large+2];
	v8c[8] = v8c[0] = NeighborhoodLarge[0]= -Nx-1;
	v8c[9] = v8c[1] = NeighborhoodLarge[1]= -Nx;
	v8c[2] = NeighborhoodLarge[2]= -Nx+1;
	v8c[7] = NeighborhoodLarge[3]=  -1; //1;
    v8c[3] = NeighborhoodLarge[4]=  1; //Nx+1;
    v8c[6] = NeighborhoodLarge[5]=  Nx-1;//Nx;
    v8c[5] = NeighborhoodLarge[6]=  Nx; //Nx-1;
    v8c[4] = NeighborhoodLarge[7]=  Nx+1;//-1;
    //-----------------------------------------------------------
    Simplicies1 = new int   [connectivity_large];
    Simplicies2 = new int   [connectivity_large];
    signsX      = new int   [connectivity_large];
    signsY      = new int   [connectivity_large];
    h1          = new float [connectivity_large];
    h2          = new float [connectivity_large];
    h1_h2       = new float [connectivity_large];
    h22         = new float [connectivity_large];
    //-----------------------------------------------------------
    Simplicies1[0] = 1; Simplicies2[0] = 1+Nx;
    signsX[0] = -1; signsY[0] = -1;
    h1[0] = hy; h2[0] = hx; h1_h2[0] = hy/hx; h22[0] = hx*hx;
    Simplicies1[1] = 1; Simplicies2[1] =  1-Nx;
    signsX[1] = -1; signsY[1] = 1;
    h1[1] = hy; h2[1] = hx; h1_h2[1] = hy/hx; h22[1] = hx*hx;
    Simplicies1[2] = -1; Simplicies2[2] = -1+Nx;
    signsX[2] = 1; signsY[2] = -1;
    h1[2] = hy; h2[2] = hx; h1_h2[2] = hy/hx; h22[2] = hx*hx;
    Simplicies1[3] = -1; Simplicies2[3] = -1-Nx;
    signsX[3] = 1; signsY[3] = 1;
    h1[3] = hy; h2[3] = hx; h1_h2[3] = hy/hx; h22[3] = hx*hx;
    Simplicies1[4] = Nx; Simplicies2[4] =  1+Nx;
    signsX[4] = -1; signsY[4] = -1;
    h1[4] = hx; h2[4] = hy; h1_h2[4] = hx/hy; h22[4] = hy*hy;
    Simplicies1[5] = Nx; Simplicies2[5] = -1+Nx;
    signsX[5] = 1; signsY[5] = -1;
    h1[5] = hx; h2[5] = hy; h1_h2[5] = hx/hy; h22[5] = hy*hy;
    Simplicies1[6] = -Nx; Simplicies2[6] = 1-Nx;
    signsX[6] = -1; signsY[6] = 1;
    h1[6] = hx; h2[6] = hy; h1_h2[6] = hx/hy; h22[6] = hy*hy;
    Simplicies1[7] = -Nx; Simplicies2[7] = -1-Nx;
    signsX[7] = 1; signsY[7] = 1;
    h1[7] = hx; h2[7] = hy; h1_h2[7] = hx/hy; h22[7] = hy*hy;
    kDead = -1;
    kOpen = -2;
    kFar = -3;
    kBorder = -4;
};

//================================================================
void InitializeArrays()
//================================================================
{
  int x, y, point;
  Seeds = new std::map<float,int>;
  _edges_to_split = new std::queue<Edge*>;
  _inverted_triangles = new VoronoiQueue;
  _triangles = new std::list<Triangle>;
  //copy the weight list and initialize
  S = new short[size];
  Tree = new int[size];
  Trial = new int[size];
  M1 = new float[size];
  M2 = new float[size];
  M3 = new float[size];
  q_gradient = new float[3];
  _C = new Edge*[size];
  //------------------------------------------------------------
  for (x = 0; x < nx; x++)
  {
    for (y = 0; y < ny; y++)
    {
      point = (x+1) + (y+1)*Nx;
      M1[point] = T[x + y*nx];
      M2[point] = T[x + y*nx + 3*nxny];
      M3[point] = T[x + y*nx + nxny];
      V[point] = -1; S[point] = kFar;
      if ((hx != 1.0) || (hy != 1.0))
      {
        M1[point] = hx2* T[x + y*nx];
        M2[point] = hy2* T[x + y*nx + 3*nxny];
        M3[point] = hxhy*T[x + y*nx +   nxny];
      }
    }
  }
  for (x = 0; x < size; x++)
  {
    _C[x] = NULL;  
    U[x] = INFINITE;
    Tree[x] = -1;
  }
  //------------------------------------------------------------
  PtrToTrial = Trial - 1;
  //------------------------------------------------------------
  // Initialize Borders
  for (x = 0; x < Nx; x++) {
    y = 0;
    point = x + y*Nx;
    V[point] = kBorder; S[point] = kBorder;
    y = Ny-1;
    point = x + y*Nx;
    V[point] = kBorder; S[point] = kBorder;
  }
  for (y = 0; y < Ny; y++) {
    x = 0;
    point = x + y*Nx;
    V[point] = kBorder; S[point] = kBorder;
    x = Nx-1;
    point = x + y*Nx;
    V[point] = kBorder; S[point] = kBorder;
  }
};

//================================================================
bool insertSeed(int p, float l)
//================================================================
{
    if (p >= size) return false;
    if (S[p] == kSeed) return false;
    U[p] = 0.0;
    V[p] = l;
    S[p] = kSeed;
    Tree_PushIn(p);
    Seeds->insert(std::make_pair(l,p));
}

//================================================================
void insertSeed(int x, int y, float l)
//================================================================
{
    x += 1; y += 1;
    insertSeed(x+y*Nx, l);
}

//================================================================
void addGridBoundaryAsConstraints()
//================================================================
{
  Edge *e = NULL;
  int i, end = 2*(Nx-1);
  //--------------------------------------------------------------
  // inserting grid boundary as edges
  e = new Edge(_labels, _labels+1);
  //Edges.insert(std::make_pair(e, 0.));
  for (i = Nx+2; i < end; i++)
    {
      _C[i] = e;
      e->points.push(i);
    }
  //--------------------------------------------------------------
  // 
  e = new Edge(_labels, _labels+2);
  //Edges.insert(std::make_pair(e, 0.));
  end = (Ny-3)*Nx+1;
  for (i = 2*Nx+1; i <= end; i += Nx)
    {
      _C[i] = e;
      e->points.push(i);
    }
  //--------------------------------------------------------------
  e = new Edge(_labels+1, _labels+3);
  //Edges.insert(std::make_pair(e, 0.));
  end = (Ny-2)*Nx-1;
  for (i = 3*Nx-2; i < end; i += Nx)
    {
      _C[i] = e;
      e->points.push(i);
    }
  //--------------------------------------------------------------
  e = new Edge(_labels+2, _labels+3);
  //Edges.insert(std::make_pair(e, 0.));
  end = (Ny-1)*Nx-2;
  for (i = (Ny-2)*Nx+2; i < end; i++)
    {
      _C[i] = e;
      e->points.push(i);
    }
  //--------------------------------------------------------------
  insertSeed(0, 0, _labels);
  _labels+=1;
  insertSeed(nx-1, 0, _labels);
  _labels+=1;
  insertSeed(0, ny-1, _labels);
  _labels+=1;
  insertSeed(nx-1, ny-1, _labels);
  _labels+=1;
}

//================================================================
void InitializeOpenHeap()
//================================================================
{
  int point, i, j;
  
  for (int s = 0; s < nb_start_points; ++s)
    {
      i = (int)round(start_points[2*s]);
      j = (int)round(start_points[1+2*s]);
      point = i + j*Nx;
      //--------------------------------------------------------
      if (point >= size) mexErrMsgTxt("start_points should be in the domain.");
      //--------------------------------------------------------
      if (U[point] == 0) mexErrMsgTxt("start_points should not contain duplicates.");
      //--------------------------------------------------------
      U[point] = 0.0; S[point] = kSeed; V[point] = _labels;
      // add to heap
      Tree_PushIn(point);
      Seeds->insert(std::make_pair(_labels,point));
      _labels+=1;
      //--------------------------------------------------------
    }
};

//================================================================
void TsitsiklisQuadrantLength(float *result, float Pc,float Ua,float Ub, int QuadNb)
//================================================================
{
    float ha, hb, hb2, ha_hb;
    ha = h1[QuadNb];
    hb = h2[QuadNb];
    hb2= h22[QuadNb];
    ha_hb = h1_h2[QuadNb];
    //-----------------------------------------------------------
    if (Ua <= Ub) *result = Ua + Pc*ha;
	//-----------------------------------------------------------
	else
        if (Pc*hb <= (Ua-Ub)*sqrt_hx2_plus_hy2/hb)
            *result = Ub + Pc*sqrt_hx2_plus_hy2;
        //-----------------------------------------------------------
        else *result = Ua + ha_hb*sqrt(Pc*Pc*hb2 - (Ua-Ub)*(Ua-Ub));
};

//================================================================
float TsitsiklisQuadrantGradient(float* result_gradient, float m1, float m2, float m3,
                                 float Ua, float Ub, const float &Va, const float &Vb, int QuadNb)
//================================================================
{
	float ha, hb, hb2, ha_hb;
    float k1, k2, alpha, r11, r12, r22, R;
    //-----------------------------------------------------------
    k1 = Ua - Ub;
    k2 = Ub;
    if(QuadNb < 4) {
        r11 = m2;
        r12 = m2 + m3*(signsX[QuadNb]*signsY[QuadNb]);
    }
    else {
        r11 = m1;
        r12 = m1 + m3*(signsX[QuadNb]*signsY[QuadNb]);
    }
    r22 = m2 + 2.0*m3*(signsX[QuadNb]*signsY[QuadNb]) + m1;
    R = r11*r22 - r12*r12; // toujours positif
    //-----------------------------------------------------------
    if (k1 >= sqrt(r11)) {
		*result_gradient = Ub + sqrt(r22);
        //result_gradient[1] = signsX[QuadNb]/sqrt(r22);
        //result_gradient[2] = signsY[QuadNb]/sqrt(r22);
        return Vb;
	}
    //-----------------------------------------------------------
    else
        if(k1 <= -sqrt(r11)) {
            if(QuadNb < 4) {
                *result_gradient = Ua + sqrt(m1);
                //result_gradient[1] = signsX[QuadNb]/sqrt(m1);
                //result_gradient[2] = 0.0;
            }
            else {
                *result_gradient = Ua + sqrt(m2);
                //result_gradient[1] = 0.0;
                //result_gradient[2] = signsY[QuadNb]/sqrt(m2);
            }
            return Va;
        }
        //-----------------------------------------------------------
        else{
            //-------------------------------------------------------
            if (r12 <= k1*sqrt(R/(r11-k1*k1))) {
                *result_gradient = Ub + sqrt(r22);
                //result_gradient[1] = signsX[QuadNb]/sqrt(r22);
                //result_gradient[2] = signsY[QuadNb]/sqrt(r22);
                return Vb;
            }
            //-------------------------------------------------------
            else if (r12 > (r11 + k1*sqrt(R/(r11-k1*k1)))) {
                if (QuadNb < 4) {
                    *result_gradient = Ua + sqrt(m1);
                    //result_gradient[1] = signsX[QuadNb]/sqrt(m1);
                    //result_gradient[2] = 0.0;
                }
                else {
                    *result_gradient = Ua + sqrt(m2);
                    //result_gradient[1] = 0.0;
                    //result_gradient[2] = signsY[QuadNb]/sqrt(m2);
                }
                return Va;
            }
            //-------------------------------------------------------
            else {
                alpha = (r12 - k1 * sqrt(R/(r11-k1*k1))) / r11;
                *result_gradient = alpha*k1 + k2 + sqrt(R/(r11-k1*k1));
                /*if (QuadNb < 4) {
                    result_gradient[1] = signsX[QuadNb] / sqrt(R/(r11-k1*k1));
                    result_gradient[2] = (1-alpha) * signsY[QuadNb] / sqrt(R/(r11-k1*k1));
                }
                else {
                    result_gradient[1] = (1-alpha) * signsX[QuadNb] / sqrt(R/(r11-k1*k1));
                    result_gradient[2] = signsY[QuadNb] / sqrt(R/(r11-k1*k1));
                }*/
                return ((Ua+sqrt(r11)) < (Ub+sqrt(r22)) ? Va : Vb);
                //return (Ua < Ub ? Va : Vb);
            }
        }
    return -1;
};

//================================================================
bool TsitsiklisUpdate(int point)
/*
COMMENTS : 
*/
//================================================================
{
	int npoint1, npoint2;
    float Ua, Ub, La, Lb;
	float Ur = U[point], Utmp;
	float Vr, Vtmp;
	//float Lr = L[point];
	float dUxr = dUx[point];
	float dUyr = dUy[point];
	bool is_updated = false;
	//--------------------------------------------------------------
	// Get the U & L values for each neighbor.
	for (int i = 0; i < connectivity_large; i++)
    {
		npoint1 = point + Simplicies1[i];
        npoint2 = point + Simplicies2[i];
		if (S[npoint1] == kDead || S[npoint1] == kSeed) { Ua = U[npoint1]; /*La = L[npoint1];*/ }
        else { Ua = INFINITE; La = INFINITE; }
        if (S[npoint2] == kDead || S[npoint1] == kSeed) { Ub = U[npoint2]; /*Lb = L[npoint2];*/ }
        else { Ub = INFINITE; Lb = INFINITE; }
        Vtmp = TsitsiklisQuadrantGradient(&Utmp, M1[point], M2[point], M3[point],
                                          Ua, Ub, V[npoint1], V[npoint2], i);
        if (Utmp < Ur)
        {
            Ur   =  Utmp; //q_gradient[0];
    		//dUxr =  q_gradient[1]; dUyr =  q_gradient[2];
            //TsitsiklisQuadrantLength(&Lr, 1, La, Lb, i);
            Vr = Vtmp;
            is_updated = true;
    	}
	}
	//--------------------------------------------------------------
	if (is_updated)
    {
		U[point] = Ur;
		//L[point] = Lr;
		V[point] = Vr;
		//dUx[point] = dUxr; dUy[point] = dUyr;
	}
	//--------------------------------------------------------------
	return is_updated;
}

//================================================================
void CorrectMaps()
//================================================================
{
	int point;
	for (point = 0; point < size; point++) {
		if (V[point] < 0) { U[point] = 0; /*L[point] = 0;*/ }
    }
}

//================================================================
int NextPoint()
//================================================================
{
  float Umax = 0.0;
  int fps = -1;
  for (int x = Nx+1; x < size;  x++)
  {
      if (S[x] != kBorder)
      {
          if (S[x] != kSeed)
          {
                if (Umax < U[x]) { Umax = U[x]; fps = x; }
                S[x] = kFar;
          }
          Tree[x] = -1;
      }
  }
  PtrToTrial = Trial - 1;
  return fps;
}

//================================================================
struct Backup
{
    int point;
    float u;
    float v;
    Backup(const int &p, const float &uu, const float &vv)
    : point(p), u(uu), v(vv)
    {
    
    }
    ~Backup() { }
};

//================================================================
bool RunLocalPropagation(const int &seed)   // for seed not on a constrained edge
//================================================================
{
  std::list<Backup> backup_list;
  backup_list.push_back(Backup(seed, U[seed], V[seed]));
  U[seed] = 0.0;
  S[seed] = kSeed;
  V[seed] = _labels;
  Tree_PushIn(seed);
  //--------------------------------------------------------------
  int point, npoint, k;
  float v;
  bool is_updated = false;
  float u;
  Edge *e = NULL;
  //--------------------------------------------------------------
  while (Tree_GetSize() > 0)
    {
      point = Tree_PopHead();
      //--------------------------------------------------------------
      if ((e = _C[point]) != NULL)  // front touch a constrained edge
      {
          std::list<Backup>::iterator it = backup_list.begin(), end = backup_list.end();
          for (; it != end; ++it)
          {
              U[it->point] = it->u;
              V[it->point] = it->v;
              S[it->point] = kFar;
              Tree[it->point] = -1;
          }
          PtrToTrial = Trial - 1;
          if (e->midpoint == 0) e->midpoint = point;
          _edges_to_split->push(e);
          //std::cerr << "runnot " << e->label1 << " " << e->label2 << " seeds=";
          //print_point((*Seeds)[e->label1]);
          //std::cerr << " ";
          //print_point((*Seeds)[e->label2]);
          //std::cerr  << " midpoint=";
          //print_point(e->midpoint);
          //std::cerr << " point=";
          //print_point(point);
          //std::cerr << std::endl << std::flush;
          return false;
      }
      //--------------------------------------------------------------
      //if (_umax > U[point]) std::cerr << _umax << " non-monotone" << U[point] << std::endl;
      //else _umax = U[point];
      //--------------------------------------------------------------
      S[point] = (S[point] == kSeed ? kSeed : kDead);
      //Lmax = MAX(Lmax, L[point]);
      //if ((Dmax > 0) && (Lmax > Dmax)) break;
      //--------------------------------------------------------------
      for (k = 0; k < connectivity_large; k++)
	  {
        npoint = point + NeighborhoodLarge[k];
        //--------------------------------------------------------------
        
        //--------------------------------------------------------------
        if (S[npoint] == kOpen)
	    {
	      is_updated = TsitsiklisUpdate(npoint);
	      if (is_updated) Tree_UpdateChange(Tree[npoint]);
          continue;
	    }
        //--------------------------------------------------------------
        if (S[npoint] == kFar)
	    {
            u = U[npoint];
            v = V[npoint];
	        is_updated = TsitsiklisUpdate(npoint);
            if (is_updated)
            {
              S[npoint] = kOpen;
              Tree_PushIn(npoint);
              backup_list.push_back(Backup(npoint,u,v));
            }
        }
      }
  }
  Seeds->insert(std::make_pair(_labels,seed));
  _labels+=1;
  return true;
}

//================================================================
int split(Edge *e)
//================================================================
{
    //--------------------------------------------------------------
    if ((e->points.size() < 1) || (e->midpoint == 0))
    {
        //std::cerr << "ERROR : edge cannot be splitted or non-midpoint\n";
        return -1;
    }
    //--------------------------------------------------------------
    int seed = e->midpoint, p;
    Edge *e1 = new Edge(e->label1, _labels);
    Edge *e2 = new Edge(_labels, e->label2);
    //--------------------------------------------------------------
    // first subedge
    while (!e->points.empty())
    {
        
       p = e->points.front();
       e->points.pop();
	   if (p == seed) // go to next subedge
	   {
           _C[p] = NULL;
           break;
       }
       _C[p] = e1;
       e1->points.push(p);
    }
    //--------------------------------------------------------------
    // second subedge
    while (!e->points.empty())
    {
       p = e->points.front();
       _C[p] = e2;
       e2->points.push(p);
       e->points.pop();
    }
    // delete e;
    return seed;
}


//================================================================
bool insertOnEdge(Edge *e)   // for a seed on a constrained edge
//================================================================
{
  // split the edge in two subedges
  int seed = split(e);
  if (seed == -1)
  {
      //std::cerr << "ERROR" << std::endl;
      return false;
  }
  U[seed] = 0.0;
  S[seed] = kSeed;
  V[seed] = _labels;
  Tree_PushIn(seed);
  //--------------------------------------------------------------
  int point, npoint, k;
  float v;
  bool is_updated = false;
  Edge *ed = NULL;
  //--------------------------------------------------------------
  while (Tree_GetSize() > 0)
    {
      point = Tree_PopHead();
      //--------------------------------------------------------------
      if ((ed = _C[point]) != NULL)  // front touch a constrained edge
      {
          if ((ed->mark == false) && (ed->label1 != _labels) && (ed->label2 != _labels))
          {
              //std::cerr << "split push " << ed->label1 << " " << ed->label2 << " ";
              //print_point(ed->midpoint);
              //std::cerr << " ";
              //print_point(point);
              //std::cerr << " at ";
              //print_point(seed);
              //std::cerr << std::endl;
              _edges_to_split->push(ed);
              ed->mark = true;
          }
      }
      //--------------------------------------------------------------
      S[point] = (S[point] == kSeed ? kSeed : kDead);
      v = V[point];
      //--------------------------------------------------------------
      for (k = 0; k < connectivity_large; k++)
	  {
        npoint = point + NeighborhoodLarge[k];
        //--------------------------------------------------------------
        if (S[npoint] == kOpen)
	    {
	      is_updated = TsitsiklisUpdate(npoint);
	      if (is_updated) Tree_UpdateChange(Tree[npoint]);
          continue;
	    }
        //--------------------------------------------------------------
        if (S[npoint] == kFar)
	    {
            is_updated = TsitsiklisUpdate(npoint);
            if (is_updated) // update
            {
              S[npoint] = kOpen;
              Tree_PushIn(npoint);
            }
            else  // try to find a constrained edge midpoint
            {
                if ((ed != NULL) && (ed->midpoint == 0) && (_C[npoint] != NULL) && (V[npoint] != v))
                {
                    ed->midpoint = point;
                    //std::cerr << "midpoint detected ";
                    //print_point(ed->midpoint);
                    //std::cerr << std::endl;
                    
                }
            }
            
        }
      }
  }
  Seeds->insert(std::make_pair(_labels,seed));
  _labels+=1;
  return true;
}

//================================================================
int orient2D(int ax, int ay, int bx, int by, int cx, int cy)
//================================================================
{
    int d = (ax*(by-cy) - ay*(bx-cx) + (bx*cy - cx*by));
    return (d > 0 ? 1 : (d < 0 ? -1 : 0));
}

//================================================================
int orient2D(int p1, int p2, int p3)
//================================================================
{
    int y1 = p1 / Nx, y2 = p2 / Nx, y3 = p3 / Nx;
    return orient2D(p1-Nx*y1, y1, p2-Nx*y2, y2, p3-Nx*y3, y3);
}

//================================================================
bool isInvertedTriangle(int p1, int p2, int p3,
                        int s1, int s2, int s3)
//================================================================
{
    return (orient2D(p1, p2, p3) != orient2D(s1, s2, s3));
}

//================================================================
bool isInvertedVoronoiVertex(const int &p1, const int &p2, const int &p3,
                             const float &v1, const float &v2, const float &v3,
                             const float &vd)
//================================================================
{
    if (vd == v2) return true;
    if ((vd == v1) || (vd == v3))
    {
        if (isInvertedTriangle(p1, p2, p3, (*Seeds)[v1], (*Seeds)[v2], (*Seeds)[v3]))
        {
            float u1, u2, u3;
            u1 = U[p1];
            u2 = U[p2];
            u3 = U[p3];
            std::cerr << "\n";
            std::cerr << "yo1 " << u1 << " " << u2 << " " << u3 << " ,";
            print_point(p1);
            std::cerr << " ";
            print_point(p2);
            std::cerr << " ";
            print_point(p3);
            std::cerr << ", ";
            print_point((*Seeds)[v1]);
            std::cerr << " ";
            print_point((*Seeds)[v2]);
            std::cerr << " ";
            print_point((*Seeds)[v3]);
            std::cerr << std::endl;
            
            if (u1 > 0.0 && u1 >= u2 && u1 >= u3)
            {
                _inverted_triangles->push(VoronoiVertex(p1, v1, u1));
                return true;
            }
            if (u2 > 0.0 && u2 >= u1 && u2 >= u3)
            {
                _inverted_triangles->push(VoronoiVertex(p2, v2, u2));
                return true;
            }
            if (u3 > 0.0 && u3 >= u1 && u3 >= u2)
            {
                _inverted_triangles->push(VoronoiVertex(p3, v3, u3));
                return true;
            }
            std::cerr << "shot1 " << u1 << " " << u2 << " " << u3 << std::endl;
            int y = p1 / Nx;
            int x = p1 - Nx * y;
            std::cerr << "[" << x << ";" << y << "] ";
            y = p2 / Nx;
            x = p2 - Nx * y;
            std::cerr << "[" << x << ";" << y << "] ";
            y = p3 / Nx;
            x = p3 - Nx * y;
            std::cerr << "[" << x << ";" << y << "] " << std::endl << std::flush;
        }
        //--------------------------------------------------------------
        // all three points are seeds, insert triangle
        _triangles->push_back(Triangle(v1, v2, v3, p1));
        return true;
    }
    return false;
}

//================================================================
bool addVoronoiInverted(const int &p1, const int &p2, const int &p3,
                        const float &v1, const float &v2, const float &v3)
//================================================================
{
    float u1 = U[p1], u2 = U[p2], u3 = U[p3];
    int o1 = orient2D(p1,p2,p3);
    int o2 = orient2D((*Seeds)[v1],(*Seeds)[v2],(*Seeds)[v3]);
    //--------------------------------------------------------------
    if (o1 == o2)
    {
        _triangles->push_back(Triangle(v1, v2, v3, p1));
        return false;
    }
    //--------------------------------------------------------------
    if ((S[p1] != kSeed) && (u1 >= u2) && u1 >= u3) 
    {
        _inverted_triangles->push(VoronoiVertex(p1, v1, u1));
        return true;
    }
    //--------------------------------------------------------------
    if ((S[p2] != kSeed) && (u2 >= u1) && (u2 >= u3))
    {
        _inverted_triangles->push(VoronoiVertex(p2, v2, u2));
        return true;
    }
    //--------------------------------------------------------------
    if ((S[p3] != kSeed) && (u3 >= u1) && (u3 >= u2))
    {
        _inverted_triangles->push(VoronoiVertex(p3, v3, u3)); 
        return true;
    }
    //--------------------------------------------------------------
    //std::cerr << " error ";
    //print_point(p1);std::cerr << " ";
    //print_point(p2);std::cerr << " ";
    //print_point(p3);std::cerr << " ; ";
    //print_point((*Seeds)[v1]);std::cerr << " ";
    //print_point((*Seeds)[v2]);std::cerr << " ";
    //print_point((*Seeds)[v3]);
    //std::cerr << " " << orient2D(p1,p2,p3) << " " << orient2D((*Seeds)[v1],(*Seeds)[v2],(*Seeds)[v3]);
    //std::cerr << std::endl;
    //--------------------------------------------------------------
    return false;
}

//================================================================
bool addQuadVoronoi(const int &p, const int &p1, const int &p2, const int &p3,
                    const float &v, const float &v1, const float &v2, const float &v3)
//================================================================
{
    //--------------------------------------------------------------
    if ((S[p] == kSeed) && (S[p1] == kSeed) &&
        (S[p2] == kSeed) && (S[p3] == kSeed))
    {
        _triangles->push_back(Triangle(v, v1, v2, p));
        _triangles->push_back(Triangle(v2, v3, v, p2));
        return false;
    }
    //--------------------------------------------------------------
    float u = U[p], u1 = U[p1], u2 = U[p2], u3 = U[p3];
    float uu2 = u - u2, u1u3 = u1 - u3;
    if (uu2 < 0.0) uu2 *= -1.0;
    if (u1u3 < 0.0) u1u3 *= -1.0;
    int o1, o2, o3, o4;
    //--------------------------------------------------------------
    if (uu2 < u1u3)
    {
        o1 = orient2D(p,p1,p2);
        o2 = orient2D((*Seeds)[v],(*Seeds)[v1],(*Seeds)[v2]);
        o3 = orient2D(p2,p3,p);
        o4 = orient2D((*Seeds)[v2],(*Seeds)[v3],(*Seeds)[v]);
        //--------------------------------------------------------------
        if ((o1 == o2) && (o3 == o4))
        {
            _triangles->push_back(Triangle(v, v1, v2, p));
            _triangles->push_back(Triangle(v2, v3, v, p2));
            return false;
        }
        //--------------------------------------------------------------
        if (o1 != o2)
        {
            if (S[p1] != kSeed)
            {
                _inverted_triangles->push(VoronoiVertex(p1, v1, u1));
                return true;
            }
            if (S[p2] != kSeed)
            {
                _inverted_triangles->push(VoronoiVertex(p2, v2, u2));
                return true;
            }
            if (S[p] != kSeed)
            {
                _inverted_triangles->push(VoronoiVertex(p, v, u));
                return true;
            }
        }
        //--------------------------------------------------------------
        if (o3 != o4)
        {
            if (S[p3] != kSeed)
            {
                _inverted_triangles->push(VoronoiVertex(p3, v3, u3));
                return true;
            }
            if (S[p2] != kSeed)
            {
                _inverted_triangles->push(VoronoiVertex(p2, v2, u2));
                return true;
            }
            if (S[p] != kSeed)
            {
                _inverted_triangles->push(VoronoiVertex(p, v, u));
                return true;
            }
        }
        //std::cerr << "error\n";
        return false;
    }
    //--------------------------------------------------------------
    o1 = orient2D(p1,p2,p3);
    o2 = orient2D((*Seeds)[v1],(*Seeds)[v2],(*Seeds)[v3]);
    o3 = orient2D(p3,p,p1);
    o4 = orient2D((*Seeds)[v3],(*Seeds)[v],(*Seeds)[v1]);
    //--------------------------------------------------------------
    if ((o1 == o2) && (o3 == o4))
    {
        _triangles->push_back(Triangle(v1, v2, v3, p1));
        _triangles->push_back(Triangle(v3, v, v1, p3));
        return false;
    }
    //--------------------------------------------------------------
    if (o1 != o2)
    {
        if (S[p2] != kSeed)
        {
            _inverted_triangles->push(VoronoiVertex(p2, v2, u2));
            return true;
        }
        if (S[p1] != kSeed)
        {
            _inverted_triangles->push(VoronoiVertex(p1, v1, u1));
            return true;
        }
        if (S[p3] != kSeed)
        {
            _inverted_triangles->push(VoronoiVertex(p3, v3, u3));
            return true;
        }
    }
    //--------------------------------------------------------------
    if (o3 != o4)
    {
        if (S[p] != kSeed)
        {
            _inverted_triangles->push(VoronoiVertex(p, v, u));
            return true;
        }
        if (S[p3] != kSeed)
        {
            _inverted_triangles->push(VoronoiVertex(p3, v3, u3));
            return true;
        }
        if (S[p1] != kSeed)
        {
            _inverted_triangles->push(VoronoiVertex(p1, v1, u1));
            return true;
        }
    }
    
    //std::cerr << " error quad ";
    /*print_point(p1);std::cerr << " ";
    print_point(p2);std::cerr << " ";
    print_point(p3);std::cerr << " ; ";
    print_point((*Seeds)[v1]);std::cerr << " ";
    print_point((*Seeds)[v2]);std::cerr << " ";
    print_point((*Seeds)[v3]);
    std::cerr << " " << orient2D(p1,p2,p3) << " " << orient2D((*Seeds)[v1],(*Seeds)[v2],(*Seeds)[v3]);
    std::cerr << std::endl;*/
    return false;
}

//================================================================
bool isInvertedVoronoiVertex(const int &point)
//================================================================
{
    float v = V[point];
    int npoint1 = point + v8c[3];
    float vnp1 = V[npoint1];
    int npoint2 = point + v8c[4];
    float vnp2 = V[npoint2];
    int npoint3 = point + v8c[5];
    float vnp3 = V[npoint3];
    bool quad_inv = false;
    //--------------------------------------------------------------
    bool a = (vnp1 != vnp2), b = (vnp1 != vnp3), c = (vnp2 != vnp3), d = (v != vnp1),
        e = (v != vnp2), f = (v != vnp3);
    //--------------------------------------------------------------
    bool i1 = isInvertedTriangle(npoint1, npoint2, npoint3, (*Seeds)[vnp1], (*Seeds)[vnp2], (*Seeds)[vnp3]);
    bool i2 = isInvertedTriangle(npoint2, npoint3, point, (*Seeds)[vnp2], (*Seeds)[vnp3], (*Seeds)[v]);
    bool i3 = isInvertedTriangle(npoint3, point, npoint1, (*Seeds)[vnp3], (*Seeds)[v], (*Seeds)[vnp1]);
    bool i4 = isInvertedTriangle(point, npoint1, npoint2, (*Seeds)[v], (*Seeds)[vnp1], (*Seeds)[vnp2]);
    //--------------------------------------------------------------
    // quad detection
    if (a && b && c && d && e && f)
    {
        return addQuadVoronoi(point, npoint1, npoint2, npoint3, v, vnp1, vnp2, vnp3);
    }
    //--------------------------------------------------------------
    if (a && b && c)
    {
        if (!e) return false;
        return addVoronoiInverted(npoint1, npoint2, npoint3, vnp1, vnp2, vnp3);
    }
    //--------------------------------------------------------------
    if (c && f && e)
    {
        if (!b) return false;
        return addVoronoiInverted(npoint2, npoint3, point, vnp2, vnp3, v);
    }
    //--------------------------------------------------------------
    if (f && d && b)
    {
        if (!e) return false;
        return addVoronoiInverted(npoint3, point, npoint1, vnp3, v, vnp1);
    }
    //--------------------------------------------------------------
    if (d && a && e)
    {
        if (!b) return false;
        return addVoronoiInverted(point, npoint1, npoint2, v, vnp1, vnp2);
    }
    return false;
}

//================================================================
bool insertVoronoiInverted(const int &seed)   // for VoronoiVertex as seed when dual triangle is inverted
//================================================================
{
  std::list<Backup> backup_list;
  backup_list.push_back(Backup(seed, U[seed], V[seed]));
  //std::cerr << "inverted insertion ";
  //print_point(seed);
  //std::cerr << S[seed] << std::endl << std::flush;
  U[seed] = 0.0;
  S[seed] = kSeed;
  V[seed] = _labels;
  Tree_PushIn(seed);
  //--------------------------------------------------------------
  int point, npoint, k;
  float v;
  bool is_updated = false;
  float u;
  Edge *e = NULL;
  //--------------------------------------------------------------
  while (Tree_GetSize() > 0)
    {
      point = Tree_PopHead();
      //--------------------------------------------------------------
      if ((e = _C[point]) != NULL)  // front touch a constrained edge
      {
          std::list<Backup>::iterator it = backup_list.begin(), end = backup_list.end();
          for (; it != end; ++it)
          {
              U[it->point] = it->u;
              V[it->point] = it->v;
              S[it->point] = kFar;
              Tree[it->point] = -1;
          }
          PtrToTrial = Trial - 1;
          if (e->midpoint == 0) e->midpoint = point;
          _edges_to_split->push(e);
          return false;
      }
      //--------------------------------------------------------------
      S[point] = (S[point] == kSeed ? kSeed : kDead);
      //--------------------------------------------------------------
      for (k = 0; k < connectivity_large; k++)
	  {
        npoint = point + NeighborhoodLarge[k];
        //--------------------------------------------------------------
        
        //--------------------------------------------------------------
        if (S[npoint] == kOpen)
	    {
	      is_updated = TsitsiklisUpdate(npoint);
	      if (is_updated) Tree_UpdateChange(Tree[npoint]);
          continue;
	    }
        //--------------------------------------------------------------
        if (S[npoint] == kFar)
	    {
            u = U[npoint];
            v = V[npoint];
	        is_updated = TsitsiklisUpdate(npoint);
            if (is_updated)
            {
              S[npoint] = kOpen;
              Tree_PushIn(npoint);
              backup_list.push_back(Backup(npoint,u,v));
            }
        }
      }
  }
  //--------------------------------------------------------------
  Seeds->insert(std::make_pair(_labels,seed));
  _labels+=1;
  //--------------------------------------------------------------
  // prepare next insertion
  std::list<Backup>::iterator it = backup_list.begin(), end = backup_list.end();
  for (; it != end; ++it)
  {
      if (S[it->point] != kSeed) S[it->point] = kFar;   // ??? to optimize
      Tree[it->point] = -1;
  }
  PtrToTrial = Trial - 1;
  //--------------------------------------------------------------
  return true;
}

//================================================================
bool validateTriangulation()
//================================================================
{
    //std::cerr << "validate " << _labels << std::endl << std::flush;
    int x, y, point, endx = Nx - 2, endy = Ny - 2;
    Edge *e = NULL;
    //--------------------------------------------------------------
    // find Voronoi vertices and detect inverted triangles
    for (x = 1; x < endx; x++)
    {
        for (y = 1; y < endy; y++)
        {
            point = x + y*Nx;
            if ((S[point] == kBorder))// || (U[point] == 0.0))
                continue;
            isInvertedVoronoiVertex(point);
        }
    }
    //std::cerr << _labels << " " << _inverted_triangles->size() << std::endl;
    
    //return true;
    //--------------------------------------------------------------
    // no inverted triangles : triangulation is valid
    if (_inverted_triangles->empty()) return true;
    //--------------------------------------------------------------
    // insert Voronoi vertices to break inverted triangles
    while (!_inverted_triangles->empty())
    {
        const VoronoiVertex &vv = _inverted_triangles->top();
        //--------------------------------------------------------------
        // check validity
        if (vv.label != V[vv.point]) { _inverted_triangles->pop(); continue; }
        //y = vv.point / Nx;
        //x = vv.point - Nx * y;
        //std::cerr << "[" << x << ";" << y << "],";
        //--------------------------------------------------------------
        // try to insert
        if (insertVoronoiInverted(vv.point) == true)
        {
            //std::cerr << "inserted" << std::endl << std::flush;
            _inverted_triangles->pop();
            continue;
        }
        //--------------------------------------------------------------
        // encroached edge detected, insert edge midpoint
        while (!_edges_to_split->empty())
        {
            e = _edges_to_split->front();
            _edges_to_split->pop();
            insertOnEdge(e);
            //--------------------------------------------------------------
            // prepare next insertion
            for (point = 0; point < size;  point++)
            {
                if (S[point] != kBorder)
                {
                    if (S[point] != kSeed) S[point] = kFar;
                    Tree[point] = -1;
                }
                
            }
            PtrToTrial = Trial - 1;
        }
        //--------------------------------------------------------------
        // check if inverted triangle is still alive
        //if (vv.label == V[vv.point]) insertVoronoiInverted(vv.point);
        //--------------------------------------------------------------
        _inverted_triangles->pop();
    }
    return false;
}

//================================================================
void RunPropagation(int nb_inserted_points, float Ustop = 0.0)
//================================================================
{
  if (nb_inserted_points > nx*ny) nb_inserted_points = nx*ny;
  int point, npoint, k, pmax;
  float v, vnp;
  bool is_updated = false;
  Edge *e = NULL;
  //--------------------------------------------------------------
  while (Tree_GetSize() > 0)
    {
        
      point = Tree_PopHead();
      //--------------------------------------------------------------
      S[point] = (S[point] == kSeed ? kSeed : kDead);
      e = _C[point];
      v = V[point];
      //--------------------------------------------------------------
      for (k = 0; k < connectivity_large; k++)
	  {
        npoint = point + NeighborhoodLarge[k];
        //--------------------------------------------------------------
        if (S[npoint] == kOpen)
	    {
	      is_updated = TsitsiklisUpdate(npoint);
	      if (is_updated) Tree_UpdateChange(Tree[npoint]);
          continue;
	    }
        //--------------------------------------------------------------
        if (S[npoint] == kFar)
	    {
	      S[npoint] = kOpen;
	      TsitsiklisUpdate(npoint);
	      Tree_PushIn(npoint);
          continue;
	    }
        //--------------------------------------------------------------
        if (S[npoint] == kDead)
        {
            if (V[npoint] != v)
            {
                //--------------------------------------------------------------
                // midpoint of constrained edge
                if ((e != NULL) && (e->midpoint == 0) && (_C[npoint] != NULL)) 
                {
                    e->midpoint = point;
                }
                //--------------------------------------------------------------
                // for Voronoi vertex detection
                
            }
            
        }
      }
  }
  //--------------------------------------------------------------
  nb_inserted_points -= Seeds->size();
  float upmax;
  while (nb_inserted_points >= 0)
  {
      pmax = NextPoint();
      if (pmax == -1) break;
      upmax = U[pmax];
      if ((upmax == 0.0) || (upmax < Ustop)) break;
      //--------------------------------------------------------------
      if (RunLocalPropagation(pmax) == true) nb_inserted_points--;
      else // insert on the encroached constrained edge e
      {
          while (!_edges_to_split->empty())
          {
              e = _edges_to_split->front();
              _edges_to_split->pop();
              insertOnEdge(e);
              //std::cerr << "Insertion on edge ";
              //print_point(e->midpoint);
              //std::cerr << " " << e->label1 << " " << e->label2 << std::endl;
              for (point = 0; point < size;  point++)
              {
                  if (S[point] != kBorder)
                  {
                      if (S[point] != kSeed) S[point] = kFar;
                      Tree[point] = -1;
                  }
                  
              }
              PtrToTrial = Trial - 1;
              nb_inserted_points--;
          }
      }
  }
  pmax = NextPoint();
  while (!validateTriangulation()) _triangles->clear();
}

//================================================================
void exportTriangulation(float *mat)
//================================================================
{
    if (mat == NULL) return;
    std::list<Triangle>::iterator it = _triangles->begin(), end = _triangles->end();
    int i = 0;
    for (; it != end; ++it)
    {
      mat[i++] = (float)it->l1;
      mat[i++] = (float)it->l2;
      mat[i++] = (float)it->l3;
    }
}
    
//--------------------------------------------------------------
void pos_to_point(int p, int *x, int *y)
//--------------------------------------------------------------
{
    *y = p / Nx;
    *x = p - Nx * (*y);
}

//---------------------------------------------------------
int NbSeeds()
//--------------------------------------------------------------
{
    return Seeds->size();
}
//--------------------------------------------------------------
int NbTriangles()
//--------------------------------------------------------------
{
    return _triangles->size();
}

//--------------------------------------------------------------
void GetSeeds(double *tab)
//--------------------------------------------------------------
{
    std::map<float,int>::iterator it = Seeds->begin(), end = Seeds->end();
    int i = 0, x, y, p;
    while (it != end)
    {
        p = it->second;
        pos_to_point(p, &x, &y);
        tab[i++] = (double)x;
        tab[i++] = (double)y;
        ++it;
    }
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
            U[point] = U[Point];
            //L[point] = L[Point];
            dUx[point] = dUx[Point];
            dUy[point] = dUy[Point];
            V[point] = V[Point];
        }
}

