// FastMarching2D.cpp : Defines the entry point for the DLL application.
//

#include "fm.h"
#include<iostream>
short kDead = -1;
short kOpen = -2;
short kFar = -3;
short kBorder = -4;

/* Global variables */
int nx;			// real size on X
int ny;			// real size on Y
int nxny;
int Nx, Ny;     // size for computing
int size;       // number of points
float hx, hy;   // spacing
float hx2, hy2, hxhy;
float hx2hy2;
float hx2_plus_hy2;
float sqrt_hx2_plus_hy2;
float* U = NULL;// action map
float*  dUx = NULL;
float*  dUy = NULL;
float* L = NULL; // distance map
short* S = NULL; // states
short* V = NULL; // voronoi
float* M1 = NULL; // M1
float* M2 = NULL; // M1
float* M3 = NULL; // M1
double* T = NULL;
//================================================================
int connectivity_large;
int* NeighborhoodLarge = NULL;
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

//================================================================
void InitializeNeighborhoods()
//================================================================
{
	connectivity_large = 8;
    // Freeman code of the 8-adjacency pixel graph
	NeighborhoodLarge = new int[connectivity_large+1];
	NeighborhoodLarge[0]= -Nx-1;
	NeighborhoodLarge[1]= -Nx;
	NeighborhoodLarge[2]= -Nx+1;
	NeighborhoodLarge[3]=  -1; //1;
    NeighborhoodLarge[4]=  1; //Nx+1;
    NeighborhoodLarge[5]=  Nx-1;//Nx;
    NeighborhoodLarge[6]=  Nx; //Nx-1;
    NeighborhoodLarge[7]=  Nx+1;//-1;
    NeighborhoodLarge[8]=  -Nx-1;
    //-----------------------------------------------------------
    Simplicies1 = new int  [connectivity_large];
    Simplicies2 = new int  [connectivity_large];
    signsX      = new int  [connectivity_large];
    signsY      = new int  [connectivity_large];
    h1          = new float[connectivity_large];
    h2          = new float[connectivity_large];
    h1_h2       = new float[connectivity_large];
    h22         = new float[connectivity_large];
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
  //copy the weight list and initialize
  S = new short[size];   // ok
  Tree = new int[size];  // ok
  Trial = new int[size]; // ok
  M1 = new float[size];  // ok
  M2 = new float[size];  // ok
  M3 = new float[size];
  q_gradient = new float[3];
  //------------------------------------------------------------
  for (x = 0; x < nx; x++)
  {
    for (y = 0; y < ny; y++)
    {
      point = (x+1) + (y+1)*Nx;
      M1[point] = T[x + y*nx];
      M2[point] = T[x + y*nx + 3*nxny];
      M3[point] = T[x + y*nx + nxny];
      if (M1[point] == 0 || M2[point] == 0)
      { V[point] = kBorder; S[point] = kBorder; }
      else { V[point] = kDead; S[point] = kFar; }
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
    U[x] = INFINITE;
    L[x] = INFINITE;
    Tree[x]=-1;
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
void InitializeOpenHeap()
//================================================================
{
  int point, i, j;
  
  for (int s = 0; s < nb_start_points; ++s)
    {
      i = round(start_points[2*s]);
      j = round(start_points[1+2*s]);
      point = i + j*Nx;
      //--------------------------------------------------------
      if (point >= size)
	mexErrMsgTxt("start_points should be in the domain.");
      //--------------------------------------------------------
      if (U[point] == 0)
	mexErrMsgTxt("start_points should not contain duplicates.");
      //--------------------------------------------------------
      U[point] = 0.0; L[point] = 0.0; S[point] = kOpen; V[point] = s;
      // add to heap
      Tree_PushIn(point);
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
short TsitsiklisQuadrantGradient(float* result_gradient, float m1, float m2, float m3,
                                 float Ua, float Ub, short Va, short Vb, int QuadNb)
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
		result_gradient[0] = Ub + sqrt(r22);
        result_gradient[1] = signsX[QuadNb]/sqrt(r22);
        result_gradient[2] = signsY[QuadNb]/sqrt(r22);
        return Vb;
	}
    //-----------------------------------------------------------
    else
        if(k1 <= -sqrt(r11)) {
            if(QuadNb < 4) {
                result_gradient[0] = Ua + sqrt(m1);
                result_gradient[1] = signsX[QuadNb]/sqrt(m1);
                result_gradient[2] = 0.0;
            }
            else {
                result_gradient[0] = Ua + sqrt(m2);
                result_gradient[1] = 0.0;
                result_gradient[2] = signsY[QuadNb]/sqrt(m2);
            }
            return Va;
        }
        //-----------------------------------------------------------
        else{
            //-------------------------------------------------------
            if (r12 <= k1*sqrt(R/(r11-k1*k1))) {
                result_gradient[0] = Ub + sqrt(r22);
                result_gradient[1] = signsX[QuadNb]/sqrt(r22);
                result_gradient[2] = signsY[QuadNb]/sqrt(r22);
                return Vb;
            }
            //-------------------------------------------------------
            else if (r12 > (r11 + k1*sqrt(R/(r11-k1*k1)))) {
                if (QuadNb < 4) {
                    result_gradient[0] = Ua + sqrt(m1);
                    result_gradient[1] = signsX[QuadNb]/sqrt(m1);
                    result_gradient[2] = 0.0;
                }
                else {
                    result_gradient[0] = Ua + sqrt(m2);
                    result_gradient[1] = 0.0;
                    result_gradient[2] = signsY[QuadNb]/sqrt(m2);
                }
                return Va;
            }
            //-------------------------------------------------------
            else {
                alpha = (r12 - k1 * sqrt(R/(r11-k1*k1))) / r11;
                result_gradient[0] = alpha*k1 + k2 + sqrt(R/(r11-k1*k1));
                if (QuadNb < 4) {
                    result_gradient[1] = signsX[QuadNb] / sqrt(R/(r11-k1*k1));
                    result_gradient[2] = (1-alpha) * signsY[QuadNb] / sqrt(R/(r11-k1*k1));
                }
                else {
                    result_gradient[1] = (1-alpha) * signsX[QuadNb] / sqrt(R/(r11-k1*k1));
                    result_gradient[2] = signsY[QuadNb] / sqrt(R/(r11-k1*k1));
                }
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
	float Ur = U[point];
	short Vr, Vtmp;
	float Lr = L[point];
	float dUxr = dUx[point];
	float dUyr = dUy[point];
	bool is_updated = false;
	//--------------------------------------------------------------
	// Get the U & L values for each neighbor.
	for (int i = 0; i < connectivity_large; i++)
    {
		npoint1 = point + Simplicies1[i];
        npoint2 = point + Simplicies2[i];
		if (S[npoint1] == kDead) { Ua = U[npoint1]; La = L[npoint1]; }
        else { Ua = INFINITE; La = INFINITE; }
        if (S[npoint2] == kDead) { Ub = U[npoint2]; Lb = L[npoint2]; }
        else { Ub = INFINITE; Lb = INFINITE; }
        Vtmp = TsitsiklisQuadrantGradient(q_gradient, M1[point], M2[point], M3[point],
                                          Ua, Ub, V[npoint1], V[npoint2], i);
        if (q_gradient[0] < Ur)
        {
            Ur   =  q_gradient[0];
    		dUxr =  q_gradient[1]; dUyr =  q_gradient[2];
            TsitsiklisQuadrantLength(&Lr, 1, La, Lb, i);
            Vr = Vtmp;
            is_updated = true;
    	}
	}
	//--------------------------------------------------------------
	if (is_updated)
    {
		U[point] = Ur;
		L[point] = Lr;
		V[point] = Vr;
		dUx[point] = dUxr; dUy[point] = dUyr;
	}
	//--------------------------------------------------------------
	return is_updated;
};

//================================================================
void CorrectMaps()
//================================================================
{
	int point;
	for (point = 0; point < size; point++) {
		if (V[point] < 0) { U[point] = 0; L[point] = 0; }
    }
};

//================================================================
void RunPropagation()
//================================================================
{
  int point,npoint,k;
  bool is_updated = false;
  double Lmax = 0.0, _umax = 0.0;
  
  //--------------------------------------------------------------
  while (Tree_GetSize() > 0)
    {
        
      point = Tree_PopHead();
      //if (_umax > U[point]) std::cerr << _umax << " " << U[point] << std::endl;
      //else _umax = U[point];
      if(S[point] != kOpen) mexErrMsgTxt("err");
      S[point] = kDead;
      Lmax = MAX(Lmax, L[point]);
      if ((Dmax > 0) && (Lmax > Dmax)) break;
      //--------------------------------------------------------------
      for (k = 0; k < connectivity_large; k++)
	  {
        npoint = point + NeighborhoodLarge[k];
        //--------------------------------------------------------------
        if (S[npoint] == kOpen)
	    {
	      is_updated = TsitsiklisUpdate(npoint);
	      if (is_updated) Tree_UpdateChange(Tree[npoint]);
	    }
        //--------------------------------------------------------------
        else if (S[npoint] == kFar)
	    {
	      S[npoint] = kOpen;
	      TsitsiklisUpdate(npoint);
	      Tree_PushIn(npoint);
	    }
	  //--------------------------------------------------------------
      }
  }
  _umax=0;
  int pp;
  for (int x = 0; x < size;  x++)
  {
      if (S[x] != kBorder && _umax < U[x]) {_umax=U[x]; pp = x;}
      
  }
  int y = pp / Nx;
  int x = pp - Nx * y;
  //std::cerr << _umax << " [" << x << ";" << y << "]" << std::endl;
  //--------------------------------------------------------------
  //CorrectMaps();
};

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
            L[point] = L[Point];
            dUx[point] = dUx[Point];
            dUy[point] = dUy[Point];
            V[point] = V[Point];
        }
};
