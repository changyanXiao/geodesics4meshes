//================================================================
//================================================================
// File: fm2dSubGradiend.h
// (C) 02/2010 by Fethallah Benmansour
//================================================================
//================================================================

#include "fm.h"

#define kDead -1
#define kOpen -2
#define kFar -3
#define kBorder -4

/* Global variables */
int nx;			// real size on X
int ny;			// real size on Y
int Nx, Ny; // size for computing
int size;
float hx, hy;// spacing
float hx2, hy2;
float hx2hy2;
float hx2_plus_hy2;
float* U = NULL;// action map
float*  dUx = NULL;
float*  dUy = NULL;
short* S = NULL; // states
float* W = NULL; // potential
bool* Obstacle = NULL;
float* OutputGradient = NULL;
double* WW = NULL;
int   connectivity_small;
int*    NeighborhoodSmall = NULL;
float w; // regularisation parameter
double* start_points = NULL;
int START_POINT;
int* END_POINTS = NULL;
double* end_points    = NULL;
int nb_start_points;
int nb_end_points;


//================================================================
struct POINT{
    float* Grad;
    unsigned short nb_dead_neighbors;
};

POINT* Gradients = NULL;
//================================================================
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
  if(position)
  {
	if (ceil((float)position/2)==((float)position/2))
	  return (position/2 -1);
	else 
	  return((position-1)/2);
  }
  else
	return -1;
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
  while((position!=0)&&(U[s_idx]<U[f_idx]))
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
  else 
	return -1;
};

//================================================================
int Tree_GetLeftSon(int position)
//================================================================
{
  if ((2*position+1) < (int)(PtrToTrial - Trial +1))
	return 2*position+1 ;
  else 
	return -1;
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
  while((position >= 0)&&(!stop))
  {		
	if((Tree_GetRightSon(position)>0)&&(Tree_GetLeftSon(position)>0))
	{
	  int ls_idx = Trial[Tree_GetLeftSon(position)];
	  int rs_idx = Trial[Tree_GetRightSon(position)];
	  if( U[ls_idx] <= U[rs_idx] )
	  {
		Trial[position] = Trial[Tree_GetLeftSon(position)];
		Tree[Trial[position]] = position;
		position = Tree_GetLeftSon(position);
	  }
	  else
	  {
		Trial[position] = Trial[Tree_GetRightSon(position)];				
		Tree[Trial[position]] = (position);
		position = Tree_GetRightSon(position);				
	  }
	}
	else
	if(Tree_GetLeftSon(position)>0)
	{
	  Trial[position] = Trial[Tree_GetLeftSon(position)];				
	  Tree[Trial[position]] = (position);
	  position = Tree_GetLeftSon(position);
	}
	else 
	  stop = true;
  }
  if(position != (PtrToTrial - Trial))
  {
	Tree[*PtrToTrial] = position;
	Trial[position]=*PtrToTrial;
	int s_idx = Trial[position];
	int f_idx = Trial[Tree_GetFather(position)];
	while((position!=0)&&(U[s_idx]<U[f_idx]))
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
  else 
	return NULL;
};

//================================================================
void Tree_UpdateChange(int position)
//================================================================
{
  int s_idx = Trial[position];
  int f_idx = Trial[Tree_GetFather(position)];
  while((position!=0)&&(U[s_idx]<U[f_idx]))
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
  U[Tree_PopHead()]=Uv;
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
	connectivity_small = 4;
	NeighborhoodSmall = new int[connectivity_small];
	NeighborhoodSmall[ 0]= -Nx;
	NeighborhoodSmall[ 1]= -1;
	NeighborhoodSmall[ 2]=  1;
	NeighborhoodSmall[ 3]=  Nx;
};

//================================================================
void InitializeArrays()
//================================================================
{
	int x, y, point;
	//copy the weight list and initialize
	W = new float[size];
	S = new short[size];
    Gradients = new POINT[size];
    dUx = new float[size];
    dUy = new float[size];
    Tree = new int[size];
	Trial = new int[size];
    //------------------------------------------------------------
    for(x = 0; x < nx; x++){
		for(y = 0; y < ny; y++){
			point = (x+1) + (y+1)*Nx;
			W[point] = WW[x + y*nx] + w; S[point] = kFar;
		}
	}
	for(x = 0; x < size; x++){
		U[x] = INFINITE;
        Tree[x]=-1;
        Gradients[x].nb_dead_neighbors = 0;
        if(Obstacle[x])
            S[x] = kDead;
	}
    //------------------------------------------------------------
	PtrToTrial = Trial - 1;
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
};

//================================================================
void InitializeOpenHeap()
//================================================================
{
	int k, i, j, s, npoint;

    for(  s=0; s<nb_start_points; ++s ){
		i = round(start_points[2*s]);
		j = round(start_points[1+2*s]);
		START_POINT = i + j*Nx;
		//--------------------------------------------------------
		if(START_POINT >=size)
			mexErrMsgTxt("start_points should in the domaine.");
		//--------------------------------------------------------
		if( U[START_POINT]==0 )
			mexErrMsgTxt("start_points should not contain duplicates.");
		//--------------------------------------------------------
		U[START_POINT] = 0.0; S[START_POINT] = kOpen;
        Gradients[START_POINT].Grad = new float[size];
        Gradients[START_POINT].nb_dead_neighbors = 0;
        for(k = 0; k < size; k++)
            Gradients[START_POINT].Grad[k] = 0.0;
        for(k=0; k<connectivity_small; k++){
            npoint = START_POINT+NeighborhoodSmall[k];
            if( (S[npoint]==kDead) || (S[npoint]==kBorder) )
                Gradients[START_POINT].nb_dead_neighbors++;
        }
        //--------------------------------------------------------
        // add to heap
		Tree_PushIn(START_POINT);
		//--------------------------------------------------------
	}
    END_POINTS = new int[nb_end_points];
    //--------------------------------------------------------
    for(  s=0; s<nb_end_points; s++ ){
		i = round(end_points[2*s]);
		j = round(end_points[1+2*s]);
		END_POINTS[s] = i + j*Nx;
        if(END_POINTS[s] >=size)
			mexErrMsgTxt("start_points should be in the domaine.");
    }
};

//================================================================
float SethianQuadrant(float Pc,float Ux,float Uy)
//================================================================
{
    float Ua,Ub,qa,qb,Delta;
    
	float result;
	if (Ux<Uy){
		Ua = Ux;  qa = hx2;
		Ub = Uy;  qb = hy2;
	}
	else{
		Ua = Uy;  qa = hy2;
		Ub = Ux;  qb = hx2;
	}
	result = INFINITE;
	if ((sqrt(qa)*Pc)>(Ub-Ua)){
		Delta = (qa*qb)*((qa+qb)*Pc*Pc-(Ua-Ub)*(Ua-Ub));
		if (Delta>=0)
			result = ((qb*Ua+qa*Ub)+sqrt(Delta))/ hx2_plus_hy2;
	}
	else
		result = Ua+sqrt(qa)*Pc;
	return result;
};

//================================================================
void SethianQuadrantGradient(float* result_gradient, float Pc,float Ux,float Uy)
//================================================================
{
	float qx, qy, Delta;
	bool  is_quadratic = 0;
	float Gx = Ux+hx*Pc;
	float Gy = Uy+hy*Pc;
	if (Gx<Gy){
		result_gradient[0]=Gx;
		result_gradient[1]=Pc;
		result_gradient[2]=0;
	}
	else{
		result_gradient[0]=Gy;
		result_gradient[1]=0;
		result_gradient[2]=Pc;
	}
	if (Ux<Uy)
		is_quadratic= (Gx>Uy);
	else
		is_quadratic= (Gy>Ux);
	if (is_quadratic){
		qx=hx2;
		qy=hy2;
		Delta = hx2hy2*(hx2_plus_hy2*Pc*Pc-(Ux-Uy)*(Ux-Uy));
		if (Delta>=0){
			float result = ((qx*Uy+qy*Ux)+sqrt(Delta)) / hx2_plus_hy2;
			result_gradient[0] = result;
			result_gradient[1] = (result-Ux)/hx;
			result_gradient[2] = (result-Uy)/hy;
		}
	}
};

//================================================================
bool SethianUpdate(int point)
/*
COMMENTS : 
*/
//================================================================
{
	int	  npoint, npoint1, npoint2;
	float*  neighborU= new float [connectivity_small];
	float*  q_gradient = new float [3];
	float	  Pc = W[point];
	float	  Ur = U[point];
	float   dUxr = dUx[point];
	float   dUyr = dUy[point];
	bool is_updated = false;
	//--------------------------------------------------------------
	// Get the U & L values for each neighbor.
	for (int i=0;i<connectivity_small;i++){
		npoint=point+NeighborhoodSmall[i];
		if (S[npoint]==kDead){
			neighborU[i]=U[npoint];
		}
		else{			
			neighborU[i]=INFINITE;
		}
	}
	//--------------------------------------------------------------
	// Quadrant 1 : (x-1,y);(x,y);(x,y-1).
	SethianQuadrantGradient(q_gradient, Pc, neighborU[1], neighborU[0]);
	if (q_gradient[0]<Ur){
		Ur   =  q_gradient[0];
		dUxr =  q_gradient[1];
		dUyr =  q_gradient[2];
		is_updated = true;
	}
	//--------------------------------------------------------------
	// Quadrant 2 : (x+1,y);(x,y);(x,y-1).
	SethianQuadrantGradient(q_gradient, Pc, neighborU[2], neighborU[0]);
	if (q_gradient[0]<Ur){
		Ur   =  q_gradient[0];
		dUxr = -q_gradient[1];
		dUyr =  q_gradient[2];;
		is_updated = true;
	}
	//--------------------------------------------------------------
	// Quadrant 3 : (x-1,y);(x,y);(x,y+1).
	SethianQuadrantGradient(q_gradient, Pc, neighborU[1], neighborU[3]);
	if (q_gradient[0]<Ur){
		Ur   =  q_gradient[0];
		dUxr =  q_gradient[1];
		dUyr = -q_gradient[2];
		is_updated = true;
	}
	//--------------------------------------------------------------
	// Quadrant 4 : (x+1,y);(x,y);(x,y+1).
	SethianQuadrantGradient(q_gradient, Pc, neighborU[2], neighborU[3]);
	if (q_gradient[0]<Ur){
		Ur   =  q_gradient[0];
		dUxr = -q_gradient[1];
		dUyr = -q_gradient[2];
		is_updated = true;
	}
	//--------------------------------------------------------------
	if (is_updated){
		U[point]=Ur;
		dUx[point]=dUxr;
		dUy[point]=dUyr;
	}
	//--------------------------------------------------------------
	delete [] q_gradient; 
	delete [] neighborU;
	//--------------------------------------------------------------
	return is_updated;
};

//================================================================
void ComputeGradient(int point)
//================================================================
{
    //--------------------------------------------------------------
    int k, s, npointX, npointY;
    bool parentX, parentY;
    float alpha, beta;
    float dUxr = dUx[point];
    float dUyr = dUy[point];
    bool delete_condi = false;
    //--------------------------------------------------------------
    Gradients[point].Grad = new float[size];
    for(k = 0; k < size; k++)
        Gradients[point].Grad[k] = 0.0;
    //--------------------------------------------------------------
    if(dUxr > 0){
        npointX = point + NeighborhoodSmall[1];
        parentX = true;
    }
    else if(dUxr < 0){
        npointX = point + NeighborhoodSmall[2];
        parentX = true;
    }
    else parentX = false;
    //--------------------------------------------------------------
    if(dUyr > 0){
        npointY = point + NeighborhoodSmall[0];
        parentY = true;
    }
    else if(dUyr < 0){
        npointY = point + NeighborhoodSmall[3];
        parentY = true;
    }
    else parentY = false;
    //--------------------------------------------------------------
    Gradients[point].nb_dead_neighbors = 0;
    for(k=0; k < connectivity_small; k++)
        if(S[point + NeighborhoodSmall[k]] == kBorder)
            Gradients[point].nb_dead_neighbors++;
    //--------------------------------------------------------------
    if(parentX && parentY){
        if( (S[npointX]!=kDead) || (S[npointY]!=kDead) )
            mexErrMsgTxt("Parents must be kDead");
        alpha = (U[point] - U[npointX])/hx2;
        beta  = (U[point] - U[npointY])/hy2;
        for(k=0; k < size; k++)
            if(S[k] == kDead)
                Gradients[point].Grad[k] = ((k==point)*W[point] + alpha*Gradients[npointX].Grad[k]+ beta*Gradients[npointY].Grad[k])
                    / (alpha + beta);
        Gradients[point].nb_dead_neighbors = Gradients[point].nb_dead_neighbors + 2;
        Gradients[npointX].nb_dead_neighbors++;
        Gradients[npointY].nb_dead_neighbors++;
        delete_condi = (Gradients[npointX].nb_dead_neighbors == 4);
        for(s=0; s <nb_end_points; s++)
             delete_condi = delete_condi && (npointX != END_POINTS[s]);
        if(delete_condi)
            DELETEARRAY(Gradients[npointX].Grad);
        delete_condi = (Gradients[npointY].nb_dead_neighbors == 4);
        for(s=0; s <nb_end_points; s++)
             delete_condi = delete_condi && (npointY != END_POINTS[s]);
        if(delete_condi)
            DELETEARRAY(Gradients[npointY].Grad);
    }
    //--------------------------------------------------------------
    else if(parentX){
        if( S[npointX]!=kDead )
            mexErrMsgTxt("Parents must be kDead");
        for(k=0; k < size; k++)
            if(S[k] == kDead)
                Gradients[point].Grad[k] = hx*(k==point) + Gradients[npointX].Grad[k];
        Gradients[point].nb_dead_neighbors++;
        Gradients[npointX].nb_dead_neighbors++;
        delete_condi = (Gradients[npointX].nb_dead_neighbors == 4);
        for(s=0; s <nb_end_points; s++)
             delete_condi = delete_condi && (npointX != END_POINTS[s]);
        if(delete_condi)
            DELETEARRAY(Gradients[npointX].Grad);
    }
    //--------------------------------------------------------------
    else if(parentY){
        if( S[npointY]!=kDead )
            mexErrMsgTxt("Parents must be kDead");
       for(k=0; k < size; k++)
            if(S[k] == kDead)
                Gradients[point].Grad[k] = hy*(k==point) + Gradients[npointY].Grad[k];
        Gradients[point].nb_dead_neighbors++;
        Gradients[npointY].nb_dead_neighbors++;
        delete_condi = (Gradients[npointY].nb_dead_neighbors == 4);
        for(s=0; s <nb_end_points; s++)
             delete_condi = delete_condi && (npointY != END_POINTS[s]);
        if(delete_condi)
            DELETEARRAY(Gradients[npointY].Grad);
    }
    //--------------------------------------------------------------
    else{
        mexErrMsgTxt("recently fixed point must have at least one parent: could not be an orphen");
    }        
};

//================================================================
void CorrectMaps()
//================================================================
{
	int point;
	for (point=0;point<size;point++)
		if (S[point]!=kDead || Obstacle[point])
			U[point]=0;
};

//================================================================
void RunPropagation()
//================================================================
{
	int point,npoint,k, s;
	bool is_updated = false;
    bool end_points_reached = false;
	//--------------------------------------------------------------
	while ( (Tree_GetSize()>0)  && (!end_points_reached) ){
		point = Tree_PopHead();
		if(S[point]!=kOpen)
			mexErrMsgTxt("err : point must be Open");
		S[point]=kDead;
        //--------------------------------------------------------------
        end_points_reached = ( S[END_POINTS[0]] == kDead );
        for(k = 1; k < nb_end_points; k++)
            end_points_reached = end_points_reached &&  (S[END_POINTS[k]] == kDead) ;
        //--------------------------------------------------------------
        if(point!=START_POINT)
            ComputeGradient(point);
        if(end_points_reached){
            for(s = 0; s < nb_end_points; s++)
                for(k=0; k < size; k++)
                    OutputGradient[k + s*size] = Gradients[END_POINTS[s]].Grad[k];
        }
        //--------------------------------------------------------------
		for (k=0;k<connectivity_small;k++){
			npoint = point+NeighborhoodSmall[k];
			//--------------------------------------------------------------
			if (S[npoint]==kOpen){
				is_updated = SethianUpdate(npoint);
				if(is_updated)
					Tree_UpdateChange(Tree[npoint]);
			}
			//--------------------------------------------------------------
			else if (S[npoint]==kFar){
				S[npoint] = kOpen;
				SethianUpdate(npoint);
				Tree_PushIn(npoint);
			}
			//--------------------------------------------------------------
		}
	}
	//--------------------------------------------------------------
    CorrectMaps();
};

//================================================================
void resize()
//================================================================
{
    int x, y, s, point, Point;
    for(y=0;y<ny;y++)
        for(x=0;x<nx;x++){
            point = x+y*nx;
            Point = (x+1)+(y+1)*Nx;
            U[point] = U[Point];
            dUx[point] = dUx[Point];
            dUy[point] = dUy[Point];
        }
    
    int ss = nx*ny;
    for(s = 0; s < nb_end_points; s++)
        for(y=0;y<ny;y++)
            for(x=0;x<nx;x++){
                point = x+y*nx + s*ss;
                Point = (x+1)+(y+1)*Nx + s*size;
                OutputGradient[point] = OutputGradient[Point];
    }
};

//================================================================
void cleanGradients(){
//================================================================
    int k;
    for(k=0; k<size; k++)
        if(Gradients[k].nb_dead_neighbors > 0)
            DELETEARRAY(Gradients[k].Grad);
};