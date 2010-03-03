/**************************************************************************
*
* Any use of this code or of the ideas it embodies should be acknowledged
* in writing in all publications describing the resulting product or
* research contribution as follows: "This solution is based on the
* Edgebreaker compression technique and software developed by Prof.
* Rossignac and his colleagues at Georgia Institute of Technology."    
*
*
* Software developed by: Alla Safonova at Georgia Tech
* Last modified: 05/11/20001
*
**************************************************************************/

//Include
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>


/***************************** Types *************************************/
struct Coord3D
{
	float x;
	float y;
	float z;
};
typedef Coord3D Vertex; 
typedef Coord3D Vector;

#define MAX_SIZE 256

enum MeshType {MANIFOLD, TPATCH};
enum FileFormat {BINARY, ASKII};

void initCompression(int c, MeshType eMeshType);
void Compress(int c);
void CheckHandle(int c);
void EncodeDelta(int c);

/***************************** Variables *********************************/

//Input arguments
static char sOVTable[MAX_SIZE];
static MeshType eMeshType = MANIFOLD;
static FileFormat eFileFormat = ASKII;
static int nStartCorner = 0;

//File Output
static FILE* fclers = NULL;			//clers file (See File Formats for detais)
static FILE* fvertices = NULL;		//vertices (See File Formats for detais)
static FILE* fhandles = NULL;		//handles pairs (See File Formats for detais)

//Variables for storing Input OVTable and geometry
static int*	O = NULL;				//Input Opposite table
static int*	V = NULL;				//Input Vertex indices table
static Vertex*	G = NULL;			//Input Geometry table
static Vertex*	G_est = NULL;			//Input Geometry table


//Compression variables
static int T = 0;					//triangles count
static int N = 0;					//vertices count
static int *M = NULL;				//Vetex marking array
static int *U = NULL;				//Triangles marking array



/****************************** Main Block *******************************/

void ClearMemoryAndFiles()
{
	//close files
	if (fclers != NULL)
		fclose(fclers);
	if (fvertices != NULL)
		fclose(fvertices);
	if (fhandles != NULL)
		fclose(fhandles);

	//disallocate memory
	if (O != NULL)
		delete [] O;
	if (V != NULL)
		delete [] V;
	if (G != NULL)
		delete [] G;
	if (U != NULL)
		delete [] U;
	if (M != NULL)
		delete [] M;

}

//Print Error string and exit.
void PrintErrorAndQuit(char* sErrorString)
{
	printf(sErrorString);
	ClearMemoryAndFiles();
	exit(0);
}

void ProcessArguments(int argc, char* argv[])
{
	if (argc != 5)
	{
		printf("Wrong number of agruments.\n\n"
			   "Usage: EBCompression OVTable MeshType FileFormat\n\n");
		exit(0);
	}

	char sMeshType[MAX_SIZE];
	char sFileFormat[MAX_SIZE];
	char sStartCorner[MAX_SIZE];

	strcpy(sOVTable, argv[1]);
	strcpy(sMeshType, argv[2]);
	strcpy(sStartCorner, argv[3]);
	strcpy(sFileFormat, argv[4]);

	//get mesh type
	if (strcmp("MANIFOLD", sMeshType)==0)
	{
		eMeshType = MANIFOLD;
		printf("MeshType - MANIFOLD\n");
	}
	else if (strcmp("TPATCH", sMeshType)==0)
	{
		eMeshType = TPATCH;
		printf("MeshType - TPATCH\n");
	}
	else
		PrintErrorAndQuit("Not supported mesh type\n");

	//get start corner
	nStartCorner = atoi(sStartCorner);

	//get output file format
	if (strcmp("BINARY", sFileFormat)==0)
	{
		eFileFormat = BINARY;
		printf("FileFormat - BINARY\n");
	}
	else if (strcmp("ASKII", sFileFormat)==0)
	{
		eFileFormat = ASKII;
		printf("FileFormat - ASKII\n");
	}
	else
		PrintErrorAndQuit("Not supported file format\n");

}

//Open filr for storing clers, geometry and handles
void OpenOutputFiles()
{
	fclers = fopen("clers.txt", "w+t");
	if (fclers == NULL)
		PrintErrorAndQuit("Can not open clers file\n");

	fvertices = fopen("vertices.txt", "w+t");
	if (fvertices == NULL)
		PrintErrorAndQuit("Can not open vertices file\n");

	fhandles = fopen("handles.txt", "w+t");
	if (fhandles == NULL)
		PrintErrorAndQuit("Can not open handles file\n");

}

//Process input file and build Vertex and Triangle Arrays.
void ProcessInputFile(char* sFileName)
{
	int i = 0;
	FILE* pInFile = NULL;

	//Number of Triangles 	
	int nNumOfTriangles;
	//Number of Vertices 			
	int nNumOfVertices;							
	
	//Open the file.
	pInFile = fopen(sFileName, "rt");
	if (pInFile == NULL)
		PrintErrorAndQuit("Can not open input file\n");
	
	//Read number of triangles from the file
	if (fscanf(pInFile, "%d", &nNumOfTriangles) == -1)
		PrintErrorAndQuit("Error reading file\n");

	
	//Allocate memory for V Table
	V = new int [nNumOfTriangles*3]; 

	//Allocate memory for O Table
	O = new int [nNumOfTriangles*3]; 

	//Read VO TAble from the file
	for (i = 0; i<nNumOfTriangles; i++)
	{
		if (fscanf(pInFile, "%d %d", &(V[i*3]), &(O[i*3]))==-1)
			PrintErrorAndQuit("Error reading file\n");
		if (fscanf(pInFile, "%d %d", &(V[i*3+1]), &(O[i*3+1])) == -1)
			PrintErrorAndQuit("Error reading file\n");
		if (fscanf(pInFile, "%d %d", &(V[i*3+2]), &(O[i*3+2])) == -1)
			PrintErrorAndQuit("Error reading file\n");
	}


	//Read number of vertices from the file
	if (fscanf(pInFile, "%d", &nNumOfVertices) == -1)
		PrintErrorAndQuit("Error reading file\n");
	
	//Allocate memory for vertex array
	G = new Vertex [nNumOfVertices]; 
	G_est = new Vertex [nNumOfVertices]; 

	//Read all vertices from the file
	for (i = 0; i<nNumOfVertices; i++)
	{
		if (fscanf(pInFile, "%f %f %f", &(G[i].x), &(G[i].y), &(G[i].z)) == -1)
			PrintErrorAndQuit("Error reading file\n");
	}

	//Allocate memory for M and U tables
	M = new int [nNumOfVertices];	//Table for marking visited vetrices
	U = new int [nNumOfTriangles];	//Table for marking visited triangles

	//init tables for marking visited vertices and triangles
	for (i = 0; i<nNumOfVertices; i++) M[i] = 0;
	for (i = 0; i<nNumOfTriangles; i++) U[i] = 0;

	//Close the file.
	fclose(pInFile);
}


/*
* Usage: EBCompression OVTable MeshType FileFormat
*		 
*	OVTable: Connectivity and geometry of the mesh in OVTable format 
*			(See File Formats for detais)
*	MeshType: 2 Mesh types are currently supported - MANIFOLD and TPATCH
*			  MANIFOLD - is a manifold mesh, consistently oriented with no holes. 		
*			  TPATCH - is a manifold mesh with boundary, consistently oriented.
*	StartCorner: Is the corner where to start EBCompression. 
*				If  MeshType is TPATCH it should be a corner corresponding to 
*					the "dummy" vertex.
*				If  MeshType is MANIFOLD it can be any corner, but since the 
*					triangles incident on  StartCorner are not stored it is 
*					advantageous to pass a corner that has maximum number of 
*					triangles incident on it as a StartCorner.
*	FileFormat: BINARY or ASKII
*			(See File Formats for detais)
*
*/
int main(int argc, char* argv[])
{

	//Process arguments
	ProcessArguments(argc, argv);

	//Open output files
	OpenOutputFiles();

	//Read OVTableFile
	ProcessInputFile(sOVTable);

	//Compress Mesh
	initCompression(nStartCorner, eMeshType);

	//disallocate memory
	ClearMemoryAndFiles();
	return 0;
}

/***************************** EB Helper Functions ***********************/


int NextEdge(int edge)
{
	return (3*(edge / 3) + (edge + 1) % 3);
}

int PrevEdge(int edge)
{
	return NextEdge(NextEdge(edge));
}


int RightTri(int c, int* O_table)
{
	//c.r = c.n.r 
	return O_table[NextEdge(c)];
}

int LeftTri(int c, int* O_table)
{
	//c.l = c.n.n.r 
	return O_table[NextEdge(NextEdge(c))];
}

int E2T(int edge)
{
	return (edge / 3);
}

/***************************** EB Compression ****************************/

/*
*	Arguments:
*		c - start compression from corner c
*		MeshType: 2 Mesh types are currently supported - MANIFOLD and TPATCH
*					MANIFOLD - is a manifold mesh, consistently oriented with no holes. 		
*					TPATCH - is a manifold mesh with boundary, consistently oriented.
*		FileFormat: BINARY or ASKII (See File Formats for detais)
*
*/
void initCompression(int c, MeshType eMeshType)
{
	int i = 0;

	//init tables for marking visited vertices and triangles
	//was done in ProcessInputFile function

	//id of the last triangle compressed so far
	T = 0;

	c = PrevEdge(c);	

	//estimate 1st  vertex
	EncodeDelta(NextEdge(c));
	
	//if we do not have a hole mark 1st  vertex as visited
	//in which case estimate function can use it for estimation
	//if we do have a hole, we do not mark 1st  vertex as visited
	//and it is not used for estimation since it is a dummy vertex
	if (eMeshType==MANIFOLD) M[V[NextEdge(c)]] = 1;


	//estimate third vertex and mark it as visited
	EncodeDelta(c);
	M[V[c]] = 1;

	//estimate second vertex and mark it as visited
	EncodeDelta(PrevEdge(c));
	M[V[PrevEdge(c)]] = 1;

	//paint the triangle 
	U[E2T(c)] = 1; // mark the triangle as visited

	//traverse triangles incident on the first vertex
	//we do not want to store clers symbols for them
	int a = O[c];

	//we keep a count of number of triangles incident on the first corner
	int count = 1;
	//first traverse 'C' triangles 
	while (a != PrevEdge(O[PrevEdge(c)]))
	{
		//increment count for number of triangles incident on the first corner
		count++;

		//paint the triangle, increment # of triangles 
		U[E2T(a)] = 1;
		T++;

		//estimate next vertex and mark it as visited
		EncodeDelta(a);
		M[V[a]] = 1;
		
		//continue with the right neighbor 
		a = O[NextEdge(a)];
	}

	//traverse 'R' triangle incident on first vertex 
	U[E2T(a)] = 1;
	T++;
	count++;
	
	//write mesh type to clers file
	if (eMeshType == MANIFOLD)
	{
		if (eFileFormat == ASKII)
			fprintf(fclers, "%s\n", "MANIFOLD");
	}
	else if (eMeshType == TPATCH)
	{
		if (eFileFormat == ASKII)
			fprintf(fclers, "%d\n", "TPATCH");
	}
	else
		PrintErrorAndQuit("Not supported mesh type\n");

	//write number of triangles incident on first vertex to clers file
	if (eFileFormat == ASKII)
		fprintf(fclers, "%d\n", (int)count);
	
	//start connectivity compression
	Compress(O[PrevEdge(a)]);
}



void Compress(int c)
{
	//start traversal for triangle tree
	do 
	{
		//mark the triangle as visited
		U[E2T(c)] = 1; 
		T++;
		
		//check for handles
		CheckHandle(c);

		//test whether tip vertex was visited
		if (M[V[c]] == 0)	
		{	
			//append encoding of C to clers
			fprintf(fclers, "%c\n", 'C');
			
			//estimate next vertex and mark it as visited
			EncodeDelta(c);
			M[V[c]] = 1;

			//continue with the right neighbor
			c = RightTri(c, O);
		}										
		else 
		//test whether right triangle was visited
		if (U[E2T(RightTri(c, O))] > 0)
		{
			//test whether left triangle was visited
			if (U[E2T(LeftTri(c, O))] > 0)
			{
				//append code for E and pop
				fprintf(fclers, "%c\n", 'E');
		        return; 
			}		
			else 
			{
				//append code for R, move to left triangle
				fprintf(fclers, "%c\n", 'R');
				c = LeftTri(c, O);
			}
		}
		else 
		//test whether left triangle was visited
		if (U[E2T(LeftTri(c, O))] > 0)
		{
			//append code for L, move to right triangle
			fprintf(fclers, "%c\n", 'L');
			c = RightTri(c, O);
		}	
		else 
		{
			//store corner number in decompression, to support handles
			U[E2T(c)] = T*3+2;

			//append code for S
			fprintf(fclers, "%c\n", 'S');

			//recursive call to visit right branch first
			Compress(RightTri(c, O));

			//move to left triangle
			c = LeftTri(c, O);

			//if the triangle to the left was visited, then  return
			if (U[E2T(c)]>0)
				return;
		} 
	}while(true);

}


void CheckHandle(int c)
{
	//check for handles from the right
	if (U[E2T(O[NextEdge(c)])] >1)
	{
		//write opposite corners for handle triangles into file
		fprintf(fhandles, "%d %d\n", U[E2T(O[NextEdge(c)])], T*3+1);
	}

	//check for handles from the left
	if (U[E2T(O[PrevEdge(c)])] >1)
	{
		//write opposite corners for handle triangles into file
		fprintf(fhandles, "%d %d\n", U[E2T(O[PrevEdge(c)])], T*3+2);
	}
}


/******************* Vector Operations for Estimate functions ************/

//Returns v1 - v2
Vector VMinus(Vertex v1, Vertex v2)
{
	Vector tempVector;
	tempVector.x = v1.x - v2.x;
	tempVector.y = v1.y - v2.y;
	tempVector.z = v1.z - v2.z;

	return tempVector;
}


//Returns v1 + v2
Vector VPlus(Vertex v1, Vector v2)
{
	Vector tempVector;
	tempVector.x = v2.x + v1.x;
	tempVector.y = v2.y + v1.y;
	tempVector.z = v2.z + v1.z;

	return tempVector;
}

//Returns v1*k
Vector VMult(Vertex v1, float k)
{
	Vector tempVector;
	tempVector.x = v1.x*k;
	tempVector.y = v1.y*k;
	tempVector.z = v1.z*k;

	return tempVector;
}


/***************************** Estimate functions ************************/
/*
* This function does not do any prediction, it just writes vertices into array
*/
void EncodeNoPrediction(int c)
{
	//Store vertex coordinates into file
	fprintf(fvertices, "%f %f %f\n", G[V[c]].x,
									 G[V[c]].y,	
									 G[V[c]].z);	
}


void EncodeWithPrediction(int c)
{
	Vector vPred, delta;
	Vertex zeroV = {0.0, 0.0, 0.0};

	if (M[V[O[c]]] > 0 && M[V[PrevEdge(c)]] > 0) 
	{
		vPred = VPlus(G_est[V[NextEdge(c)]], G_est[V[PrevEdge(c)]]);
		vPred = VMinus(vPred, G_est[V[O[c]]]);
		delta = VMinus(G[V[c]], vPred);
		//return vPred;
	}
	else if (M[V[O[c]]] > 0) 
	{
		vPred = VMult(G_est[V[NextEdge(c)]], 2);
		vPred = VMinus(vPred, G_est[V[O[c]]]);
		delta = VMinus(G[V[c]], vPred);
		//return vPred;
	}
	else if (M[V[NextEdge(c)]] > 0 && M[V[PrevEdge(c)]] > 0) 
	{
		vPred = VPlus(G_est[V[NextEdge(c)]], G_est[V[PrevEdge(c)]]);
		vPred = VMult(vPred, 0.5f);
		delta = VMinus(G[V[c]], vPred);
		//return vPred;
	}
	else if (M[V[NextEdge(c)]] > 0) 
	{
		vPred = G_est[V[NextEdge(c)]];
		delta = VMinus(G[V[c]], vPred);
		//return vPred;
	}
	else if (M[V[PrevEdge(c)]] > 0) 
	{
		vPred = G_est[V[PrevEdge(c)]];
		delta = VMinus(G[V[c]], vPred);
		//return vPred;
	}
	else
	{
		vPred = zeroV;
		delta = VMinus(G[V[c]], vPred);
	}

	G_est[V[c]] = VPlus(delta, vPred);

	fprintf(fvertices, "%f %f %f\n", delta.x,
									 delta.y,	
									 delta.z);	
}


void EncodeDelta(int c)
{
	EncodeNoPrediction(c);
	//EncodeWithPrediction(c);
	
}

