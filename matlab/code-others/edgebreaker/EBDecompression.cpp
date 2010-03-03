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

#include <vector>

// Added by GP
#define SEPARATOR "/"

using namespace std;

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

void initDecompression();
void DecompressConnectivity(int c);
bool CheckHandle(int c);
void Zip(int c);
void DecompressVertices(int c);
Vertex DecodeDelta(int c);

/***************************** Variables *********************************/

//Input arguments
static char sInputDirectory[MAX_SIZE];
static char sOutFileName[MAX_SIZE];
static MeshType eMeshType = MANIFOLD;
static FileFormat eFileFormat = ASKII;

//File Input
static FILE* fclers = NULL;			//clers file (See File Formats for detais)
static FILE* fvertices = NULL;		//vertices (See File Formats for detais)
static FILE* fhandles = NULL;		//handles pairs (See File Formats for detais)

//Variables for storing Input OVTable and geometry
static int nNumOfTriangles;			//Number of Triangles 		
static int nNumOfVertices;			//Number of Vertices 				

static int*	O = NULL;				//Output Opposite table
static int*	V = NULL;				//Output Vertex indices table
static Vertex*	G = NULL;			//Output Geometry table


//Compression variables
static int T = 0;					//triangles count
static int N = 0;					//vertices count
static int I = 0;					//number of triangles on first vertex
static int A = 0;					//handles count
static int *M = NULL;				//Vetex marking array
static int *U = NULL;				//Triangles marking array

//Handle corners array
static vector <int> H;

//Clers array
static vector <char> C;

//Geometry array
static vector <Vertex> G_in;

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
	if (argc != 4)
	{
		printf("Wrong number of agruments.\n\n"
			   "EBDecompression InputFileDir OutFileName FileFormat\n\n");
		exit(0);
	}


	strcpy(sInputDirectory, argv[1]);
	strcpy(sOutFileName, argv[2]);

	char sFileFormat[MAX_SIZE];
	strcpy(sFileFormat, argv[3]);

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
void OpenInputFiles()
{
	char sTempStr[MAX_SIZE];

	
	strcpy(sTempStr, sInputDirectory);
	strcat(sTempStr, SEPARATOR);
	strcat(sTempStr, "clers.txt");

	fclers = fopen(sTempStr, "rt");
	if (fclers == NULL)
		PrintErrorAndQuit("Can not open clers file\n");

	strcpy(sTempStr, sInputDirectory);
	strcat(sTempStr, SEPARATOR);
	strcat(sTempStr, "vertices.txt");

	fvertices = fopen(sTempStr, "rt");
	if (fvertices == NULL)
		PrintErrorAndQuit("Can not open vertices file\n");

	strcpy(sTempStr, sInputDirectory);
	strcat(sTempStr, SEPARATOR);
	strcat(sTempStr, "handles.txt");

	fhandles = fopen(sTempStr, "rt");
	if (fhandles == NULL)
		PrintErrorAndQuit("Can not open handles file\n");

}


//Write OV TAble.
void WriteVOTableFile(char* sFileName)
{
	int i = 0;
	int nNumOfTriangles = C.size() + 1;
	int nNumOfVertices = G_in.size();

	FILE* pOutFile = fopen(sFileName, "w+t");
	if (pOutFile == NULL)
		PrintErrorAndQuit("Can not open output file\n");
	
	//Write the number of vertices to the file
	if (!fprintf(pOutFile, "%d\n", nNumOfTriangles))
		PrintErrorAndQuit("Error writing to the file\n");
	
	//Write VO TAble from the file
	for (i = 0; i<nNumOfTriangles; i++)
	{
		if (!fprintf(pOutFile,
				  "%d %d\n", 
				  V[i*3],
				  O[i*3]))
			PrintErrorAndQuit("Error reading file\n");
		if (!fprintf(pOutFile,
				  "%d %d\n", 
				  V[i*3+1],
				  O[i*3+1]))
			PrintErrorAndQuit("Error reading file\n");
		if (!fprintf(pOutFile,
				  "%d %d\n", 
				  V[i*3+2],
				  O[i*3+2]))
			PrintErrorAndQuit("Error reading file\n");
	}


	//Write the number of vertices from the file
	if (!fprintf(pOutFile, "%d\n", nNumOfVertices))
		PrintErrorAndQuit("Error reading file\n");
	
	//Write all vertices to the file
	for (i = 0; i<nNumOfVertices; i++)
	{
		if (!fprintf(pOutFile,
				  "%f %f %f\n", 
				  G[i].x,
				  G[i].y,
				  G[i].z))
			PrintErrorAndQuit("Error reading file\n");
	}

	//Close the file.
	fclose(pOutFile);
}

/*
* Usage: EBDecompression InputFileDir OutFileName FileFormat
*		 
*	InputFileDir: Directory with clers.txt handles.txt vertices.txt 
*	OutFileName: Name of output file 
*	FileFormat: BINARY or ASKII (See File Formats for detais)
*
*/
int main(int argc, char* argv[])
{

	//Process arguments
	ProcessArguments(argc, argv);

	//Open input files
	OpenInputFiles();

	//Compress Mesh
	initDecompression();

	//write reconstructed model to the file
	WriteVOTableFile(sOutFileName);

	//clear memory
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

/***************************** File Read/Write functions *****************/
void myScanf(FILE* fFile, char* format, void* var)
{
	fscanf(fFile,  format, var);
//	if (fscanf(fFile,
//			   format, 
//			  var) == -1)
//		PrintErrorAndQuit("Error reading file\n");

}

/*
* Read input vertives into G_in array
*/
void ReadInputVertices()
{
	//Read input vertices
	if (eFileFormat == ASKII)
	{
		Vertex vTemp;
		myScanf(fvertices, "%f", (void*)&vTemp.x);
		myScanf(fvertices, "%f", (void*)&vTemp.y);
		myScanf(fvertices, "%f", (void*)&vTemp.z);
		do
		{
			G_in.push_back(vTemp);
			myScanf(fvertices, "%f", (void*)&vTemp.x);
			myScanf(fvertices, "%f", (void*)&vTemp.y);
			myScanf(fvertices, "%f", (void*)&vTemp.z);
		}while(!feof(fvertices));
	}
}


/*
* Read Handles pairs into H array
*/
void ReadHandlesPairs()
{
	//Read handle corners
	if (eFileFormat == ASKII)
	{
		int nTemp;
		myScanf(fhandles, "%d", (void*)&nTemp);
		do
		{
			H.push_back(nTemp);
			myScanf(fhandles, "%d", (void*)&nTemp);
		}while(!feof(fhandles));
	}
}

/*
* Read Mesh Type
*/
void ReadMeshType()
{
	char sMeshType[MAX_SIZE];
	if (eFileFormat == ASKII)
	{
		//read mesh type
		myScanf(fclers, "%s", (void*)sMeshType);

		//remember mesh type
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
	}
}

/*
* Read numbet of triangles incident on first vertex
*/
void ReadNumTriOnFirstVertex()
{
	if (eFileFormat == ASKII)
	{
		//read numbet of triangles 
		myScanf(fclers, "%d", (void*)&I);
	}
	
}

/*
* Read clers into C array,
* Append (I-2) Cs and 1R to the beginning of the array
*/
void ReadClearsAndUppendClersForFirstVertexTri()
{
	if (eFileFormat == ASKII)
	{
		//Write (I-2) Cs and 1R to the beginning of the array

		//First write I-2 C's
		for (int i = 0; i < (I-2); i++)
		{
			C.push_back('C');
		}

		//now write one R
		C.push_back('R');

		//now read clers from input file
		char cTemp; 
		//read new line character
		myScanf(fclers, "%c", (void*)&cTemp);
		//read clers symbol
		myScanf(fclers, "%c", (void*)&cTemp);

		do
		{
			C.push_back(cTemp);

			//read new line character
			myScanf(fclers, "%c", (void*)&cTemp);

			//read clers symbol
			myScanf(fclers, "%c", (void*)&cTemp);

		}while (!feof(fclers));
	}

	
}

void InitDecompressionModule()
{
	int nNumOfTriangles = C.size() + 1;
	int nNumOfVertices = G_in.size();
	int i;

	//Allocate Memory for Tables
	V = new int[3*nNumOfTriangles];		//table of vertex Ids for each corner
	O = new int[3*nNumOfTriangles];		//table of opposite corner Ids for each corner or free edge orientation

	//Init memory
	for (i = 0; i<3*nNumOfTriangles; i++) {O[i] = -3; V[i] = 0;}

	//Initialize the first triangle 
	V[0] = 0; V[1] = 2; V[2] = 1;
//	V[0] = 0; V[1] = 1; V[2] = 2;
	O[0] = -1; O[1] = -1;

	//Allocate memory for M and U tables
	U = new int[nNumOfTriangles];
	M = new int[nNumOfVertices];

	//init tables for marking visited vertices and triangles
	for (i = 0; i < nNumOfTriangles; i++) U[i] = 0;
	for (i = 0; i < nNumOfVertices; i++) M[i] = 0;

	//Allocate memory for geometry array
	G = new Vertex[nNumOfVertices];
}

/***************************** EB Decompression **************************/
/*
* Decompress the mesh
*/
void initDecompression()
{
	int i=0;

	//Initialization of O and V tables is done in InitDecompression function
	//aftet we read clers from the file and findout how many tryangles are in 
	//the mesh.
	
	//id of the last triangle decompressed so far
	T = 0;								
	//id of the last vertex encountered
	N = 2;							
	//id of the last handle encountered
	A = 0;
	
	//read input vertices
	ReadInputVertices();

	//read handles pairs into H array
	ReadHandlesPairs();

	//read meshType from clers file
	ReadMeshType();

	//read number of incident triangles on first vertex
	ReadNumTriOnFirstVertex();
	
	//read clers into C array,
	//append (I-2) Cs and 1R to the beginning of the array
	ReadClearsAndUppendClersForFirstVertexTri();

	//Init Memory and some variables
	InitDecompressionModule();

	//start connectivity decompression
	DecompressConnectivity(2);						


	//Initialization of tables for marking visited vertices and triangles
	//Done in InitDecompression function

	//estimate 1st  vertex
	G[0] = DecodeDelta(0);	
	//if we do not have a hole mark 1st  vertex as visited
	if (eMeshType == MANIFOLD) M[0] = 1;

	//estimate third vertex and mark it as visited
	G[1] = DecodeDelta(2);	
	//G[1] = DecodeDelta(1);	
	M[1] = 1;

	//estimate second vertex and mark it as visited
	G[2] = DecodeDelta(1);	
	//G[2] = DecodeDelta(2);	
	M[2] = 1;

	//id of the last vertex encountered
	N=2;
	//paint the triangle and go to opposite corner
	U[0] = 1; 
	
	//start vertices decompression
	DecompressVertices(O[2]);

}

void DecompressConnectivity(int c) 
{
	//Loop builds triangle tree and zips it up
	do
	{
		//new triangle
		T++;							

		//attach new triangle, link opposite corners
		O[c] = 3*T;	O[3*T] = c;			

		//enter vertex Ids for shared vertices
		V[3*T+1] = V[PrevEdge(c)];		
		V[3*T+2] = V[NextEdge(c)];

		//move corner to new triangle
		c = NextEdge(O[c]);	

		//select operation based on next symbol
		switch (C[T-1])		
		{
		case 'C':						
			//C: left edge is free, store ref to new vertex			
			O[NextEdge(c)] = -1;
			V[3*T] = ++N;
			break;
		case 'L':						
			//L: orient free edge
			O[NextEdge(c)] = -2;
			//check for handles, if non, try to zip
			if (!CheckHandle(NextEdge(c)))
				Zip(NextEdge(c));
			break;
		case 'R':						
			//R: orient free edge, check for handles, go left 
			O[c] = -2;
			CheckHandle(c);
			c = NextEdge(c);
			break;
		case'S':						
			//O[NextEdge(c)] = -10;
			//S: recursion going right, then go left
			DecompressConnectivity(c);
			c = NextEdge(c);
			//if the triangle to the left was visited, then return	
			if (O[c] >= 0)
				return;
			break;
		case 'E':						
			//E: left and right edges are  free
			O[c] = -2;
			O[NextEdge(c)] = -2;
			//check for handles on the right
			CheckHandle(c);
			//check for handles on the left, if non, try to zip
			if (!CheckHandle(NextEdge(c)))
				Zip(NextEdge(c));			
			//pop 	
			return;
			break;
		}
	}while(true);
}


bool CheckHandle(int c)
{
	//check if this is a handle
	//if (A < H.size() && c == H[A+1])
	if (A >= H.size() || c != H[A+1])
		return false;
	else
	{
		//link opposite corners
		O[c] = H[A];
		O[H[A]] = c;

		//find corner of next free edge if any 
		int a = PrevEdge(c);
		while( (O[a] >=0) && (a != H[A]) )
			a = PrevEdge(O[a]);

		//zip if found cw edge
		if (O[a] == -2)
			Zip(a);

		//find corner of next free edge if any
		a = PrevEdge(O[c]);
		while( (O[a] >=0) && (a != c) )
			a = PrevEdge(O[a]);

		//zip if found cw edge
		if (O[a] == -2)
			Zip(a);

		//next handle
		A+=2;
		return true;
	}
}

void Zip(int c) 
{
	//tries to zip free edges opposite c
	int b = NextEdge(c);	

	//search clockwise for free edge
	while(O[b] >= 0 && O[b] != c) 
	{
		b = NextEdge(O[b]);
	}

	//pop if no zip possible
	if (O[b] != -1)
	{
		return;
	}

  	//link opposite corners
	O[c] = b; O[b] = c;

	//assign co-incident corners
	int a = NextEdge(c);
	V[NextEdge(a)]=V[NextEdge(b)];
	while(O[a]>=0 && a!=b)
	{
		a = NextEdge(O[a]);
		V[NextEdge(a)]=V[NextEdge(b)];
	}

	//find corner of next free edge on right
	c=PrevEdge(c);
	while(O[c] >=0 && c!= b)
	{	
		c=PrevEdge(O[c]);
	}

	//try to zip again
	if (O[c] == -2)
		Zip(c);
}

void DecompressVertices(int c)
{
	//start traversal for triangle tree
	do 
	{
		//mark the triangle as visited
		U[E2T(c)] = 1; 

		//test whether tip vertex was visited
		if (M[V[c]] == 0)	
		{	
			//update new vertex
			N++;
			G[N] = DecodeDelta(c);	
			//mark tip vertex as visited
			M[V[c]] = 1;
			//continue with the right neighbor
			c = RightTri(c, O);
		}										
		else 
		//test whether right triangle was visited
		if (U[E2T(RightTri(c, O))] == 1)
		{
			//test whether left triangle was visited
			if (U[E2T(LeftTri(c, O))] == 1)
			{
				//E, pop
		        return; 
			}									
			else 
			{
				//R,move to left triangle
				c = LeftTri(c, O);
			}
		}
		else 
		//test whether left triangle was visited
		if (U[E2T(LeftTri(c, O))] == 1)
		{
			//L, move to right triangle
			c = RightTri(c, O);
		}	
		else 
		{
			//S, recursive call to visit right branch first
			DecompressVertices(RightTri(c, O));
			//move to left triangle
			c = LeftTri(c, O);
			//if the triangle to the left was visited, then return
			if (U[E2T(c)]>0)
				return;
		} 
	}while(true);

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
* This function assume vertices was not encoded using any prediction 
* It just reads the next vertex 
*/
Vertex DecodeNoPrediction(int c)
{
	return G_in[V[c]];
}

Vertex DecodeWithPrediction(int c)
{
	static int EBVcount = 0;
	
	Vertex vPred;
	Vector delta = G_in[EBVcount];
	EBVcount++;

	if (M[V[O[c]]] > 0 && M[V[PrevEdge(c)]] > 0) 
	{
		vPred = VPlus(G[V[NextEdge(c)]], G[V[PrevEdge(c)]]);
		vPred = VMinus(vPred, G[V[O[c]]]);
		vPred = VPlus(delta, vPred);
		return vPred;
	}
	if (M[V[O[c]]] > 0) 
	{
		vPred = VMult(G[V[NextEdge(c)]], 2);
		vPred = VMinus(vPred, G[V[O[c]]]);
		vPred = VPlus(delta, vPred);
		return vPred;
	}
	if (M[V[NextEdge(c)]] > 0 && M[V[PrevEdge(c)]] > 0) 
	{
		vPred = VPlus(G[V[NextEdge(c)]], G[V[PrevEdge(c)]]);
		vPred = VMult(vPred, 0.5f);
		vPred = VPlus(delta, vPred);
		return vPred;
	}
	if (M[V[NextEdge(c)]] > 0) 
	{
		vPred = G[V[NextEdge(c)]];
		vPred = VPlus(delta, vPred);
		return vPred;
	}
	else if (M[V[PrevEdge(c)]] > 0) 
	{
		vPred = G[V[PrevEdge(c)]];
		vPred = VPlus(delta, vPred);
		return vPred;
	}
	else
	{
		vPred = delta;
		return vPred;
	}

}

Vertex DecodeDelta(int c)
{
	return DecodeNoPrediction(c);
	//return DecodeWithPrediction(c);
	
}

