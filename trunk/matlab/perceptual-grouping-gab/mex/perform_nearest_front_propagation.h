#ifndef _PERFORM_NEAREST_FRONT_PROPAGATION_H_
#define _PERFORM_NEAREST_FRONT_PROPAGATION_H_

#include <math.h>
#include "config.h"
#include <stdio.h>
#include <string.h>
#include <vector>
#include <list>
#include <algorithm>
#include<iostream>

// some global variables
extern int n;			// width
extern int p;			// height
extern double* D;
extern double* S;
extern double* W;
extern double* Q;
extern double* L;
extern double* start_points;
extern double* end_points;
extern double* H;
extern double* values;
extern int nb_iter_max;
extern int nb_start_points;
extern int nb_end_points;
extern double* adj_mtrx; // adjacency graph (subgraph of the Voronoi dual)

typedef bool (*T_callback_intert_node)(int i, int j, int ii, int jj);

// main function
void perform_nearest_front_propagation(T_callback_intert_node callback_insert_node = NULL);
int NbSaddlePoints();
double* SaddlePoints(double*);

#endif // _PERFORM_NEAREST_FRONT_PROPAGATION_H_