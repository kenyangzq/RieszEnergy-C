// This is a header containting functions pertaining to points
// I have issues with this, I want to make the array a pointer, but there is a segmentation fault
// This is now general for all sets

#ifndef POINT_H_INCLUDED
#define POINT_H_INCLUDED

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>


//data structure for points; components is int array containing its components
typedef struct{
    int index;
//  double components[3];
	double* components;
} pt;

//finds the the cube index that the point belongs to
int CubeIndex(pt thisPt, double r, double* a, int numCubes, int dim);



//creates list of all combiniations from 0 to 1 for arrays of size dim
void trinaryLists(int *trinaryList,int dim);



double dist(pt a,pt b,int dim);

double cutoff(double distance, double radius);

double cutoffPrime( double distance, double radius);




#endif

