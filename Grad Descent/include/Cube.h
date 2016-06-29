//works for all dimensions now
//This is now general for all sets
//header file for the cube functions

#ifndef CUBE_H_INCLUDED
#define CUBE_H_INCLUDED

#include "Point.h"

typedef struct{
  int* pointList;
  int index;
  int pointCount;
} cube;

int findNumCubes(double r,double* a,double* b);

cube* CreateCubes(double r,double* a,double* b,cube* cubes,int* cubePtsArray,int maxNbrs,int dim);

void pt2cube(double r,double* a,cube* cubes,int numCubes,pt thisPt,int dim);

void AssignCubes(double r,double* a,double* b,pt* pts,int numpts,cube* cubes,int* cubePtsArray,int maxNbrs,int dim);









#endif
