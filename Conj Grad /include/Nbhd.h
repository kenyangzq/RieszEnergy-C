//This is now general for all sets
//works for all dimensions now

#ifndef NBHD_H_INCLUDED
#define NBHD_H_INCLUDED

#include "Cube.h"

typedef struct{
  int index;
  int neighborcount;
  int *neighbors;
} nbhd;


void CubeNbhd(cube* cubes,int numCubes,int* cubeIndex,cube* neighborhood,int* trinaryList,nbhd* cubeNbhd,int* testIndex,int dim);

void AssignNbhds(cube* cubes,double r,int numPts,int numCubes,pt* pts,nbhd* ans,int* nbrsArray,cube* neighborhood,int maxNbrs,int* trinaryList,nbhd* cubeNbhd,int* testIndex,int dim);

#endif
