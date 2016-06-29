//works for all dimensions now
//This is now general for all sets
#include "Cube.h"



int findNumCubes(double r,double* a,double* b){
  int numCubes;    
  numCubes = ceil((b[0] - a[0]) / r);
  return numCubes;
}

//creates empty n-dim array
//array is actually 1 dim, access [a,b,c...] by [a + b*numCubes + c*numCubes*numCubes...]
//where numCubes is the number of cubes in one side length
//r is cutoff radius, a and b are opposite corners of cube
cube* CreateCubes(double r,double* a,double* b,cube* cubes,int* cubePtsArray,int maxNbrs,int dim){
  int numCubes;
  numCubes = findNumCubes(r, a, b);
  int j;
  for(j = 0; j<pow(numCubes, dim);j++){
    cube cube1;
    cube1.pointList = &cubePtsArray[j*maxNbrs];
    cube1.index = j;
    cube1.pointCount = 0;
    cubes[j] = cube1;
  }
  return cubes;
}

//adds point to correct cube
void pt2cube(double r,double* a,cube* cubes,int numCubes,pt thisPt,int dim){
  int cubeInd = CubeIndex(thisPt, r, a, numCubes,dim);
  cubes[cubeInd].pointList[cubes[cubeInd].pointCount] = thisPt.index;
  cubes[cubeInd].pointCount++;
}

void AssignCubes(double r,double* a,double* b,pt* pts,int numpts,cube* cubes,int* cubePtsArray,int maxNbrs,int dim){
  cubes = CreateCubes(r, a, b, cubes, cubePtsArray, maxNbrs,dim);
  int numCubes;
  numCubes = findNumCubes(r, a, b);
  //go through each point, assign it to the cube
  int i;
	int cubeInd;
  for(i = 0; i < numpts; i++){
    pt2cube(r, a, cubes, numCubes, pts[i],dim);
  }
}
