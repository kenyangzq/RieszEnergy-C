//This is now general for all sets
//works for all dimensions now

#include "Nbhd.h"

//This creates all of the nbhds, empty


// questions:
// what's trinaryList?
// what do these two methods do?
// what are the fields of Nbhd class?


void CubeNbhd(cube* cubes,int numCubes,int* cubeIndex,cube* neighborhood,int* trinaryList, nbhd* cubeNbhd,int* testIndex,int dim){
  int count = 1;
  int flag;
  int i;
  int j;
  int k;
  int a,b,c;
  //    int testIndex0;
  //    int testIndex1;
  //    int testIndex2;
  int step;
  int index = 0;
  int zeroCount = 0;
  for(i = 0; i < dim; i++){
    index+= cubeIndex[i]*pow(numCubes, i);
  }
  a=cubeIndex[1];
  b=cubeIndex[2];
  c=cubeIndex[0];
  neighborhood[0] = cubes[index];
  for(i = 0; i < pow(3, dim); i++){
    flag = 0;
    zeroCount = 0;
    for( j = 0; j < dim; j++){
      if(trinaryList[i*dim + j] == 0){
        zeroCount++;
      }
      step = cubeIndex[j] + trinaryList[i*dim + j];
      testIndex[j] = step;
      if(step < 0 || step >= numCubes || zeroCount==dim) flag++;
    }
    if(flag==0){
      index = 0;
      for(k = 0; k < dim; k++){
        index+= testIndex[k] * pow(numCubes, k);
      }
      neighborhood[count] = cubes[index];
      count++;
    }
  }
  int upperBound = 0;
  for(i = 0; i < count; i++){
    upperBound += neighborhood[i].pointCount;
  }
  int* allPts = cubeNbhd->neighbors;
  int l;
  int count2 = 0;
  for(l = 0; l < count; l++){
	  int* root;
    root = neighborhood[l].pointList;
    for( i = 0; i < neighborhood[l].pointCount; i++){
      allPts[count2] = neighborhood[l].pointList[i];
      count2++;
    }
  }
  cubeNbhd->neighborcount = count2;
}




//this assigns the corrects points to each neighborhood by looping through each point, and finding all points within a certain radius of it
void AssignNbhds(cube* cubes,double r,int numPts,int numCubes,pt* pts,nbhd* ans,int* nbrsArray,cube* neighborhood,int maxNbrs,int* trinaryList,nbhd* cubeNbhd,int* testIndex,int dim){

  int i=0;    
  int j=0;
  int k=0;
  int count;
  int pointCount;
  int* allPts;
  pt centerPoint; 
  pt measurePoint;
  nbhd pointNeighborhood;
  int coords[dim];
  int cubePtsNum;
  double disty ;
  int abc = 0;
  int step;
  trinaryLists(trinaryList,dim);
  //go through each cube
  for(i = 0; i < pow(numCubes, dim); i++){
    step = i;
    for(abc = 0; abc < dim; abc++){
      coords[abc] = step % numCubes;
      step = step / numCubes;
    }
    CubeNbhd(cubes, numCubes, coords, neighborhood, trinaryList, cubeNbhd, testIndex,dim);
    pointCount = cubeNbhd->neighborcount;
    allPts = cubeNbhd->neighbors;
    cubePtsNum = cubes[i].pointCount;
    //go through each point in the cube
    for(j = 0; j < cubePtsNum; j++){
      count = 0;
      centerPoint = pts[allPts[j]];
      //go through each point in the cube neighborhood of the point
      for( k = 0; k < pointCount; k++){
        measurePoint = pts[allPts[k]];
        disty = dist(centerPoint, centerPoint,dim);
        if(dist(centerPoint, measurePoint,dim) < r && centerPoint.index != measurePoint.index){
          nbrsArray[centerPoint.index*maxNbrs+count] = measurePoint.index;
          count++;     
          if(count > maxNbrs){
            perror("\nWarning: maxNbrs is too small.\n");
            count = maxNbrs;
            break;                        
          }
        }
			//ended the k loop here, wasn't ended here before, pretty sure this is correct
      }
      pointNeighborhood.neighbors = NULL;
      if(count == 0){
        pointNeighborhood.neighborcount = 0;
        pointNeighborhood.neighbors = NULL;
        pointNeighborhood.index = centerPoint.index;
      }
			//included below in else, before it was an empty else, didn't make sense
      else{
	      pointNeighborhood.neighborcount = count;                
        pointNeighborhood.neighbors= &nbrsArray[centerPoint.index*maxNbrs];                
        pointNeighborhood.index = centerPoint.index;            
        ans[centerPoint.index] = pointNeighborhood;
      }  
    }
  }
}









