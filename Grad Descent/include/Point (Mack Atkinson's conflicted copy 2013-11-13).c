//This is now general for all sets
#include "Point.h"

//given a point returns the index of the cube it belongs in
//"index" is an array of indeces, 0,0,0 representing the "lower-left" corner
//r is cutoff radius, a is array of components of "lower-left" corner of cube
int CubeIndex(pt thisPt,double r,double* a,int numCubes,int dim){
  int ans[dim];
  int i;
  for(i = 0; i < dim; i++){
    ans[i] = floor(((thisPt.components)[i] - a[i]) / r);
  }
  int realAns = 0;
  int k;
  for(k = 0; k < dim; k++){
    realAns+= ans[k] * pow(numCubes, k);
  }
  return realAns;
}

//finds the the cube index that the point belongs to
void trinaryLists(int *trinaryList,int dim){
  int i, j, count = 0, step = 0;
  for(i = 0; i < pow(3, dim); i++){
    step = i;
    for(j = 0; j < dim; j++){
      trinaryList[count] = step%3 - 1;
      count++;
      step = step/3;
    }
  }
}

double dist(pt a,pt b,int dim){
  int i;
  double totalsq =0;
  for(i = 0; i < dim; i++){
    totalsq += pow(a.components[i] - b.components[i], 2);
  }
  return sqrt(totalsq);
}

double cutoff(double distance, double radius){
  double t, value;
  t = distance/radius;
  if(t<1){
    value= pow(1-pow(t,4),3);
  }
  else{
    value = 0;
  }
  return value;
}

double cutoffPrime( double distance, double radius ){
  double t,value;
  t=distance/radius;
  if (t<1){value=3*pow(1-pow(t,4),2)*(-2*t*t);}
  else {value=0;}
  return value;
}









