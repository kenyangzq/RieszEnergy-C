


// This is now general for all sets
#include "Point.h"

// int CubeIndex(pt thisPt, double r, double* a, int numCubes, int dim)
// Description: given a point returns the index of the cube it belongs in
//              "index" is an array of indeces, 0,0,0 representing the "lower-left" corner
//              r is cutoff radius, a is array of components of "lower-left" corner of cube
// Parameter: thisPt -- the point
//            r -- the cutoff radius
//            a -- the lower left corner of the cube
//            numCubes -- the number of cubes in the cutoff
//            dim -- the dimension
// return: the index of the cube,
//
int CubeIndex(pt thisPt, double r, double* a, int numCubes, int dim){
    int ans[dim], realAns = 0;
    for(int i = 0; i < dim; i++){
        ans[i] = floor((thisPt.components[i] - a[i]) / r);
        realAns += ans[i] * pow(numCubes, i);
    }
    return realAns;
}

//creates list of all combiniations from 0 to 1 for arrays of size dim
void trinaryLists(int*trinaryList,int dim){
    int count = 0;
    for(int i = 0; i < pow(3, dim); i++){
        int step = i;
        for(int j = 0; j < dim; j++){
            trinaryList[count] = step%3 - 1;
            count++;
            step = step/3;
        }
    }
}

//finds the distance between two points
double dist(pt a,pt b,int dim){
    double totalsq =0;
    for(int i = 0; i < dim; i++){
        totalsq += pow(a.components[i] - b.components[i], 2);
    }
    return sqrt(totalsq);
}

//defines cutoff when calculating energy
double cutoff(double distance, double radius){
    double t, value;
    t = distance/radius;
    if(t<1){
        value= pow(1-pow(t,4),3);
    }else{
        value = 0;
    }
    return value;
}


double cutoffPrime( double distance, double radius ){
    double t, value;
    t= distance / radius;
    if(t<1){
        value=3*pow(1-pow(t,4),2)*(-2*t*t);
    }else{
        value=0;
    }
    return value;
}










