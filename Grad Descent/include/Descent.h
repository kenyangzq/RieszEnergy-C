//This is now general for all sets

#ifndef DESCENT_H_INCLUDED
#define DESCENT_H_INCLUDED

#include "Energy.h"

void updateBall(pt *newpts, pt *oldpts, int numpts, double *pvec, double t,int dim,double inner);

void updateShell(pt *newpts, pt *oldpts, int numpts, double *pvec, double t,int dim,double inner);

void updateSphere(pt *newpts, pt *oldpts, int numpts, double *pvec, double t,int dim,double inner);

double GradDescent(pt *points, pt *newpoints, int numpts, double *theGrad, double t0, double radius, double* cornera, double* cornerb, cube *cubes, int numCubes, int *neighborArray, int maxNbrs, double s, cube *neighborhood, nbhd *nbhds, FILE *writeTo, int *trinaryList, nbhd* cubeNbhd, int* testIndex, int* cubePtsArray,int dim,double inner,void(*updateMethod)(pt *,pt *,int,double *,double,int,double),double(*gradMethod)(pt *,double *,double,double,double *,double *,int,cube *,int,nbhd *,int *,cube *,int,int *,nbhd *,int *,int *,int,double));

#endif
