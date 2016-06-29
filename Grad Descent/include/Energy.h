//Works for all dimensions now
//This is now general for all sets

#ifndef ENERGY_H_INCLUDED
#define ENERGY_H_INCLUDED

#include "Nbhd.h"

double ptenergy(pt *pts, pt thePt, nbhd theNbhd, double s, double radius,int dim);

double ptenergies(pt *pts, int numPts, nbhd *nbhds, double s, double radius, double* ptenarray,int dim);

double energy(pt *pts, int numPts, nbhd *nbhds, double s, double radius,int dim);

void ptgrad( pt *pts, pt thePt, nbhd theNbhd, double s, double radius , double *theGrad,int dim);

void grad( pt *pts, int numpts, nbhd *ptsnbhd, double *gradarray, double s, double radius,int dim);

void gradSphere( pt *pts, int numpts, nbhd *ptsnbhd, double *gradarray, double s, double radius,int dim);

void gradBall( pt *pts, int numpts, nbhd *ptsnbhd, double *gradarray, double s, double radius,int dim);

void gradShell( pt *pts, int numpts, nbhd *ptsnbhd, double *gradarray, double s, double radius,int dim,double inner_rad);

void randptBall(double* coordinates, int dim);

void randptShell(double* coordinates, int dim, double inner_rad);

void randptSphere(double* coordinates, int dim);

void randptsBall(double* coords, int numPts,int dim);

void randptsShell(double* coords, int numPts,int dim, double inner_radius);

void randptsSphere(double* coords, int numPts, int dim);

double sEnergy(pt *points, double s, double radius, double* cornera, double* cornerb, int numpts, cube *cubes, int numCubes, nbhd *nbhds, int *neighborArray, cube *neighborhood, int maxNbrs, int *trinaryList, nbhd* cubeNbhd, int* testIndex, int* cubePtsArray,int dim);

double sEnGradBall(pt *points, double *theGrad, double s, double radius, double* cornera, double* cornerb, int numpts, cube *cubes, int numCubes, nbhd *nbhds, int *neighborArray, cube *neighborhood, int maxNbrs, int *trinaryList, nbhd* cubeNbhd, int* testIndex, int* cubePtsArray,int dim, double inner_rad);

double sEnGradSphere(pt *points, double *theGrad, double s, double radius, double cornera[3], double cornerb[3], int numpts, cube *cubes, int numCubes, nbhd *nbhds, int *neighborArray, cube *neighborhood, int maxNbrs, int *trinaryList, nbhd* cubeNbhd, int* testIndex, int* cubePtsArray,int dim, double inner_rad);

double sEnGradShell(pt *points, double *theGrad, double s, double radius, double cornera[3], double cornerb[3], int numpts, cube *cubes, int numCubes, nbhd *nbhds, int *neighborArray, cube *neighborhood, int maxNbrs, int *trinaryList, nbhd* cubeNbhd, int* testIndex, int* cubePtsArray,int dim,double inner_rad);

int roundcount;

int done;

#endif
