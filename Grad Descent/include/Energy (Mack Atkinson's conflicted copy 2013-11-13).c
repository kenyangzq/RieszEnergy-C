//This works for any dimension now 
//This is now general for all sets

#include "Energy.h"

double ptenergy(pt *pts, pt thePt, nbhd theNbhd, double s, double radius,int dim){
  int i;
  double energy = 0;
  double ptdist;
  for(i = 0; i < theNbhd.neighborcount; i++){
    ptdist = dist(pts[theNbhd.neighbors[i]], thePt,dim);
    energy += cutoff(ptdist, radius)/pow(ptdist, s);
  }
  return energy;
}

double ptenergies(pt *pts, int numPts, nbhd *nbhds, double s, double radius, double* ptenarray,int dim){
  int i;   
  double pten;
  double totalEn=0;
  for(i = 0; i < numPts; i++){
    pten = ptenergy(pts, pts[i], nbhds[i], s, radius,dim);
    ptenarray[i] = pten;
    totalEn += pten;
  }
  return totalEn;
}

double energy(pt *pts, int numPts, nbhd *nbhds, double s, double radius,int dim){
  int i,j;   
  double pten;
  double ptdist;
  double totalEn=0;
  nbhd theNbhd;
  for(i = 0; i < numPts; i++){
    pten=0;
    theNbhd=nbhds[i];
    for(j = 0; j < theNbhd.neighborcount; j++){
      if (theNbhd.neighbors[j]<i) {
        ptdist = dist(pts[theNbhd.neighbors[j]], pts[i],dim);
        pten += cutoff(ptdist, radius)/pow(ptdist, s);
      }
    }
    totalEn += pten;
  }
  return 2*totalEn;
}

void initialize(int numPts, int numCubes, cube* neighborhood, double* ptenarray ,nbhd* ans, int* tempArray, cube* cubes){
  neighborhood = malloc(27 * sizeof(cube));
  ptenarray = malloc(numPts * sizeof(double));
  ans = malloc( numPts * sizeof(nbhd));
  tempArray = malloc(numPts*sizeof(int));
  cubes = malloc(numCubes * numCubes * numCubes * sizeof(cube));    
}

void ptgrad( pt *pts, pt thePt, nbhd theNbhd, double s, double radius , double *theGrad,int dim){
  int i,j;
  double ptdist;
  double gradfactor;
  pt nbhrPt;
  for(j=0;j<dim;j++){
    theGrad[j]=0;
  }
	for(i=0;i<theNbhd.neighborcount;i++){
    nbhrPt=pts[theNbhd.neighbors[i]];
    ptdist=dist(nbhrPt, thePt,dim);
    gradfactor=-1*s*cutoff(ptdist, radius)/pow(ptdist, s+2)+2*cutoffPrime(ptdist,radius)/(radius*radius*pow(ptdist,s));
    for(j=0;j<dim;j++){
      theGrad[j]+=gradfactor*(thePt.components[j]-nbhrPt.components[j]);
    }
  }
}

void grad( pt *pts, int numpts, nbhd *ptsnbhd, double *gradarray, double s, double radius,int dim){
  int i;
  for(i=0;i<numpts;i++){
	  ptgrad(pts,pts[i],ptsnbhd[i], s, radius, &gradarray[i*dim],dim);
  }
}

void ptgradSphere( pt *pts, pt thePt, nbhd theNbhd, double s, double radius , double *theGrad,int dim){
  int i,j;
  double ptdist;
  double gradfactor;
  double graddotpt;
  pt nbhrPt;
  for(j=0;j<dim;j++){
    theGrad[j]=0;
  }
  for(i=0;i<theNbhd.neighborcount;i++){
    nbhrPt=pts[theNbhd.neighbors[i]];
    ptdist=dist(nbhrPt, thePt,dim);
    gradfactor=-1*s*cutoff(ptdist, radius)/pow(ptdist, s+2)+2*cutoffPrime(ptdist,radius)/(radius*radius*pow(ptdist,s));
    for(j=0;j<dim;j++){
      theGrad[j]+=gradfactor*(thePt.components[j]-nbhrPt.components[j]);
    }
    graddotpt=0;
    for(j=0;j<dim;j++){
      graddotpt+=theGrad[j]*thePt.components[j];
    }
    for(j=0;j<dim;j++){
      theGrad[j]+=-1.*graddotpt*thePt.components[j];
    }
  }
}

void gradSphere( pt *pts, int numpts, nbhd *ptsnbhd, double *gradarray, double s, double radius,int dim){
  int i;
  for(i=0;i<numpts;i++){
    ptgradSphere(pts,pts[i],ptsnbhd[i], s, radius, &gradarray[i*dim],dim);
  }
}

void randptBall(double* coordinates, int dim){
	int i;
	double normsq=2;
	double z;
	while(normsq>1){
		normsq=0;
		for(i=0;i<dim;i++){
			z=1-(2*(double)rand()/(double)RAND_MAX);
			normsq+=z*z;
  		coordinates[i] = z;
		}
	}
}

void randptShell(double* coordinates, int dim, double inner_rad){
	int i;
	double normsq=2.0;
	double inradsq;
	double z;
	int count=0;
	inradsq=inner_rad*inner_rad;
	while(normsq>1 || normsq<inradsq){
		normsq=0;
		for(i=0;i<dim;i++){
			z=1-(2*(double)rand()/(double)RAND_MAX);
			normsq+=z*z;
  		coordinates[i] = z;
		}
	}
}

void randptSphere(double* coordinates, int dim){
	int i;
	double normsq=2;
	double norm;
	double z;
	while(normsq>1 || normsq==0){
		normsq=0;
		for(i=0;i<dim;i++){
			z=1-(2*(double)rand()/(double)RAND_MAX);
			normsq+=z*z;
  		coordinates[i] = z;
		}	
	}
	norm=sqrt(normsq);
	for(i=0;i<dim;i++){
  	coordinates[i] = coordinates[i]/norm;
	}
}

void randptsBall(double* coords, int numPts,int dim){
	int i;
	for(i = 0; i < numPts; i++){
    randptBall(&coords[i*dim], dim);
  }
}

void randptsShell(double* coords, int numPts,int dim, double inner_radius){
  int i;
  for(i = 0; i < numPts; i++){
    randptShell(&coords[i*dim],dim, inner_radius);
  }
}

void randptsSphere(double* coords, int numPts, int dim){
  int i;
  for(i = 0; i < numPts; i++){
		randptSphere(&coords[i*dim],dim);
  }
}

double sEnergy(pt *points, double s, double radius, double* cornera, double* cornerb, int numpts, cube *cubes, int numCubes, nbhd *nbhds, int *neighborArray, cube *neighborhood, int maxNbrs, int *trinaryList, nbhd* cubeNbhd, int* testIndex, int* cubePtsArray,int dim){

  AssignCubes(radius, cornera, cornerb, points, numpts, cubes, cubePtsArray, maxNbrs,dim);
  AssignNbhds(cubes, radius, numpts, numCubes, points, nbhds, neighborArray, neighborhood, maxNbrs, trinaryList, cubeNbhd, testIndex,dim);
  return energy(points, numpts, nbhds, s, radius,dim);
}

double sEnGradBall(pt *points, double *theGrad, double s, double radius, double* cornera, double* cornerb, int numpts, cube *cubes, int numCubes, nbhd *nbhds, int *neighborArray, cube *neighborhood, int maxNbrs, int *trinaryList, nbhd* cubeNbhd, int* testIndex, int* cubePtsArray,int dim){
    
  AssignCubes(radius, cornera, cornerb, points, numpts, cubes, cubePtsArray, maxNbrs,dim);
  AssignNbhds(cubes, radius, numpts, numCubes, points, nbhds, neighborArray, neighborhood, maxNbrs, trinaryList, cubeNbhd, testIndex,dim);
  grad(points, numpts, nbhds, theGrad, s, radius,dim);
  return energy(points, numpts, nbhds, s, radius,dim);
}

double sEnGradSphere(pt *points, double *theGrad, double s, double radius, double cornera[3], double cornerb[3], int numpts, cube *cubes, int numCubes, nbhd *nbhds, int *neighborArray, cube *neighborhood, int maxNbrs, int *trinaryList, nbhd* cubeNbhd, int* testIndex, int* cubePtsArray,int dim){
  AssignCubes(radius, cornera, cornerb, points, numpts, cubes, cubePtsArray, maxNbrs,dim);
  AssignNbhds(cubes, radius, numpts, numCubes, points, nbhds, neighborArray, neighborhood, maxNbrs, trinaryList, cubeNbhd, testIndex,dim);
  gradSphere(points, numpts, nbhds, theGrad, s, radius,dim);
  return energy(points, numpts, nbhds, s, radius,dim);
}

