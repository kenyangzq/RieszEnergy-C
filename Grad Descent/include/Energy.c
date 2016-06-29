//This works for any dimension now 
//This is now general for all sets

#include "Energy.h"

//finds the energy of a single point
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

//finds energies of all the points, and adds it up
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

//finds the total energy by including a cutoff
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

//finds the gradient for a given point
void ptgrad( pt *pts, pt thePt, nbhd theNbhd, double s, double radius , double *theGrad,int dim){
  int i,j;
  double ptdist;
  double gradfactor;
	double norm;
	double err=pow(10,-12);
	double dot;
  pt nbhrPt;
  for(j=0;j<dim;j++){
    theGrad[j]=0;
  }
	for(i=0;i<theNbhd.neighborcount;i++){
    nbhrPt=pts[theNbhd.neighbors[i]];
    ptdist=dist(nbhrPt, thePt,dim);
    gradfactor=-1*s*cutoff(ptdist, radius)/pow(ptdist, s+2)+2*cutoffPrime(ptdist,radius)/(radius*radius*pow(ptdist,s));
		norm=0;
		dot=0;
    for(j=0;j<dim;j++){
      theGrad[j]+=gradfactor*(thePt.components[j]-nbhrPt.components[j]);
			norm=norm+thePt.components[j];
    }
  }
}

//finds gradient for all points in a general set
void grad( pt *pts, int numpts, nbhd *ptsnbhd, double *gradarray, double s, double radius,int dim){
  int i;
  for(i=0;i<numpts;i++){
	  ptgrad(pts,pts[i],ptsnbhd[i], s, radius, &gradarray[i*dim],dim);
  }
}

//finds gradient for all points in the sphere
void gradSphere( pt *pts, int numpts, nbhd *ptsnbhd, double *gradarray, double s, double radius,int dim){
  int i,j ;
 double graddotpt;
 
  for(i=0;i<numpts;i++){
    ptgrad(pts,pts[i],ptsnbhd[i], s, radius, &gradarray[i*dim],dim);
    
    graddotpt=0;
    for(j=0;j<dim;j++){
      graddotpt+=gradarray[i*dim+j]*pts[i].components[j];
    }
    for(j=0;j<dim;j++){
      gradarray[i*dim+j]+=-1.*graddotpt*pts[i].components[j];
    }

  }
}

//finds the gradient for all points in the ball
void gradBall( pt *pts, int numpts, nbhd *ptsnbhd, double *gradarray, double s, double radius,int dim){
  int i,j ;
  double graddotpt;
  double err;
  double normsq;

  err=pow(10,-12);
 
  for(i=0;i<numpts;i++){
    ptgrad(pts,pts[i],ptsnbhd[i], s, radius, &gradarray[i*dim],dim);
    
    normsq=0;
    for(j=0;j<dim;j++){
	normsq+=pts[i].components[j]*pts[i].components[j];
	}
	
    if(normsq>1-err){
    		graddotpt=0;
    		for(j=0;j<dim;j++){
      		graddotpt+=gradarray[i*dim+j]*pts[i].components[j];
    	}
    	for(j=0;j<dim;j++){
      		gradarray[i*dim+j]+=-1.*graddotpt*pts[i].components[j];
    	}
    }
  }
}

//finds the gradient for all points in the shell
void gradShell( pt *pts, int numpts, nbhd *ptsnbhd, double *gradarray, double s, double radius,int dim,double inner_rad){
  int i,j ;
  double graddotpt;
  double err;
  double normsq;

  err=pow(10,-12);
 
  for(i=0;i<numpts;i++){
    ptgrad(pts,pts[i],ptsnbhd[i], s, radius, &gradarray[i*dim],dim);
    
    normsq=0;
    for(j=0;j<dim;j++){
	normsq+=pts[i].components[j]*pts[i].components[j];
	}
	
    if(normsq>1-err){
    	graddotpt=0;
    	for(j=0;j<dim;j++){
      	graddotpt+=gradarray[i*dim+j]*pts[i].components[j];
    	}
    	for(j=0;j<dim;j++){
      	gradarray[i*dim+j]+=-1.*graddotpt*pts[i].components[j];
    	}
    }
		else if(normsq<inner_rad*inner_rad+err){
			for(j=0;j<dim;j++){
				graddotpt+=gradarray[i*dim+j]*pts[i].components[j];
			}
			for(j=0;j<dim;j++){
				gradarray[i*dim+j]+=-1.*graddotpt*pts[i].components[j]/(inner_rad*inner_rad);
			}
		}
  }
}

//finds one random point in the ball
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

//finds one random point in the shell
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

//finds one random point in the sphere
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

//finds all random points in the ball
void randptsBall(double* coords, int numPts,int dim){
	int i;
	for(i = 0; i < numPts; i++){
    randptBall(&coords[i*dim], dim);
  }
}

//finds all random points in the shell
void randptsShell(double* coords, int numPts,int dim, double inner_radius){
  int i;
  for(i = 0; i < numPts; i++){
    randptShell(&coords[i*dim],dim, inner_radius);
  }
}

//finds all random points in the sphere
void randptsSphere(double* coords, int numPts, int dim){
  int i;
  for(i = 0; i < numPts; i++){
		randptSphere(&coords[i*dim],dim);
  }
}

//finds the energy of a set by creating the cube and neighborhood structure, and then finding the energy
double sEnergy(pt *points, double s, double radius, double* cornera, double* cornerb, int numpts, cube *cubes, int numCubes, nbhd *nbhds, int *neighborArray, cube *neighborhood, int maxNbrs, int *trinaryList, nbhd* cubeNbhd, int* testIndex, int* cubePtsArray,int dim){

  AssignCubes(radius, cornera, cornerb, points, numpts, cubes, cubePtsArray, maxNbrs,dim);
  AssignNbhds(cubes, radius, numpts, numCubes, points, nbhds, neighborArray, neighborhood, maxNbrs, trinaryList, cubeNbhd, testIndex,dim);
  return energy(points, numpts, nbhds, s, radius,dim);
}

//finds the energy and gradient of the ball by initializing the structure, then finding the gradient and energy
double sEnGradBall(pt *points, double *theGrad, double s, double radius, double* cornera, double* cornerb, int numpts, cube *cubes, int numCubes, nbhd *nbhds, int *neighborArray, cube *neighborhood, int maxNbrs, int *trinaryList, nbhd* cubeNbhd, int* testIndex, int* cubePtsArray,int dim, double inner_rad){
    
  AssignCubes(radius, cornera, cornerb, points, numpts, cubes, cubePtsArray, maxNbrs,dim);
  AssignNbhds(cubes, radius, numpts, numCubes, points, nbhds, neighborArray, neighborhood, maxNbrs, trinaryList, cubeNbhd, testIndex,dim);
  gradBall(points, numpts, nbhds, theGrad, s, radius,dim);
  return energy(points, numpts, nbhds, s, radius,dim);
}

//finds the energy and gradient of the sphere by initializing the structure, then finding the gradient and energy
double sEnGradSphere(pt *points, double *theGrad, double s, double radius, double cornera[3], double cornerb[3], int numpts, cube *cubes, int numCubes, nbhd *nbhds, int *neighborArray, cube *neighborhood, int maxNbrs, int *trinaryList, nbhd* cubeNbhd, int* testIndex, int* cubePtsArray,int dim, double inner_rad){
  AssignCubes(radius, cornera, cornerb, points, numpts, cubes, cubePtsArray, maxNbrs,dim);
  AssignNbhds(cubes, radius, numpts, numCubes, points, nbhds, neighborArray, neighborhood, maxNbrs, trinaryList, cubeNbhd, testIndex,dim);
  gradSphere(points, numpts, nbhds, theGrad, s, radius,dim);
  return energy(points, numpts, nbhds, s, radius,dim);
}

//finds the energy and gradient of the shell by initializing the structure, then finding the gradient and energy
double sEnGradShell(pt *points, double *theGrad, double s, double radius, double cornera[3], double cornerb[3], int numpts, cube *cubes, int numCubes, nbhd *nbhds, int *neighborArray, cube *neighborhood, int maxNbrs, int *trinaryList, nbhd* cubeNbhd, int* testIndex, int* cubePtsArray,int dim,double inner_rad){
  AssignCubes(radius, cornera, cornerb, points, numpts, cubes, cubePtsArray, maxNbrs,dim);
  AssignNbhds(cubes, radius, numpts, numCubes, points, nbhds, neighborArray, neighborhood, maxNbrs, trinaryList, cubeNbhd, testIndex,dim);
  gradShell(points, numpts, nbhds, theGrad, s, radius,dim,inner_rad);
  return energy(points, numpts, nbhds, s, radius,dim);
}

