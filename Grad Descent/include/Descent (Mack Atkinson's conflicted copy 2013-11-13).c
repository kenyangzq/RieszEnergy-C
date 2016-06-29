//This is now general for all sets

#include "Energy.h"

void updateBall(pt *newpts, pt *oldpts, int numpts, double *pvec, double t,int dim,double inner){
  int i, j;
  double distSQ, newcomponent, scalefac;
  for(i=0;i<numpts;i++){
    newpts[i].index=oldpts[i].index;
    distSQ=0;
    for(j=0;j<dim;j++){
      newcomponent=oldpts[i].components[j]-t*pvec[i*dim+j];
      newpts[i].components[j]=newcomponent;
      distSQ+=newcomponent*newcomponent;
    }
    if (distSQ>1){
      scalefac=1/sqrt(distSQ);     
      for(j=0;j<dim;j++){
        newpts[i].components[j]=scalefac*newpts[i].components[j];
      }
    }
  }     
}

void updateShell(pt *newpts, pt *oldpts, int numpts, double *pvec, double t,int dim,double inner){
    int i, j;
    double distSQ, newcomponent, scalefac;
    
    for(i=0;i<numpts;i++){
        newpts[i].index=oldpts[i].index;
        distSQ=0;
        for(j=0;j<dim;j++){
            newcomponent=oldpts[i].components[j]-t*pvec[i*dim+j];
            newpts[i].components[j]=newcomponent;
            distSQ+=newcomponent*newcomponent;
        }
        if (distSQ>1){
            scalefac=1/sqrt(distSQ);     
            for(j=0;j<dim;j++){
                newpts[i].components[j]=scalefac*newpts[i].components[j];
            }
        }
        if (distSQ<inner){
						//not sure what this scale factor does
            scalefac=.55/sqrt(distSQ);
            for(j=0;j<dim;j++){
                newpts[i].components[j]=scalefac*newpts[i].components[j];
            }
        }
    }     
}

void updateSphere(pt *newpts, pt *oldpts, int numpts, double *pvec, double t,int dim,double inner){
    int i, j;
    double distSQ, newcomponent, scalefac;
    
    for(i=0;i<numpts;i++){
        newpts[i].index=oldpts[i].index;
        distSQ=0;
        for(j=0;j<dim;j++){
            newcomponent=oldpts[i].components[j]-t*pvec[i*dim+j];
            newpts[i].components[j]=newcomponent;
            distSQ+=newcomponent*newcomponent;
        }
        scalefac=1/sqrt(distSQ);     
        for(j=0;j<dim;j++){
            newpts[i].components[j]=scalefac*newpts[i].components[j];
        }
    }     
}

double GradDescent(pt *points, pt *newpoints, int numpts, double *theGrad, double t0, double radius, double* cornera, double* cornerb, cube *cubes, int numCubes, int *neighborArray, int maxNbrs, double s, cube *neighborhood, nbhd *nbhds, FILE *writeTo, int *trinaryList, nbhd* cubeNbhd, int* testIndex, int* cubePtsArray,int dim,double inner,void(*updateMethod)(pt *,pt *,int,double *,double,int),double(*gradMethod)(pt *,double *,double,double,double *,double *,int,cube *,int,nbhd *,int *,cube *,int,int *,nbhd *,int *,int *,int)){
  double en0, en1, en2;
  int count = 0; 
  double t1,t_cal;

  en0=gradMethod(points, theGrad, s, radius, cornera, cornerb, numpts, cubes, numCubes, nbhds, neighborArray, neighborhood, maxNbrs, trinaryList, cubeNbhd, testIndex,cubePtsArray,dim);

////////////////////////NEW CODE STARTS HERE//////////////////////////////////////////////
	//I am inputting a method that finds t0 by finding the min distance and the max gradient
	double t00;
	int i,j,d;
	double delta=100000.0;
	double maxV=0.0;
	double p_dist;
	double g_dist;
	for(i=0;i<numpts;i++){
		for(j=i+1;j<numpts;j++){
			p_dist=dist(points[i],points[j],dim);
			g_dist=0.0;
			for(d=0;d<dim;d++){
				g_dist+=(pow((theGrad[i*dim+d]-theGrad[j*dim+d]),2));
			}
			g_dist=pow(g_dist,0.5);
			if(p_dist<delta){
				delta=p_dist;
			}
			if(g_dist>maxV){
				maxV=g_dist;
			}
		}
	}
	t00=0.5*delta/maxV;
/*	printf("Old t0 = %lf\n",t0);*/
/*	printf("New t0 = %lf\n",t00);*/
	t0=t00;

////////////////////////////END OF NEW SEGMENT OF CODE//////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

	t1=t0;
  updateMethod(newpoints, points, numpts, theGrad, t0,dim);
  en1 = sEnergy(newpoints, s, radius, cornera, cornerb, numpts, cubes, numCubes, nbhds, neighborArray, neighborhood, maxNbrs, trinaryList, cubeNbhd, testIndex, cubePtsArray,dim);
	printf("en0 = %lf\n",en0);
  while(en1 > en0 && count < 10){
    t1 = 0.1* t1;
    updateMethod(newpoints, points, numpts, theGrad, t1,dim);
    en1 = sEnergy(newpoints, s, radius, cornera, cornerb, numpts, cubes, numCubes, nbhds, neighborArray, neighborhood, maxNbrs, trinaryList, cubeNbhd, testIndex,cubePtsArray,dim);
    count++;
  }   
  if(en1>en0){
    perror("choose smaller t0");
    return 0;
  }
	//I don't understand why this is included
  if(t1 != t0){
    updateMethod(newpoints, points, numpts, theGrad, t1,dim);
    en1 = sEnergy(newpoints, s, radius, cornera, cornerb, numpts, cubes, numCubes, nbhds, neighborArray, neighborhood, maxNbrs, trinaryList, cubeNbhd, testIndex,cubePtsArray,dim);
  }
  updateMethod(newpoints, points, numpts, theGrad, 2*t1,dim);
  en2 = sEnergy(newpoints, s, radius, cornera, cornerb, numpts, cubes, numCubes, nbhds, neighborArray, neighborhood, maxNbrs, trinaryList, cubeNbhd, testIndex,cubePtsArray,dim);
  int k = 1;
  while(en2 < en1 && k < 10){
    en0 = en1;
    en1 = en2;
    updateMethod(newpoints, points, numpts, theGrad, t1 * (k+2),dim);
    en2 = sEnergy(newpoints, s, radius, cornera, cornerb, numpts, cubes, numCubes, nbhds, neighborArray, neighborhood, maxNbrs, trinaryList, cubeNbhd, testIndex,cubePtsArray,dim);
    k++;
  }
	//is this correct?? Should it be a while loop?
	if(en2>en1){
    t_cal = t1*(k + (en0 - en2)/(2*(en0+en2-2*en1)));
    updateMethod(newpoints, points, numpts, theGrad, t_cal,dim);
		en2 = sEnergy(newpoints, s, radius, cornera, cornerb, numpts, cubes, numCubes, nbhds, neighborArray, neighborhood, maxNbrs, trinaryList, cubeNbhd, testIndex,cubePtsArray,dim);
		if(en2>en1){
		//this is one step between the high en2 and low en1, its not (k+1) because of how k increments
			t_cal=t1*k;
			updateMethod(newpoints,points,numpts,theGrad,t_cal,dim);
			en2 = sEnergy(newpoints, s, radius, cornera, cornerb, numpts, cubes, numCubes, nbhds, neighborArray, neighborhood, maxNbrs, trinaryList, cubeNbhd, testIndex,cubePtsArray,dim);
		}
  }
	else{
		t_cal=(k+1)*t1;
	}

/////////////////NEW CODE THAT PRINTS OUT IF t_n*v_n///////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

	FILE *epsFile = fopen("epsFile.txt","a+");
	if(roundcount==0){
		fprintf(epsFile,"n\tt\t\tmaxV\t\tt*maxV\t\tenergy\n");
	}
	maxV=0;
	for(i=0;i<numpts;i++){
		for(j=i+1;j<numpts;j++){
			p_dist=dist(newpoints[i],newpoints[j],dim);
			g_dist=0.0;
			for(d=0;d<dim;d++){
				g_dist+=(pow((theGrad[i*dim+d]-theGrad[j*dim+d]),2));
			}
			g_dist=pow(g_dist,0.5);
			if(g_dist>maxV){
				maxV=g_dist;
			}
		}
	}
	if(roundcount%5==0){
		fprintf(epsFile,"%d\t%lf\t%lf\t%lf\t%lf\n",roundcount,t_cal,maxV,t_cal*maxV,en2);
	}	
	fclose(epsFile);
	
////////////////////////////END OF NEW CODE//////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	//why is t_calE14 output?
  fprintf(writeTo, "%d\t%lf\t%lf\t", roundcount, t_cal*pow(10, 14), en1);
  fprintf(writeTo, "\t%lf\n", en2);
  roundcount++;
  return t_cal;
}

