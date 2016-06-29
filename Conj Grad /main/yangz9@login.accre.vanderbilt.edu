#include <stdio.h>

//This file works for any dimension now

//currents trying to fix the point problem: making it dynamic

#include "Descent.h"

int main(int argc, char* argv[]){
  srand(time(NULL));

//creating method pointers
//these pointers allow for the set to be chosen at control.inp
	void (*update)(pt *,pt *,int,double *,double,int,double);
	double (*sEnGrad)(pt *points, double *theGrad, double s, double radius, double cornera[3], double cornerb[3], int numpts, cube *cubes, int numCubes, nbhd *nbhds, int *neighborArray, cube *neighborhood, int maxNbrs, int *trinaryList, nbhd* cubeNbhd, int* testIndex, int* cubePtsArray,int dim, double inner);

//declaring variables to read from control.inp
	int mani;
	int maxNbrs;
	double inner;
	double c;
	double s;
	int dim;
	int iterations;
	int quest=0;
	int partitionNum;
	int numpts;
	char fullName[80];
	char name[20];
	char inFile[20];
	char type[20];

    
    
    char zeros[1024], zeros_var[1024];
    memset(zeros, '0', sizeof zeros);
    int numAddedZeros = 0;

//reading control.inp
  FILE *input= fopen("control.inp", "r");
	char stuff[20];
	char eq;
	int val; 
	int count=1;
  if (input==NULL) 
	  perror ("Error opening file");    
  else {
	  int i;
		char line[200];
		fgets ( line, sizeof line, input );
		fgets ( line, sizeof line, input );
		fscanf(input, "%s",type);
    while (fscanf(input, "%s", stuff) != EOF) {
			count++;
			if(quest==1&&count==6){
				fscanf(input,"%c",&eq);
				fscanf(input,"%c",&eq);
				fscanf(input,"%s",inFile);
			}	
			if(count==4){
				fscanf(input,"%c",&eq);
				fscanf(input,"%c",&eq);
				fscanf(input,"%lf",&c);
			}
			if(count==10){
				fscanf(input,"%c",&eq);
				fscanf(input,"%c",&eq);
				fscanf(input,"%s",name);
			}		
			if(count==11){
				fscanf(input,"%c",&eq);
				fscanf(input,"%c",&eq);
				fscanf(input,"%lf",&inner);
			}
        if(count==2){
            fscanf(input,"%c",&eq);
            fscanf(input,"%c",&eq);
            fscanf(input,"%lf",&s);
        }
			else{
				fscanf(input,"%c",&eq);
				fscanf(input,"%c",&eq);
				fscanf(input,"%d",&val);
				if(count==3){
					dim=val;
				}
				if(count==5){
					quest=val;
				}
				if(quest==0&&count==6){
					numpts=val;
				}
				if(count==7){
					iterations=val;
				}
				if(count==8){
					partitionNum=val;
				}	
				if(count==9){
					maxNbrs=val;
				}														
			}
		}
	}
/*	printf("%s\n",inFile);*/

	if(type[1]=='a'){
		update=updateBall;
		sEnGrad=sEnGradBall;
		mani=dim;
	}
	else if(type[1]=='p'){
		update=updateSphere;
		sEnGrad=sEnGradSphere;
		mani=dim-1;
	}
	else if(type[1]=='h'){
		update=updateShell;
		//this is because there is no difference in the gradients for ball and shell
		sEnGrad=sEnGradShell;
		mani=dim;
	}

//creating initial points either randomly, or from an input file
  int numFilesWritten = 1;
  int i;
	int j;
	int k;



//counting numpts in the input file
	if(quest==1){
		FILE *inPoints= fopen(inFile, "r");
		if (inPoints==NULL) 
      perror ("Error opening file");    
    else {
			numpts=0;
			while(fscanf(inPoints, "%s", stuff) != EOF){
				numpts++;
			}
			numpts=numpts/dim;
			fclose(inPoints);
		}	
	}

	pt* points = malloc(numpts * sizeof(pt));
  double* coords = malloc(numpts * dim * sizeof(double));

//if questo==0, then must create random points
//checking the second letter of the name of the set in order to determine which set it is
  if(quest==0){
		if(type[1]=='a'){
    	randptsBall(coords,numpts,dim);
		}
		else if(type[1]=='p'){
			randptsSphere(coords,numpts,dim);
		}
		else if(type[1]=='h'){
			randptsShell(coords,numpts,dim,inner);
		}					
		for(i=0;i<numpts;i++){
			points[i].index=i;
			points[i].components=&coords[i*dim];
		}
  }
//otherwise, read the input file
	else{
		double inp;
		FILE *inPoints = fopen(inFile, "r");
    int i;
    for(i = 0; i < numpts*dim; i++){
      fscanf(inPoints, "%lf", &inp);
			coords[i]=inp;
    }
    fclose(inPoints);
    for(i=0;i<numpts;i++){
	  	points[i].index=i;
			points[i].components=&coords[i*dim];
		}
  }

//opening the output file to write points to
	sprintf(fullName, "%s.txt", name); 
	FILE *writeTo = fopen(fullName, "a+");
	sprintf(fullName, "%s.txt", name);
	FILE *writePtsTo = fopen(fullName, "a+");	
  for(j = 0; j < numpts; j++){
		for(k=0;k<dim;k++){		
			fprintf(writePtsTo,"%lf\t",points[j].components[k]);
		}
		fprintf(writePtsTo,"\n");
  }
	fclose(writePtsTo);
  sprintf(fullName, "%s%d.txt", name, 1);
  writePtsTo = fopen(fullName, "a+");


//declaring all necessary variables to run gradient descent

  double radius = c*pow(numpts,-1.0/mani);
	double* cornera = malloc(dim*sizeof(double));
	double* cornerb = malloc(dim*sizeof(double));
	for(j=0;j<dim;j++){
		cornera[j]=-1;
		cornerb[j]=1;
	}
  int numCubes= findNumCubes(radius, cornera, cornerb);    
  cube* neighborhood = malloc(pow(3, dim) * sizeof(cube));
  nbhd* nbhds = malloc( numpts * sizeof(nbhd));
  cube* cubes = malloc(pow(numCubes, dim) * sizeof(cube)); 
  double *newGrad=malloc(numpts*dim*sizeof(double));
  double *oldGrad=malloc(numpts*dim*sizeof(double));
  double *theDir=malloc(numpts*dim*sizeof(double));
  int* trinaryList = malloc(pow(3, dim) * dim * sizeof(int));
  int* testIndex = malloc(dim * sizeof(int));
  int* neighborArray = malloc(maxNbrs*numpts*sizeof(int));
  int* cubePtsArray = malloc(pow(numCubes,dim) * maxNbrs * sizeof(int));
  int* cubeNbhdArray = malloc(pow(3, dim) * maxNbrs*2* sizeof(int));
  nbhd cubeNbhd;
  cubeNbhd.neighbors = cubeNbhdArray;
  if (cubeNbhdArray==NULL || cubePtsArray==NULL) {
    printf("Not enough memory\n");
    return 0;
  }
  int fclose(FILE *a_file);
  pt* newpoints = malloc(numpts * sizeof(pt));
	double* newcoords = malloc(numpts * dim * sizeof(double));
	for(i = 0; i < numpts; i++){
    newpoints[i].index = i;
		newpoints[i].components=&newcoords[i*dim];
    }
  double t0 = 1.1*pow(10, -12);
  time_t ab, cd;
  (void) time(&ab);
	double delta;
	int bins=20;
	int * hist = malloc(bins*sizeof(int));
	delta=radius/bins;



//main part of the program, running conjugate gradient iterations
double currEnergy;
currEnergy=sEnGrad(points, oldGrad, s, radius, cornera, cornerb, numpts, cubes, numCubes, nbhds, neighborArray, neighborhood, maxNbrs, trinaryList, &cubeNbhd,testIndex,cubePtsArray,dim,inner);

for(i=0;i<numpts*dim;i++){
	theDir[i]=oldGrad[i];
}


  for (i = 0; i <= iterations; i++){
    if(i == iterations/(partitionNum)*numFilesWritten){
      numFilesWritten++;
      for(j = 0; j < numpts; j++){
				for(k=0;k<dim;k++){		
					fprintf(writePtsTo,"%lf\t",newpoints[j].components[k]);
				}
				fprintf(writePtsTo,"\n");
      }
      fclose(writePtsTo);

      /*sprintf(fullName, "%s%s%d.txt", name, "nbhds", numFilesWritten-1);*/
      /*writePtsTo = fopen(fullName, "a+");*/
			/*FILE *writeHist = fopen("hist.txt", "w");*/
///ADDING CODE TO TEST HOW MANY POINTS NEED TO BE ASSIGNED TO NBHDS//////////////
			/*sprintf(fullName, "%s%s%d.txt", name, "nbhdcounts", numFilesWritten-1);*/
			/*FILE *counts = fopen(fullName,"a+");*/
			/*int var;*/
			/*int maxNb=0;*/
			/*double avgNb=0;*/
      /*for(j = 0; j < numpts; j++){*/
        /*fprintf(ritePtsTo, "\n%d", j);*/
        /*[>for(k = 0; k < nbhds[j].neighborcount; k++){<]*/
          /*[>fprintf(writePtsTo, " %d \t", nbhds[j].neighbors[k]);<]*/
					/*[>//inputting code to write to histogram <]*/
					/*[>//this actually counts it to make a graph<]*/
					/*[>fprintf(writeHist,"%d\n",(int)(floor(dist(points[j],points[nbhds[j].neighbors[k]],dim)/delta)));<]*/
        /*[>}<]*/
	/*///added for experimenting////////////////*/
				/*var=nbhds[j].neighborcount;*/
				/*fprintf(counts,"%d \n",var);*/
				/*if(var>maxNb){*/
					/*maxNb=var;*/
				/*}*/
				/*avgNb=avgNb+var;*/
				/*fprintf(writePtsTo,"\n");*/
      /*}*/
			/*avgNb=avgNb/numpts;*/
			/*fprintf(counts,"maxNb = %d\n",maxNb);*/
			/*fprintf(counts,"avgNb = %lf",avgNb);*/
			/*fclose(counts);*/
      /*fclose(writePtsTo);	*/
			/*fclose(writeHist);*/

/////THIS IS ADDED CODE TO TEST HOW MANY POINTS NEED TO BE ASSIGNED TO EACH CUBE///////////
            /*sprintf(fullName, "%s%s%d.txt", name, "cubes", numFilesWritten-1);*/
			/*writePtsTo=fopen(fullName,"a+");*/
			/*int maxCube=0;*/
			/*int cubeCount=0;*/
			/*//avgCube finds average points per cube where denominator is all cubes*/
			/*double avgCube=0;*/
			/*//avgCubeCon finds average points per cube where denominator is only the cubes containing points*/
			/*double avgCubeCon=0;*/
			/*for(j=0;j<pow(numCubes,dim);j++){*/
				/*var=cubes[j].pointCount;*/
				/*fprintf(writePtsTo,"%d \n",var);*/
				/*if(var>maxCube){*/
					/*maxCube=var;*/
				/*}*/
				/*if(var>0){*/
					/*avgCube=avgCube+var;*/
					/*cubeCount++;*/
				/*}*/
			/*}*/
			/*avgCubeCon=avgCube/cubeCount;*/
			/*avgCube=avgCube/pow(numCubes,dim);*/
			/*fprintf(writePtsTo,"maxCube = %d\n",maxCube);*/
			/*fprintf(writePtsTo,"avgCubeCon = %lf\n",avgCubeCon);*/
			/*fprintf(writePtsTo,"avgCube = %lf\n",avgCube);*/
			/*fclose(writePtsTo);*/
///////////END OF TEST CODE////////////////////////////////////////////////////////////
        numAddedZeros = (int) (ceil(log10(iterations+1)) - ceil(log10(numFilesWritten+1)));
        strncpy(zeros_var, zeros, numAddedZeros);
        zeros_var[numAddedZeros] = '\0';
      sprintf(fullName, "%s%d.txt", name,  numFilesWritten);
      writePtsTo = fopen(fullName, "a+");
    }


//time consuming part
    t0 = ConjGrad(points, newpoints, numpts, currEnergy, theDir, newGrad, oldGrad, t0, radius, cornera, cornerb, cubes, numCubes, neighborArray, maxNbrs, s, neighborhood, nbhds, writeTo, trinaryList, &cubeNbhd, testIndex,cubePtsArray,dim,inner,update,sEnGrad);

		//determines if v*t0 has reached a low enough value to stop
		if(done==1){
			printf("\n done \n");
			break;
		}

currEnergy=sEnergy(newpoints, s, radius, cornera, cornerb, numpts, cubes, numCubes, nbhds, neighborArray, neighborhood, maxNbrs, trinaryList, &cubeNbhd, testIndex, cubePtsArray,dim);

    pt *swap = newpoints;
    newpoints = points;
    points = swap;

    double *swapGrad = newGrad;
    newGrad = oldGrad;
    oldGrad = swapGrad;

  }

//writing the final points   
  (void) time(&cd);
	char finalName[80];
	sprintf(finalName, "%sfinal.txt", name);
	writePtsTo=fopen(finalName,"a+");
  for(j = 0; j < numpts; j++){
		for(k=0;k<dim;k++){		
			fprintf(writePtsTo,"%lf\t",newpoints[j].components[k]);
		}
		fprintf(writePtsTo,"\n");
  }
  sprintf(finalName, "%s%sfinal.txt", name, "nbhds");
  writePtsTo = fopen(finalName, "a+");
  for(j = 0; j < numpts; j++){
    fprintf(writePtsTo, "\n%d", j);
    for(k = 0; k < nbhds[j].neighborcount; k++){
      fprintf(writePtsTo, " %d \t", nbhds[j].neighbors[k]);
    }
  }
  //printf("took %ld seconds", cd-ab);
  fclose(writePtsTo);
  fclose (writeTo);
  return 0;
}
