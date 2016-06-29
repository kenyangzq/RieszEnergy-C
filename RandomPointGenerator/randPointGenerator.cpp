//
//  randPointGenerator.cpp
//  
//
//  Created by ken on 6/14/16.
//
//

#include <iomanip>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;


void randptSphere(double coordinates[], int dim){
    
    double z;
    double norm;
    double normsq=2;
    
    while(normsq>1 || normsq==0){
        normsq=0;
        
        for(int i=0;i<dim;i++){
            z=1-(2*(double)rand()/(double)RAND_MAX);
            normsq += z*z;
            coordinates[i] = z;
        }
    }
    
    norm=sqrt(normsq);
    
    for(int i=0;i<dim;i++){
        coordinates[i] = coordinates[i]/norm;
    }
    
}

int main(int argc, char* argv[]){
    
    int numpts = 1000;
    int dim = 3;
    int numfiles = 5;
    string filename = "randPt" + to_string(numpts)+ "-";
    
    double apoint[3];
    
    ofstream output_file;
    
    for (int k = 0; k < numfiles; k++) {
        
        string name = filename + to_string(k+1) + ".txt";
        
        output_file.open(name.c_str());
        if (output_file.fail()) {
            cout << "Error writing output data file" << endl;
        }
        
        output_file << setprecision(6);
        output_file << fixed;
        
        for (int i = 0; i < numpts; i++) {
            for (int j = 0; j < dim; j++) {
                randptSphere(apoint, dim);
                output_file << apoint[0] << "\t" << apoint[1] << "\t"
                << apoint[2] << "\n";
                
            }
        }
        
        output_file.close();

    }
    
    
    
    
    return 0;
}
