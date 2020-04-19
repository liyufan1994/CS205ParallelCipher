//
// Created by Yufan Li on 2020-04-01.
//

#include "Ising.h"


#include <iostream>
#include <math.h>
#include <random>
#include <vector>
#include <fstream>
#include <map>
#include <string>
#include <algorithm>
#include <mpi.h>
#include <omp.h>
#include <chrono>

#include "Ising.h"
#include "ArrayUtilities.h"
#include "Randomize.h"




int main(int argc, char** argv){

    // The Ising lattice will be of size Nd x Nd
    int Nd = (argc > 1) ? atoi(argv[1]) : 32;

    // Each MPI will run S chains
    int S = (argc > 2) ? atoi(argv[2]) : 1;    

    // What decomposition to use: 0-strip, 1-checkerboard
    int D = (argc > 3) ? atoi(argv[3]) : 1; 

    
    // ID of MPI process and number of MPI processes respectively
    int rank, size;

    // Initialize MPI and get rank and size
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Number of steps each iteration
    int T=1;

    // Pool of chains
    int totalS=size*S;

    // Total number of iterations
    int iterNum=100;

    // Specify temperature of each chain
    // If totalS=1, the sole chain will only simulate the high temperature
    double lowtemp=0.3;
    double hightemp=0.6;

    double temps[totalS];
    double increment=(totalS==1)? 0 : ((hightemp-lowtemp)/(totalS-1));

    for (int i=0; i<totalS; ++i)
    {
        temps[i]=hightemp-increment*i;
    }
    temperedChainsIsing(iterNum, totalS, Nd, T, temps,  rank, size,D);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}


