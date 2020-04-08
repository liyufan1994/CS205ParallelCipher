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

    // ID of MPI process and number of MPI processes respectively
    int rank, size;

    // Initialize MPI and get rank and size
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Number of steps each iteration
    int T=5;

    // Pool of chains
    int totalS=size*2;

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
    temperedChainsIsing(iterNum, totalS, Nd, T, temps,  rank, size,1);

    MPI_Finalize();
    return 0;
}


