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

    /*int **x;
    int Nd=32;
    create2Dmemory(x,Nd,Nd);
    for (int i=0; i<Nd; ++i)
    {
        for (int j=0;j<Nd; ++j)
        {
            if (unifrnd(0,1)<0.5)
                x[i][j]=1;
            else{
                x[i][j]=0;
            }
        }
    }
    oneChainIsing(x, 10,Nd, 0.38);
    print2Darray(x,Nd,Nd);
*/


    int rank, size;

    /* Initialize MPI and get rank and size */
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Dimension of the key
    int Nd=5;

    // Number of steps each iteration
    int T=500;

    // Pool of chains
    int totalS=size;

    // Total number of iterations
    int iterNum=100;

    // Specify temperature of each chain
    double lowtemp=1;
    double temps[totalS];
    // double increment=(totalS==1)? 0 : ((1-lowtemp)/(totalS-1));
    for (int i=0; i<totalS; ++i)
    {
        temps[i]=lowtemp;//1-increment*i;
    }
    temps[0]=1;


    // Decipher the text using api temperedChains and store output in [result] variable below
    int **result;
    create2Dmemory(result,Nd,Nd);
    temperedChainsIsing(iterNum, totalS, Nd, T, temps, result, rank, size);

    // Print the result



    free2Dmemory(result, Nd, Nd);
    MPI_Finalize();
    return 0;
}


