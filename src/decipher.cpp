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
#include "decipher.h"

int main(int argc, char** argv){

    int rank, size;

    /* Initialize MPI and get rank and size */
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    // Dimension of the key
    int Nd=95;

    // Total number of samples we want
    int T=250;

    // Pool of chains
    int totalS=size;

    // Total number of iterations
    int iterNum=150;

    // Specify temperature of each chain
    double lowtemp=0.1;
    double temps[totalS];
    // double increment=(totalS==1)? 0 : ((1-lowtemp)/(totalS-1));
    for (int i=0; i<totalS; ++i)
    {
        temps[i]=lowtemp;//1-increment*i;
    }
    temps[0]=1;

    // Generate the ciphered file
    if (rank==0)
    {
        // Randomly generate a cipher key
        int cipherkey[Nd];
        for (int i=0; i<Nd-1; ++i)
        {
            cipherkey[i]=Nd-2-i;
        }
        cipherkey[Nd-1]=Nd-1;
        rndpermutation(cipherkey,Nd,cipherkey);

        // Use the cipher key to cipher the original file at this path: file2cipher
        std::string file2cipher="../data/code.txt";
        std::string cipheredfile="../data/ciphered.txt";
        buildCiphered(file2cipher, cipheredfile, cipherkey);
    }


    // Build transition matrix for reference text--count frequency of character pairs in reference text
    int **R;
    create2Dmemory(R, Nd, Nd);
    std::string referencetxt="../data/wap.txt";
    buildTransitionMat(R, Nd,referencetxt);


    // Count frequency of character pairs in the ciphered text (located at path: cipheredtxt)
    int **C;
    create2Dmemory(C, Nd, Nd);
    std::string cipheredtxt="../data/ciphered.txt";
    buildTransitionMat(C, Nd,cipheredtxt);

    
    // Decipher the text using api temperedChains and store output in [result] variable below
    int result[Nd];
    temperedChains(iterNum, totalS, Nd, T, R, C, temps, result, rank, size);

    if (rank==0)
    {
        // Print the result
        print1Darray(result,Nd);
        printf("Target:%f\n", logtarget(result,Nd,R,C,1));

        // Use the key found to decipher the ciphered text and store it at this path: decipheredtext
        std::string decipheredtext="../data/deciphered.txt";
        cipherkey2decipherkey(result,result,Nd);
        buildDeciphered(cipheredtxt,decipheredtext,result);
    }

    free2Dmemory(C, Nd, Nd);
    free2Dmemory(R, Nd, Nd);
    MPI_Finalize();
    return 0;
}


