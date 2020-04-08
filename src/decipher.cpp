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

#include "decipher.h"
#include "ArrayUtilities.h"




int main(int argc, char** argv){

    int rank, size;

    // Initialize MPI and get rank and size
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Dimension of the key
    int Nd=95;

    // Number of steps each iteration
    int T=500;

    // Pool of chains
    int totalS=size;

    // Total number of iterations
    int iterNum=100;

    // Specify temperature of each chain
    double lowtemp=100;
    double temps[totalS];

    for (int i=0; i<totalS; ++i)
    {
        temps[i]=lowtemp;//1-increment*i;
    }
    temps[0]=10000;

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
        //rndpermutation(cipherkey,Nd,cipherkey);

        // Use the cipher key to cipher the original file at this path: file2cipher
        std::string file2cipher="../data/code.txt";
        std::string cipheredfile="../data/ciphered.txt";
        buildCiphered(file2cipher, cipheredfile, cipherkey);
    }

    // Count frequency of character pairs in reference text
    int **R;
    create2Dmemory(R, Nd, Nd);
    std::string referencetxt="../data/ak.txt";
    buildTransitionMat(R, Nd,referencetxt);

    // Count frequency of character pairs in the ciphered text (located at path: cipheredtxt)
    int **C;
    create2Dmemory(C, Nd, Nd);
    std::string cipheredtxt="../data/ciphered.txt";
    buildTransitionMat(C, Nd,cipheredtxt);

    // Decipher the text using api temperedChains and store output in [result] variable below
    int result[Nd];
    temperedChains(iterNum, totalS, Nd, T, R, C, temps, result, rank, size);

    // Print the result
    if (rank==0) print1Darray(result, Nd);
    if (rank==0) printf("Target:%f\n", logtarget(result, Nd, R, C, 1));


    // Use the key found to decipher the ciphered text and store it at this path: decipheredtext
    std::string decipheredtext = "../data/deciphered.txt";
    cipherkey2decipherkey(result, result, Nd);
    buildDeciphered(cipheredtxt, decipheredtext, result);


    // Put deciphered file into a string for finer processing
    std::string g_cipheredstring=readcodefile(cipheredtxt);
    std::string decipheredstring=buildDecipheredstring(g_cipheredstring, result);
    if (rank==0) printf("%s\n",decipheredstring.c_str());

    std::map<std::string, int> g_dict;
    buildWordsFreqMap("../data/google-10000-english-usa.txt", g_dict);

    int CWFS;
    double percwords;
    int wordcount;
    CWFScore(decipheredstring, CWFS, percwords, g_dict,wordcount);

    int resultfine[Nd];
    fineOptimize(result,Nd,10000,g_cipheredstring,resultfine,g_dict,rank);
    decipheredstring=buildDecipheredstring(g_cipheredstring, resultfine);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==0) printf("%s\n",decipheredstring.c_str());


    free2Dmemory(C, Nd, Nd);
    free2Dmemory(R, Nd, Nd);
    MPI_Finalize();
    return 0;
}


