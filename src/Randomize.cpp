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




/*
 * This function receives an array of input of length N and uniformly select k elements (no replacement) and fill output with it
 * */
void rndkofn(int *input, int *output, int n, int k)
{

    int permuted[n];
    rndpermutation(input,n,permuted);

    for (int i=0; i<k; ++i)
    {
        output[i]=permuted[i];
    }

}


/*
 * The function receives an array of current state of Markov chain and fills the proposed state
 *
 * x: the current state of the chain
 * x_swapped: the proposed state of the chain
 *
 * */

int rndswap(int *x, int x_len, int *x_swapped)
{

    int idxarr[x_len];

    for (int j=0; j<x_len; ++j)
    {
        idxarr[j]=j;
    }

    // Permute the idxarr
    int idx_pered[x_len];
    rndpermutation(idxarr, x_len, idx_pered);
    int swap1=idx_pered[0];
    int swap2=idx_pered[1];


    for (int j=0; j<x_len; ++j)
    {
        if (j==swap1)
        {
            x_swapped[j]=x[swap2];
        }
        else if (j==swap2)
        {
            x_swapped[j]=x[swap1];
        }
        else
        {
            x_swapped[j]=x[j];
        }
    }

    return 0;

}

int rndswapwtidx(int *x, int x_len, int *x_swapped, int &swap1, int &swap2)
{

    int idxarr[x_len];

    for (int j=0; j<x_len; ++j)
    {
        idxarr[j]=j;
    }

    // Permute the idxarr
    int idx_pered[x_len];
    rndpermutation(idxarr, x_len, idx_pered);
    swap1=idx_pered[0];
    swap2=idx_pered[1];


    for (int j=0; j<x_len; ++j)
    {
        if (j==swap1)
        {
            x_swapped[j]=x[swap2];
        }
        else if (j==swap2)
        {
            x_swapped[j]=x[swap1];
        }
        else
        {
            x_swapped[j]=x[j];
        }
    }

    return 0;

}



/*
 * The function returns a real number uniformly sampled between a,b
 *
 * a: left boundary of the interval
 * b: right boundary of the interval
 * return: a real number uniformly sampled between a,b
 *
 * */

double unifrnd(double a, double b)
{
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(a,b);
    return dis(gen);
}

/*
 * The function returns an integer uniformly sampled between a,b inclusive
 *
 * a: left boundary of the interval
 * b: right boundary of the interval
 * return: an integer uniformly sampled between a,b
 *
 * */

int unifrndint(int a, int b)
{
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dis(a,b);
    return dis(gen);
}


/*
 * The function returns an integer uniformly sampled between a,b inclusive
 *
 * a: left boundary of the interval
 * b: right boundary of the interval
 * return: an integer uniformly sampled between a,b
 *
 * */

// random generator function
int rndpermutation(int *inputarr, int inputarr_len, int *outputarr)
{
    std::vector<int> arr;

    // set some values:
    for (int j = 0; j < inputarr_len; ++j)
        arr.push_back(inputarr[j]);

    std::random_device rd;
    auto rng = std::default_random_engine {rd()};
    std::shuffle(std::begin(arr), std::end(arr), rng);


    for (int j = 0; j < inputarr_len; ++j)
        outputarr[j]=arr[j];

    return 0;
}


/*
 * The function returns a real number uniformly sampled between a,b
 *
 * a: left boundary of the interval
 * b: right boundary of the interval
 * return: a real number uniformly sampled between a,b
 *
 * */

double normrnd(double mu, double std)
{
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::normal_distribution<double> dis(mu,std);
    return dis(gen);
}