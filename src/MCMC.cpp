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
using namespace std;



/*
 *
 * This function runs totalS number of parallel chains each on a temperature level defined in temps.
 * Each MPI process is responsible for 1 or more chains in the pool.
 * Each chain is run iterNum number of iterations where each iteration consists of T number of steps
 * The chains communicates via MPI send and receive. temperedChains outputs result in the result array.
 *
 * Function Arguments:
 * iterNum: number of iterations
 * totalS: total number of parallel chains (each core may run more than one chain)
 * Nd: dimension of state space
 * T: number of steps each iteration
 * R: word pair counts from reference text
 * C: word pair counts from coded text
 * result: the pointer to array we output result
 * rank: rank of current MPI process
 * size: number of concurrent MPI processes
 *
 * */
void temperedChains(int iterNum, int totalS, int Nd, int T, int **R, int **C, double *temps, int * result, int rank, int size)
{
    double maxlogtarget=0.0;

    // Each MPI process is assigned S chains to run
    int S=totalS/size;

    if (rank==size-1)
    {
        S=totalS/size+totalS%size;
    }

    // Each MPI process has a different seed
    srand(unsigned(time(0))+rank);

    // Each chain creates a new starting state from uniform sampling
    int **xs;
    create2Dmemory(xs, S, Nd);

    int x0[Nd];
    for (int chains=0; chains<S; ++chains)
    {
        for (int i=0; i<Nd; ++i)
        {
            x0[i]=i;
        }
        rndpermutation(x0,Nd,x0);
        assignRow(xs, x0, Nd, chains);
    }


    /* Define variables used in the loop */
    int exchangetimes=0; // total number of exchange that occur
    double originallog; // store value of log target in prev step
    double proplog; // store value of log target of the proposal
    int glbc1; // global idx of chain 1 to exchange
    int glbc2; // global idx of chain 2 to exchange
    int rank1; // processor that runs chain 1
    int rank2; // processor that runs chain 2
    int c1; // local idx of chain 1 to exchange
    int c2; // local idx of chain 2 to exchange
    double coin; // store value of coin toss
    double accpt; // store value of acceptance ratio

    for (int iter=0; iter<iterNum; ++iter){

        for (int chains=0; chains<S; ++chains)
        {
            oneChain(xs[chains], T, Nd, xs[chains], R, C, temps[chains+rank*S]);
        }

        // Global index of the two chain to exchange positions
        glbc1=0;
        glbc2=((iter % totalS)==0) ? (1):(iter % totalS);
        if (totalS==1)
        {
            glbc2=0;
        }

        // Which processors these two indexes belong to
        rank1=glbc1/S;
        rank2=glbc2/S;

        // Current process is not involved in exchange
        if ((rank1!=rank)&&(rank2!=rank))
        {
            continue;
        }

        // Both indexes to exchange belong to current process
        else if (rank1==rank2){

            // Convert global indexes to local indexes
            c1=glbc1%S;
            c2=glbc2%S;

            // Compute acceptance ratio
            originallog=logtarget(xs[c1],Nd,R,C,temps[c1])+logtarget(xs[c2],Nd,R,C,temps[c2]);
            proplog=logtarget(xs[c1],Nd,R,C,temps[c2])+logtarget(xs[c2],Nd,R,C,temps[c1]);

            accpt=exp(proplog-originallog);

            // Determine if accept by toss a random coin
            coin=unifrnd(0,1);
            if (coin<accpt){
                // We indeed accept the proposal
                int *imp;
                imp=xs[c1];
                xs[c1]=xs[c2];
                xs[c2]=imp;
                exchangetimes+=1;
            }
        }

        // Indexes to exchange occur on two different processes
        else {

            if (rank == rank1){

                int c1=glbc1%S;
                double logtgtc1c1=logtarget(xs[c1],Nd,R,C,temps[glbc1]);
                double logtgtc1c2=logtarget(xs[c1],Nd,R,C,temps[glbc2]);

                int xsc2[Nd];
                MPI_Recv(xsc2,Nd,MPI_INT,rank2,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                double logtgtc2c2=logtarget(xsc2,Nd,R,C,temps[glbc2]);
                double logtgtc2c1=logtarget(xsc2,Nd,R,C,temps[glbc1]);


                originallog=logtgtc1c1+logtgtc2c2;
                proplog=logtgtc1c2+logtgtc2c1;

                accpt=exp(proplog-originallog);


                coin=unifrnd(0,1);
                if (coin<accpt){
                    // We indeed accept the proposal, notify the other chain
                    int accptstatus=1;

                    MPI_Send(&accptstatus,1,MPI_INT,rank2,0,MPI_COMM_WORLD);

                    MPI_Send(xs[c1],Nd,MPI_INT,rank2,0,MPI_COMM_WORLD);

                    deepcopy1Darray(xsc2,xs[c1], Nd);

                    exchangetimes+=1;

                } else {

                    int accptstatus=0;
                    MPI_Send(&accptstatus,1,MPI_INT,rank2,0,MPI_COMM_WORLD);
                }

            } else {
                int c2=glbc2%S;
                int accptstatus;

                MPI_Send(xs[c2],Nd,MPI_INT,rank1,0,MPI_COMM_WORLD);
                MPI_Recv(&accptstatus,1,MPI_INT,rank1,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                if (accptstatus==1)
                {
                    int xsc1[Nd];
                    MPI_Recv(xsc1,Nd,MPI_INT,rank1,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    deepcopy1Darray(xsc1,xs[c2], Nd);
                    exchangetimes+=1;
                }
            }
        }

        // Keep track of the most likely state so far
        if (rank==0)
        {
            double logtargetnow=logtarget(xs[0],Nd,R,C,1);
            if (logtargetnow>maxlogtarget)
            {
                maxlogtarget=logtargetnow;
                deepcopy1Darray(xs[0],result,Nd);

            }
        }

    }

    // At this point, broadcast result from rank 0  process to all MPI processes
    MPI_Bcast(result, Nd, MPI_INT, 0, MPI_COMM_WORLD);

    free2Dmemory(xs, S, Nd);
}


/*
 * This function runs a single Markov chain started at x0 for T steps at temperature temp and
 * output the last step at xT.
 *
 * Function arguments:
 * x0: the starting state
 * T: number of steps
 * Nd: dimension of state space
 * xT: last step state
 * R: word pair counts from reference
 * C: word pair counts from  coded file
 * temp: temperature to use
 *
 * */
void oneChain(int *x0, int T, int Nd, int *xT, int **R, int **C, double temp) {

    // Store all steps of the chain
    int **samples;
    create2Dmemory(samples, T, Nd);

    // Initialize the Markov chain
    assignRow(samples,x0,Nd,0);

    // Run the Markov chain
    for(int t=1; t<T; t=t+1)
    {

        if ( t % 1000 == 0 )
        {
            //printf("%d\n",t);
        }

        // Previous state
        int *xtm1=samples[t-1];

        // Propose the next state
        int xstar[Nd];

        rndswap(xtm1,Nd,xstar);

        // Determine if accept by toss a random coin
        double coin=unifrnd(0,1);

        // Compute acceptance ratio
        double accpt=exp(logtarget(xstar,Nd,R,C,temp)-logtarget(xtm1,Nd,R,C, temp));

        if (coin<accpt)
        {
            // We indeed accept the proposal
            assignRow(samples,xstar,Nd,t);
        } else {
            assignRow(samples,xtm1,Nd,t);
        }

    }

    for (int i=0; i<Nd; ++i)
    {
        xT[i]=samples[T-1][i];
    }

    free2Dmemory(samples,T,Nd);
}




/*
 * The function receives an array of current state of Markov chain and returns tempered target function value at that state
 *
 * Function Arguments:
 * x: the current state of the chain
 * Nd: dimension of state space
 * R: word pair counts from reference
 * C: word pair counts from coded text
 * temp: temperature of current level
 * return: target function value at the input state
 *
 * */
double logtarget(int *x, int Nd, int **R, int **C, double temp)
{

    double sum=0;
    #pragma omp parallel
    {
        #pragma omp for reduction(+:sum)
        for (int i=0; i<Nd; ++i)
        {
            for (int j=0; j<Nd; ++j)
            {
                int ci=x[i];
                int cj=x[j];
                sum+=log(double(R[i][j]))*double(C[ci][cj]);
            }
        }

    }

    //std::string decipheredstring=buildDecipheredstring(g_cipheredstring, x);

    return temp*sum;
}



/*
 * Rotate out the worst U chains and re-start them at the best U chains
 * */
void rotateout(int **xs, int S, int Nd, int U, int **R, int **C, double temp)
{
    int targetvals[S];

    for (int chains=0; chains<S; ++chains)
    {
        targetvals[chains]=logtarget(xs[chains], Nd, R, C, temp);

    }

    //Assume A is a given vector with N elements
    vector<int> V(S);
    int x=0;
    iota(V.begin(),V.end(),x++); //Initializing
    sort( V.begin(),V.end(), [&](int i,int j){return targetvals[i]<targetvals[j];} );

    for (int u=0; u<U; ++u)
    {
        assignRow(xs,xs[V[S-u-1]],Nd,V[u]);
    }
}
