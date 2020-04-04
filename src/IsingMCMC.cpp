
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
#include "Randomize.h"
#include "Ising.h"
#include "ArrayUtilities.h"


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
* result: the pointer to array we output result, it consists of totalS number of columns and iterNum of rows
* rank: rank of current MPI process
* size: number of concurrent MPI processes
*
* */
void temperedChainsIsing(int iterNum, int totalS, int Nd, int T, double *temps, int ** result, int rank, int size)
{
    // Each MPI process is assigned S chains to run
    int S=totalS/size;

    if (rank==size-1)
    {
        S=totalS/size+totalS%size;
    }

    // Each MPI process stores results in partialresult array of S columns and iterNum rows
    int **partialresult;
    create2Dmemory(partialresult,iterNum,S);

    // Each MPI process has a different seed
    srand(unsigned(time(0))+rank);

    // Each chain creates a new starting state from uniform sampling
    int ***xs;
    create3Dmemory(xs, S, Nd,Nd);

    for (int chains=0; chains<S; ++chains)
    {
        for (int i=0; i<Nd; ++i)
        {
            for (int j=0;j<Nd; ++j)
            {
                if (unifrnd(0,1)<0.5)
                    xs[chains][i][j]=1;
                else{
                    xs[chains][i][j]=0;
                }
            }
        }
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
            oneChainIsing(xs[chains], T, Nd, temps[chains+rank*S]);
            partialresult[iter][chains]=t(xs[chains],Nd);
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
            originallog=logtargetIsing(xs[c1],Nd,temps[c1])+logtargetIsing(xs[c2],Nd,temps[c2]);
            proplog=logtargetIsing(xs[c1],Nd,temps[c2])+logtargetIsing(xs[c2],Nd,temps[c1]);

            accpt=exp(proplog-originallog);

            // Determine if accept by toss a random coin
            coin=unifrnd(0,1);
            if (coin<accpt){
                // We indeed accept the proposal
                int **imp;
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
                double logtgtc1c1=logtargetIsing(xs[c1],Nd,temps[glbc1]);
                double logtgtc1c2=logtargetIsing(xs[c1],Nd,temps[glbc2]);

                int **xsc2;
                create2Dmemory(xsc2,Nd,Nd);
                MPI_Recv(&(xsc2[0][0]),Nd*Nd,MPI_INT,rank2,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                double logtgtc2c2=logtargetIsing(xsc2,Nd,temps[glbc2]);
                double logtgtc2c1=logtargetIsing(xsc2,Nd,temps[glbc1]);


                originallog=logtgtc1c1+logtgtc2c2;
                proplog=logtgtc1c2+logtgtc2c1;

                accpt=exp(proplog-originallog);


                coin=unifrnd(0,1);
                if (coin<accpt){
                    // We indeed accept the proposal, notify the other chain
                    int accptstatus=1;

                    MPI_Send(&accptstatus,1,MPI_INT,rank2,0,MPI_COMM_WORLD);

                    MPI_Send(xs[c1],Nd*Nd,MPI_INT,rank2,0,MPI_COMM_WORLD);

                    deepcopy2Darray(xsc2,xs[c1], Nd, Nd);

                    exchangetimes+=1;

                } else {

                    int accptstatus=0;
                    MPI_Send(&accptstatus,1,MPI_INT,rank2,0,MPI_COMM_WORLD);
                }

            } else {
                int c2=glbc2%S;
                int accptstatus;

                MPI_Send(xs[c2],Nd*Nd,MPI_INT,rank1,0,MPI_COMM_WORLD);
                MPI_Recv(&accptstatus,1,MPI_INT,rank1,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                if (accptstatus==1)
                {
                    int **xsc1;
                    create2Dmemory(xsc1,Nd,Nd);
                    MPI_Recv(xsc1,Nd*Nd,MPI_INT,rank1,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    deepcopy2Darray(xsc1,xs[c2], Nd,Nd);
                    exchangetimes+=1;
                }
            }
        }


    }

    // At this point, broadcast result from rank 0  process to all MPI processes
    //MPI_Bcast(result, Nd, MPI_INT, 0, MPI_COMM_WORLD);

    print2Darray(partialresult,iterNum,S);
    free3Dmemory(xs, S, Nd,Nd);
}

void oneChainIsing(int **x, int T, int Nd, double temp)
{
    #pragma omp parallel shared(x)
    {
        for (int t=0; t<T; ++t)
        {
            int threadid = omp_get_thread_num();
            int numthreads = omp_get_num_threads();

            if (Nd/numthreads<=1)
            {
                throw "Too many threads! One thread must have at least 2 rows";
            }

            // Partition the matrix in strips
            int low = Nd*threadid/numthreads;
            int high=Nd*(threadid+1)/numthreads;

            if (high>Nd)
            {
                high=Nd;
            }

            // The schedule now is to update everything except the last row
            for(int i=low; i<high-1; ++i)
            {
                for(int j=0; j<Nd; ++j)
                {
                    int leftj=j-1;
                    int rightj=j+1;
                    if (j==0)
                    {
                        leftj=Nd-1;
                    }
                    else if (j==(Nd-1))
                    {
                        rightj=0;
                    }

                    int upi=i-1;
                    int downi=i+1;
                    if (i==0)
                    {
                        upi=Nd-1;
                    }else if (i==Nd-1)
                    {
                        downi=0;
                    }

                    double s=x[upi][j]+x[downi][j]+x[i][leftj]+x[i][rightj];
                    double cond_p=exp(temp*s)/(exp(temp*s)+exp(-temp*s));
                    if (unifrnd(0,1)<cond_p)
                    {
                        x[i][j]=1;
                    } else{
                        x[i][j]=0;
                    }
                }
            }

            // put a barrier here go ensure all threads update the last row on new values from neighboring regions
            #pragma omp barrier

            int i=high-1;

            for (int j=0;j<Nd;++j)
            {
                int leftj=j-1;
                int rightj=j+1;
                if (j==0)
                {
                    leftj=Nd-1;
                }
                else if (j==(Nd-1))
                {
                    rightj=0;
                }

                int upi=i-1;
                int downi=i+1;
                if (i==0)
                {
                    upi=Nd-1;
                }else if (i==Nd-1)
                {
                    downi=0;
                }

                double s=x[upi][j]+x[downi][j]+x[i][leftj]+x[i][rightj];
                double cond_p=exp(temp*s)/(exp(temp*s)+exp(-temp*s));
                if (unifrnd(0,1)<cond_p)
                {
                    x[i][j]=1;
                } else{
                    x[i][j]=0;
                }
            }

            // Put a barrier here to ensure iteration is synchronized each step
            #pragma omp barrier

        }
    }

}

double t(int **x, int Nd)
{
    double result=0;

#pragma omp parallel for reduction(+:result)
    for(int i=0; i<Nd; ++i)
    {
        for(int j=0; j<Nd; ++j)
        {
            int leftj=j-1;
            int rightj=j+1;
            if (j==0)
            {
                leftj=Nd-1;
            }
            else if (j==(Nd-1))
            {
                rightj=0;
            }

            int upi=i-1;
            int downi=i+1;
            if (i==0)
            {
                upi=Nd-1;
            }else if (i==Nd-1)
            {
                downi=0;
            }

            result+=x[i][j]*(x[upi][j]+x[downi][j]+x[i][leftj]+x[i][rightj]);
        }

    }
    return result;
}

double logtargetIsing(int **x, int Nd,double temp)
{
    double result=0;

    result=t(x,Nd);

    result=exp(temp*0.5*result);
    return result;
}