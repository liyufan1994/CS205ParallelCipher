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

void deepcopy1Darray(int *input, int *output, int Nd)
{
    for (int i=0; i<Nd; ++i)
    {
        output[i]=input[i];
    }
}


/*
 * This function fill each entry of A with filler
 * */
void init2Darray(int **A, int height, int width, int filler)
{
    for (int i=0; i<height; ++i)
    {
        for (int j=0; j<width; ++j)
        {
            A[i][j]=filler;
        }
    }
}

void print1Darray(int* x, int Nd)
{
    for (int i=0; i<Nd; ++i)
    {
        std::cout << x[i] << ' ';
    }
    std::cout << std::endl;
}



void create2Dmemory(int **&ptr, int height, int width)
{
    ptr=new int *[height];
    for(int i = 0; i <height; ++i)
        ptr[i] = new int[width];
}

void free2Dmemory(int **ptr, int height, int width)
{
    // Free each subarray
    for (int i=0; i<height; ++i)
    {
        delete[] ptr[i];
    }

    //Free the array of pointers
    delete[] ptr;
}




void print2Darray(int **myArray, int height, int width)
{
    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            std::cout << myArray[i][j] << ' ';
        }
        std::cout << std::endl;
    }
}






double GetAverage(double num[], int n)
{
    double sum = 0.0, avg;
    for(int i = 0; i < n; i++)
        sum += num[i];
    avg = sum / n;

    return avg;
}

double GetStd(double num[], int n)
{
    double avg=GetAverage(num,n);
    double sum;
    for(int i = 0; i < n; ++i)
        sum += (num[i]-avg)*(num[i]-avg);
    double stdev = std::sqrt(sum / (n-1));

    return stdev;
}


/*
 * This is function deep copies content of input array a (of length N) to the r-th row of 2D array A.
 *
 * */

int assignRow(int **A, int a[], int N, int r)
{

    for (int i=0; i<N; i++)
    {
        A[r][i]=a[i];
    }
    return 0;
}
