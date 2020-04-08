/* ArrayUtilities.cpp */
int assignRow(int **A, int a[], int N, int r);
void deepcopy1Darray(int *input, int *output, int Nd);
void print1Darray(int* x, int Nd);
void print2Darray(int **myArray, int height, int width, std::string filenm="std");
void init2Darray(int **A, int height, int width, int filler);
double GetAverage(double num[], int n);
double GetStd(double num[], int n);
void free2Dmemory(int **ptr, int height, int width);
void create2Dmemory(int **&ptr, int height, int width);
void create3Dmemory(int ***&ptr, int layers, int height, int width);
void free3Dmemory(int ***&ptr, int layers, int height, int width);
void deepcopy2Darray(int **input, int **output, int height, int width);
void conv2Dto1D(int ** input, int *output, int height, int width);
void conv1Dto2D(int *input, int **output, int height, int width);
