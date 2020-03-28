/* MCMC.cpp */
double logtarget(int *x, int Nd, int **R, int **C, double temp);
void oneChain(int *x0, int T, int Nd, int *xT, int **R, int **C, double temp);
void temperedChains(int iterNum, int totalS, int Nd, int T, int **R, int **C, double *temps, int * result, int rank, int size);
void rotateout(int **xs, int S, int Nd, int U, int **R, int **C, double temp);

/* Randomize.cpp */
double unifrnd(double a, double b);
double normrnd(double mu, double std);
int unifrndint(int a, int b);
int rndpermutation(int *inputarr, int inputarr_len, int *outputarr);
int rndswap(int *x, int x_len, int *x_swapped);
void rndkofn(int *input, int *output, int n, int k);

/* ArrayUtilities.cpp */
int assignRow(int **A, int a[], int N, int r);
void deepcopy1Darray(int *input, int *output, int Nd);
void print1Darray(int* x, int Nd);
void print2Darray(int **myArray, int height, int width);
void init2Darray(int **A, int height, int width, int filler);
double GetAverage(double num[], int n);
double GetStd(double num[], int n);
void free2Dmemory(int **ptr, int height, int width);
void create2Dmemory(int **&ptr, int height, int width);

/* CypherUtilities.cpp */
void buildTransitionMat(int **R, int numchar, std::string file);
int char2num(char ch);
char num2char(int num);
char cipher(char ch, int *cipherkey);
char decipher(char ch, int *decipherkey);
void cipherkey2decipherkey(int *cipherkey, int *decipherkey, int Nd);
void buildCiphered(std::string inputfile, std::string outputfile, int *cipherkey);
void buildDeciphered(std::string inputfile, std::string outputfile, int *decipherkey);









