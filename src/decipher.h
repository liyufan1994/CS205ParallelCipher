/* MCMC.cpp */
double logtarget(int *x, int Nd, int **R, int **C, double temp);
void oneChain(int *x0, int T, int Nd, int *xT, int **R, int **C, double temp);
void temperedChains(int iterNum, int totalS, int Nd, int T, int **R, int **C, double *temps, int * result, int rank, int size);
void rotateout(int **xs, int S, int Nd, int U, int **R, int **C, double temp);



/* CypherUtilities.cpp */
void buildTransitionMat(int **R, int numchar, std::string file);
int char2num(char ch);
char num2char(int num);
char cipher(char ch, int *cipherkey);
char decipher(char ch, int *decipherkey);
void cipherkey2decipherkey(int *cipherkey, int *decipherkey, int Nd);
void buildCiphered(std::string inputfile, std::string outputfile, int *cipherkey);
void buildDeciphered(std::string inputfile, std::string outputfile, int *decipherkey);
std::string buildDecipheredstring(std::string inputstring,  int *decipherkey);
std::string readcodefile(std::string inputfile);

/* FineSearch.cpp */
void buildWordsFreqMap(std::string dicttext, std::map<std::string, int> &dict);
void CWFScore(std::string inputstr, int &CWFS, double &percwords, std::map<std::string, int> g_dict, int &wordcount);
void fineOptimize(int * x, int Nd, int TT, std::string cipheredstring, int * output, std::map<std::string, int> g_dict, int rank);


