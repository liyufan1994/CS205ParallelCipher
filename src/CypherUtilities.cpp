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



void buildDeciphered(std::string inputfile, std::string outputfile, int *decipherkey)
{
    std::fstream fin(inputfile, std::fstream::in);
    std::fstream fout(outputfile, std::fstream::out);
    char ch;

    while (fin >> std::noskipws >> ch) {
        fout<<decipher(ch, decipherkey);
    }
}


/*
 * Given file path, read it and return it as a string
 * */
std::string readcodefile(std::string inputfile)
{
    std::fstream fin(inputfile, std::fstream::in);
    char ch;
    std::string string2ret;

    while (fin >> std::noskipws >> ch) {
        string2ret+=ch;
    }
    return string2ret;
}


/*
 * Receive an input string and decipher it using decipher key and return deciphered string
 *
 * */
std::string buildDecipheredstring(std::string inputstring,  int *decipherkey)
{
    std::string tempstr="";

    for (char const &c: inputstring) {
        tempstr+=decipher(c, decipherkey);
    }

    return tempstr;
}


void cipherkey2decipherkey(int *cipherkey, int *decipherkey, int Nd)
{
    int container[Nd];
    for (int i=0; i<Nd; ++i)
    {
        container[cipherkey[i]]=i;
    }
    deepcopy1Darray(container,decipherkey,Nd);
}


/*
 * This function receives a character and a cipher key and returns the ciphered character by the key
 *
 *
 * */
char cipher(char ch, int *cipherkey)
{
    return num2char(cipherkey[char2num(ch)]);
}

char decipher(char ch, int *decipherkey)
{
    return num2char(decipherkey[char2num(ch)]);
}


/*
 *
 * This function receives an input file of the original text and write the ciphered text to outputfile
 * according to the cipher key it receives.
 *
 * */
void buildCiphered(std::string inputfile, std::string outputfile, int *cipherkey)
{
    std::fstream fin(inputfile, std::fstream::in);
    std::fstream fout(outputfile, std::fstream::out);
    char ch;

    while (fin >> std::noskipws >> ch) {
        fout<<cipher(ch, cipherkey);
    }


}



/*
 * This function returns the numeric index of the character it receives. The indexing here is based on printable
 * ascii code from 32-126 inclusive--a total of 95 characters.
 * */
int char2num(char ch)
{
    if ((int(ch)>=32) && (int(ch)<=127))
    {
        return int(ch-32);
    }else{
        return 0;
    }
}

char num2char(int num)
{

    return char(num+32);

}

void buildTransitionMat(int **R, int numchar, std::string file)
{

    char ch1;
    char ch2;

    for (int i=0; i<numchar; ++i)
    {
        for (int j=0; j<numchar; ++j)
        {
            R[i][j]=1;
        }
    }

    std::fstream fin(file, std::fstream::in);
    bool isFirst=1;
    while (fin >> std::noskipws >> ch2) {
        if (isFirst) {
            ch1=ch2;
            isFirst = 0;
        } else {
            R[char2num(ch1)][char2num(ch2)]+=1;
            ch1=ch2;
        }
    }


}





