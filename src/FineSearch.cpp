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




bool isaphbt(char ch)
{
    return ((ch>='a' && ch<='z') || (ch>='A' && ch<='Z'));
}



/*
 * Read in the word frequency count text and convert it into a c++ map for fast access
 *
 * Function Arguments:
 * dicttext: the path of the word frequency count text
 * dict: the dictionary with word as key and frequency as value
 *
 * */
void buildWordsFreqMap(std::string dicttext, std::map<std::string, int> &dict)
{
    std::ifstream file(dicttext);
    std::string s;
    int ranking=0;

    while (std::getline(file, s)) {
        dict[s]=ranking;
        ranking=ranking+1;
    }
}


void CWFScore(std::string inputstr, int &CWFS, double &percwords, std::map<std::string, int> g_dict, int &wordcount)
{
    std::string word="";
    double numcorrect=0.0;
    double numofwords=0;
    CWFS=0;
    percwords=0;

    for (int i=0; i<(inputstr.length());++i)
    //for (int i=100; i<250;++i)
    {
        if (isaphbt(inputstr[i]))
        {
            word+=tolower(inputstr[i]);
        }
        else if (word !="")
        {
            //printf("%s\n",word.c_str());
            if (g_dict.count(word)>0)
            {
                CWFS+=g_dict.at(word);
                numcorrect+=1;
            } else{
                CWFS+=10000000;
            }
            word="";
            numofwords+=1;
        }
    }
    percwords=numcorrect/numofwords;
    wordcount=numofwords;
}



void fineOptimize(int * x, int Nd, int TT, std::string cipheredstring, int * output, std::map<std::string, int> g_dict, int rank)
{
    std::string decipheredstringprev=buildDecipheredstring(cipheredstring, x);

    int CWFSprev;
    double percwordsprev;
    int wordcountprev;
    CWFScore(decipheredstringprev, CWFSprev, percwordsprev, g_dict, wordcountprev);

    int CWFSprop;
    double percwordsprop;
    int wordcountprop;

    int xprop[Nd];
    int idx1, idx2;


    for (int t=0; t<TT; ++t)
    {


        rndswapwtidx(x,Nd,xprop,idx1,idx2);

        if (!((isaphbt(x[idx1]))&&(isaphbt(x[idx2]))))
        {
            continue;
        }

        std::string decipheredstringprop=buildDecipheredstring(cipheredstring, xprop);
        CWFScore(decipheredstringprop, CWFSprop, percwordsprop, g_dict,wordcountprop);
        printf("Prev: %f\n",percwordsprev);
        printf("Proposal: %f\n",percwordsprop);

        if ((percwordsprev<percwordsprop)&&(wordcountprop>=wordcountprev))
        {
            deepcopy1Darray(xprop,x,Nd);
            percwordsprev=percwordsprop;
            CWFSprev=CWFSprop;
            wordcountprev=wordcountprop;
            printf("Hello this is process %d\n", rank);
            printf("Number of char: %d\n",wordcountprop);
            printf("%s\n",decipheredstringprop.c_str());

            // At this point, broadcasting x to all other MPI processes' x.
            MPI_Bcast(x,Nd,MPI_INT,rank,MPI_COMM_WORLD);
        }
    }

    deepcopy1Darray(x,output,Nd);

}

