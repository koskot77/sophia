#include <stdio.h>  // sscanf
#include <stdlib.h> // atoi
#include <math.h>   // sqrt
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <map>

using namespace std;
#include "DNASequencing.cc"

#define BLOCK_SIZE (1000)

int main(int argc, char *argv[]){

    const char *refFileNames[24] = {
        "data/chromatid1.fa",
        "data/chromatid2.fa",
        "data/chromatid3.fa",
        "data/chromatid4.fa",
        "data/chromatid5.fa",
        "data/chromatid6.fa",
        "data/chromatid7.fa",
        "data/chromatid8.fa",
        "data/chromatid9.fa",
        "data/chromatid10.fa",
        "data/chromatid11.fa",
        "data/chromatid12.fa",
        "data/chromatid13.fa",
        "data/chromatid14.fa",
        "data/chromatid15.fa",
        "data/chromatid16.fa",
        "data/chromatid17.fa",
        "data/chromatid18.fa",
        "data/chromatid19.fa",
        "data/chromatid20.fa",
        "data/chromatid21.fa",
        "data/chromatid22.fa",
        "data/chromatid23.fa",
        "data/chromatid24.fa"
    }; 

    DNASequencing worker;
    worker.initTest(0); // here we may optimize for k-mers sizes, hash table parameters, etc.

{ // save some space by getting rid of the local container on leaving the scope once we hand over the results to DNASequencing worker
    vector<string> chromatidSequence[24];

    for(unsigned int chromatidId=19; chromatidId<20; chromatidId++){
        // open input file
        ifstream input( refFileNames[chromatidId] );
        if( !input ){ cout<<"Cannot open "<<refFileNames[chromatidId]<<endl; return 0; }

        // read comment line
        string tmp;
        getline(input, tmp, '\n'); 
        cout<<"Reading "<<refFileNames[chromatidId]<<" starting with: "<<endl<< tmp <<endl;

        size_t nLines=0;
        while( !input.eof() ){

            if( nLines == chromatidSequence[chromatidId].size() )
                chromatidSequence[chromatidId].resize( chromatidSequence[chromatidId].size() + BLOCK_SIZE );

            getline(input, chromatidSequence[chromatidId][nLines], '\n' );

            nLines++;
        }
        chromatidSequence[chromatidId].resize(nLines);

        cout<<"complete "<<nLines<<endl;
        input.close();
    }
/*
    worker.passReferenceGenome(0, chromatidSequence[0]);
    worker.passReferenceGenome(1, chromatidSequence[1]);
    worker.passReferenceGenome(2, chromatidSequence[2]);
    worker.passReferenceGenome(3, chromatidSequence[3]);
    worker.passReferenceGenome(4, chromatidSequence[4]);
    worker.passReferenceGenome(5, chromatidSequence[5]);
    worker.passReferenceGenome(6, chromatidSequence[6]);
    worker.passReferenceGenome(7, chromatidSequence[7]);
    worker.passReferenceGenome(8, chromatidSequence[8]);
    worker.passReferenceGenome(9, chromatidSequence[9]);
    worker.passReferenceGenome(10, chromatidSequence[10]);
    worker.passReferenceGenome(11, chromatidSequence[11]);
    worker.passReferenceGenome(12, chromatidSequence[12]);
    worker.passReferenceGenome(13, chromatidSequence[13]);
    worker.passReferenceGenome(14, chromatidSequence[14]);
    worker.passReferenceGenome(15, chromatidSequence[15]);
    worker.passReferenceGenome(16, chromatidSequence[16]);
    worker.passReferenceGenome(17, chromatidSequence[17]);
    worker.passReferenceGenome(18, chromatidSequence[18]);
*/
    worker.passReferenceGenome(19, chromatidSequence[19]);
/*
    worker.passReferenceGenome(20, chromatidSequence[20]);
    worker.passReferenceGenome(21, chromatidSequence[21]);
    worker.passReferenceGenome(22, chromatidSequence[22]);
    worker.passReferenceGenome(23, chromatidSequence[23]);
*/
}

    worker.preProcessing();


    const char *readFileNames1[1] = {
        "data/small10.fa1"
    }; 
    const char *readFileNames2[1] = {
        "data/small10.fa2"
    }; 
    vector<string> readName[1];
    vector<string> readSequence[1];

    for(unsigned int readFileId=0; readFileId<1; readFileId++){
        // open input file
        ifstream input1( readFileNames1[readFileId] );
        if( !input1 ){ cout<<"Cannot open "<<readFileNames1[readFileId]<<endl; return 0; }

        ifstream input2( readFileNames2[readFileId] );
        if( !input2 ){ cout<<"Cannot open "<<readFileNames2[readFileId]<<endl; return 0; }

        size_t nLines=0;
        while( !input1.eof() && !input2.eof() ){

            if( nLines == readName[readFileId].size() ){
                readName    [readFileId].resize( readName    [readFileId].size() + BLOCK_SIZE );
                readSequence[readFileId].resize( readSequence[readFileId].size() + BLOCK_SIZE );
            }

            getline(input1, readName    [readFileId][nLines], '\n' );
            getline(input1, readSequence[readFileId][nLines], '\n' );

//            readSequence[readFileId][nLines].erase( readSequence[readFileId][nLines].length()-1 ); // Fucking Windows eol definition!

            nLines++;

            getline(input2, readName    [readFileId][nLines], '\n' );
            getline(input2, readSequence[readFileId][nLines], '\n' );

//            readSequence[readFileId][nLines].erase( readSequence[readFileId][nLines].length()-1 ); // Fucking Windows eol definition!

            nLines++;
        }
        readName    [readFileId].resize(nLines);
        readSequence[readFileId].resize(nLines);

        cout<<"complete "<<nLines<<endl;
        input1.close();
        input2.close();

        vector<string> qwe = worker.getAlignment(18*2, 0, 0, readName[readFileId], readSequence[readFileId]);
    }


    return 0;
}
