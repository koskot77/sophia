#include <stdio.h>  // sscanf
#include <stdlib.h> // atoi
#include <math.h>   // sqrt
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <map>

#include "toolbox.h"

#define NREC (10000)

using namespace std;
vector<string> identifier(NREC);
vector<string> sequence  (NREC);
vector<string> quality   (NREC);
map<size_t,size_t> errors; // record number and position of the unrecognized symbol
map<unsigned long long,unsigned int>  counts;         // pattern occupancy

int main(int argc, char *argv[]){
    if( argc < 2 ) return 0;

    const char *fileName = argv[1]; 
    size_t barcodeWidth = 12;
    size_t viewWidth    = 10;

    if( argc>3 ){
        barcodeWidth = strtol(argv[2],NULL,0);
        viewWidth    = strtol(argv[3],NULL,0);
    }

    // open input file
    ifstream input(fileName);
    if(!input){ cout<<"Cannot open "<<fileName<<endl; return 0; }

    // read the records
    size_t nReads=0;
    while( !input.eof() ){

        if( nReads == identifier.size() ){
            identifier.resize( identifier.size() + NREC );
            sequence.  resize( sequence.  size() + NREC );
            quality.   resize( quality.   size() + NREC );
        }

        string tmp;
        getline(input, identifier[nReads], '\n' ); 
        getline(input, sequence  [nReads], '\n' );
        getline(input, tmp,                '\n' );
        getline(input, quality   [nReads], '\n' );

        if( identifier[nReads][0] != '@' || tmp[0] != '+') break;
        nReads++;
    }
    input.close();

    identifier.resize(nReads);
    sequence.  resize(nReads);
    quality.   resize(nReads);

    cout<<"Reads: "<<nReads<<endl;

    // run a density check
    for(size_t read=0; read<nReads; read++){
        const char *seq = sequence[read].c_str();
        NumericSequence numSeq(seq);

        if( numSeq.error() )
            errors[read] = numSeq.error();

//        for(size_t i=0; i<strlen(seq)-viewWidth; i++){
        for(size_t i=0; i<barcodeWidth-viewWidth; i++){
            // if there was an unknown symbol in the sequence, discard every view that includes it 
            if( numSeq.error() - i < viewWidth ) continue;

            unsigned long long view = numSeq.view(i,viewWidth);
            counts[view]++;
        }
    }

    cout<<"Found unique sequences of length "<<viewWidth<<": "<<counts.size()<<endl;

    ofstream output("output.csv");

    if( !output ){ cout<<"Cannot open output.csv"<<endl; return 0; }

    output<<"ID,label,count"<<endl;

    for(map<unsigned long long,unsigned int>::const_iterator seq = counts.begin(); seq != counts.end(); seq++)
        output<<seq->first<<","<<number2sequence(seq->first,viewWidth)<<","<<seq->second<<endl;
    output.close();

    return 0;
}
