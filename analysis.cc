#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include "splitter.cc"

const char *fileName = "../Sophia/testReads.fastq";

#define NREC (10000)

using namespace std;
vector<string> identifier(NREC);
vector<string> sequence  (NREC);
vector<string> quality   (NREC);
map<unsigned long long,unsigned int> counts;

int main(int argc, char *argv[]){
    // open input file:
    ifstream input(fileName);

    if(!input){ cout<<"Cannot open "<<inputFileName<<endl; return 0; }

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

    const int viewWidth = atoi(argv[1]);

    // run a density check
    for(size_t read=0; read<nReads; read++){
        const char *seq = sequence[read].c_str();
        NumericSequence numSeq(seq);
        for(size_t i=0; i<strlen(seq)-viewWidth; i++){
            unsigned long long view = numSeq.view(i,viewWidth);
            counts[view]++;
        }
    }

    cout<<"Found unique sequences of length "<<viewWidth<<": "<<counts.size()<<endl;

    ofstream output("output.csv");

    if( !output ){ cout<<"Cannot open "<<"output.csv"<<endl; return 0; }

    output<<"ID,label,count"<<endl;
    for(map<unsigned long long,unsigned int>::const_iterator seq = counts.begin(); seq != counts.end(); seq++)
        output<<seq->first<<","<<number2sequence(seq->first,viewWidth)<<","<<seq->second<<endl;

    output.close();

    return 0;
}
