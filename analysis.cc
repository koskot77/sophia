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

const char *fileName = "../testReads.fastq";

#define NREC (10000)

using namespace std;
vector<string> identifier(NREC);
vector<string> sequence  (NREC);
vector<string> quality   (NREC);
map<size_t,size_t> errors; // record number and position of the unrecognized symbol
map<unsigned long long,unsigned int>  counts;         // pattern occupancy
map<unsigned long long,list<size_t> > pattern2record; // list of record numbers where the pattern is found
map<unsigned long long,list<double> > averageQuality; // 
map<unsigned long long,list<double> > averageAcuracy; // 

// find an average value on the array
double mean(const char *array, size_t length){
    if( length==0 ) return -1;
    int retval = 0;
    for(size_t i=0; i<length; i++) retval += array[i];
    return double(retval) / length;
}
double mean(const list<double> &array){
    if( array.size()==0 ) return -1;
    double retval = 0;
    for(list<double>::const_iterator i = array.begin(); i != array.end(); i++) retval += *i;
    return retval / array.size();
}
// find a variance on the array
double variance(const char *array, size_t length){
    if( length==0 ) return -1;
    int retval = 0;
    for(size_t i=0; i<length; i++) retval += array[i]*array[i];
    return double(retval) / length;
}
double variance(const list<double> &array){
    if( array.size()==0 ) return -1;
    double retval = 0;
    for(list<double>::const_iterator i = array.begin(); i != array.end(); i++) retval += (*i)*(*i);
    return retval / array.size();
}
const char phredMinScore = 33;
// average Phred quality
double adjustedMean(const char *qual, size_t length){
    return mean(qual,length) - phredMinScore;
}
// standard deviation of the Phred quality
double adjustedSD(const char *qual, size_t length){
    if( length==0 ) return -1;
    double mu = mean(qual,length);
    return sqrt( variance(qual,length) - mu*mu );
}
// find probability of an error in the sequence
double accuracy(const char *qual, size_t length){
    if( length==0 ) return -1;
    double retval = 1;
    for(size_t i=0; i<length; i++){
        double invErrProb   = pow(10,(qual[i]-phredMinScore)/10.);
        double accuracyProb = (1.-1./invErrProb);
        retval *= accuracyProb; //-log(accuracyProb);
    }
    return retval;
}

int main(int argc, char *argv[]){
    if( argc < 2 ) return 0;

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

    // length of the patterns
    const unsigned int viewWidth = atoi(argv[1]);
    bool comprehensive = (argc>2 ? true : false);

    // run a density check
    for(size_t read=0; read<nReads; read++){
        const char *seq = sequence[read].c_str();
        NumericSequence numSeq(seq);

        if( numSeq.error() )
            errors[read] = numSeq.error();

        for(size_t i=0; i<strlen(seq)-viewWidth; i++){
            // if there was an unknown symbol in the sequence, discard every view that includes it 
            if( numSeq.error() - i < viewWidth ) continue;

            unsigned long long view = numSeq.view(i,viewWidth);
            counts[view]++;
            if( comprehensive ){
                pattern2record[view].push_back(i);
                averageQuality[view].push_back( adjustedMean(quality[read].c_str()+i,viewWidth) );
                averageAcuracy[view].push_back( accuracy    (quality[read].c_str()+i,viewWidth) );
            }
        }
    }

    cout<<"Found unique sequences of length "<<viewWidth<<": "<<counts.size()<<endl;

    ofstream output("output.csv");

    if( !output ){ cout<<"Cannot open "<<"output.csv"<<endl; return 0; }

    output<<"ID,label,count"<<(comprehensive?",meanQual,sdQual,acc":"")<<endl;

    for(map<unsigned long long,unsigned int>::const_iterator seq = counts.begin(); seq != counts.end(); seq++){
        output<<seq->first<<","<<number2sequence(seq->first,viewWidth)<<","<<seq->second;
        if( comprehensive ){
            list<double> &q = averageQuality[seq->first];
            double meanQual = mean(q);
            double sdQual   = sqrt(variance(q) - meanQual*meanQual);
            list<double> &a = averageAcuracy[seq->first];
            double meanAcc  = mean(a);
            output<<","<<meanQual<<","<<sdQual<<","<<meanAcc<<endl;
        } else {
            output<<endl;
        }
    }
    output.close();

    // summarize the errors
    cout<<"Found "<<errors.size()<<" errors"<<endl;

    ofstream err("errors.csv");
    err<<"id,rank,x,y,len,rec,symb,qual"<<endl;
    for(map<size_t,size_t>::const_iterator e=errors.begin(); e!=errors.end(); e++){
        char id[128], rank[128], x[128], y[128], len[128];
        sscanf( identifier[e->first].c_str(), "@%s rank=%s x=%s y=%s length=%s",id,rank,x,y,len);
        err<<id<<","<<rank<<","<<x<<","<<y<<","<<len<<","
           <<e->first<<","<<sequence[e->first][e->second-1]<<","
           <<int(quality [e->first][e->second-1]-phredMinScore)
           <<endl;
    }
    err.close();

    return 0;
}
