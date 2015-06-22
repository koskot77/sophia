#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <list>

#include <thread> // std::thread
#include <mutex>  // std::mutex

#include "./toolbox.h"

// Comile with:
// g++ -Wl,--no-as-needed -g -Wall -std=c++0x -o barcodes barcodes.cc -lpthread

using namespace std;

const char *fileName = "../Sophia/testReads.fastq";

#define NREC (10000)
size_t nReads = 0;

vector<string> identifier(NREC);
vector<string> sequence  (NREC);
vector<string> quality   (NREC);

////////////////////// union-find data structure //////////////////////
class UnionFind {
private:
    map<int,int>         node2cluster;
    map<int, list<int> > cluster2nodes;

public:

    int joinClusters(int cluster1, int cluster2){
        if( cluster1 == cluster2 ) return cluster1;
        list<int> &nodes1 = cluster2nodes[cluster1];
        list<int> &nodes2 = cluster2nodes[cluster2];
        size_t size1 = nodes1.size();
        size_t size2 = nodes2.size();
        if( size1==0 || size2==0 ) return 0;
        int newCluster = 0;
        if( size1 < size2 ){
            newCluster = cluster2;
            for(list<int>::const_iterator n = nodes1.begin(); n != nodes1.end(); n++)
                node2cluster[*n] = newCluster;
            nodes2.insert(nodes2.end(),nodes1.begin(),nodes1.end());
            cluster2nodes.erase(cluster1);
        } else {
            newCluster = cluster1;
            for(list<int>::const_iterator n = nodes2.begin(); n != nodes2.end(); n++)
                node2cluster[*n] = newCluster;
            nodes1.insert(nodes1.end(),nodes2.begin(),nodes2.end());
            cluster2nodes.erase(cluster2);
        }
        return 0;
    }

    int findCluster(int node) { return node2cluster[node]; }

    int nClusters(void) const { return cluster2nodes.size(); }

    const map<int, list<int> >& clusters(void) const { return cluster2nodes; }

    UnionFind(int maxNodes){
        for(int i=1; i<=maxNodes; i++){
             node2cluster [i] = i;
             cluster2nodes[i].push_back(i);
        }
    }
};
///////////////////////////////////////////////////////////////////////


/////////////// grouping reads with similar beginnings ////////////////

// check first 'barcodeWidth' positions and request identity of 'viewWidth' consecutive symbols
const size_t barcodeWidth = 15;
const size_t viewWidth    = 11;

// mutex for the union-find data stracture preventing data race
std::mutex mtx;
UnionFind *uf;

// check all reads after the 'read' and combine those with similar beginnings
void groupMatchesFor(size_t read){
    // safety
    if( read >= nReads ) return;

    // reference sequence
    const char *refSeq = sequence[read].c_str();

    // ignore short sequences
    if( strlen(refSeq) < barcodeWidth ) return;

    // break the read into patterns of viewWidth size
    char matchPattern[MAX_ADAPTORS][MAX_LENGTH];
    bzero(matchPattern, sizeof(matchPattern));

    for(size_t i=0; i<barcodeWidth-viewWidth; i++){
         if( i>=MAX_ADAPTORS ){ cout<<"Long record ... exiting"<<endl; exit(0); }
         strncpy(matchPattern[i],refSeq+i,viewWidth);
         matchPattern[i][viewWidth] = '\0';
    }

    // build the hash from the patterns
    LookUpTable lookUp(matchPattern);

    // check quality of the hash table
    size_t maxCollisions = 0;
    lookUp.countCollisions( maxCollisions );
    if( maxCollisions >= MAX_COLLISIONS ){
        for(size_t k=0; k<MAX_ADAPTORS; k++){
            size_t length = strlen( matchPattern[k] );
            if( length > 0 ){
                unsigned short errorPos = 0;
                unsigned long long coreCode = sequence2number(matchPattern[k],strlen(matchPattern[k]),errorPos);
                printf("%ld: %s, bucket: %lld\n",k,matchPattern[k],coreCode % BUCKETS);
            }
        }
        exit(0); // do not tolerate broken LUT
    }

    // look for other reads matching any of the patterns of the reference record
    for(size_t read2=read+1; read2<nReads; read2++){

        // probe sequence
        const char *probSeq = sequence[read2].c_str();

        // ignore short sequences
        if( strlen(probSeq) < barcodeWidth ) continue;

        // expect barcodes showing only in the beginning and ignore the rest of the sequence 
        char beginning[barcodeWidth+1];
        strncpy(beginning,probSeq,barcodeWidth);
        beginning[barcodeWidth] = '\0';

        NumericSequence numSeq(beginning);

        for(size_t i=0; i<barcodeWidth-viewWidth; i++){
            unsigned long long view = numSeq.view(i,viewWidth);
            size_t matchLength = viewWidth;
            if( lookUp.find(view,matchLength) != MAX_ADAPTORS ){ // found a match
                mtx.lock(); 
                int cluster1 = uf->findCluster(read+1);
                int cluster2 = uf->findCluster(read2+1);
                if( cluster1 != cluster2 ) uf->joinClusters(cluster1,cluster2);
                mtx.unlock(); 
            }
        }
    }
    cout<<"Read "<<read<<" done"<<endl;

}

void processReads(size_t from, size_t to){

     if( from>=to ) return;

     for(size_t read=from; read<to; read++)
         groupMatchesFor(read);

}

///////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]){
    // open input file:
    size_t iteration = 0; //atol(argv[1]);
    const size_t nCycles = 28069;

    ifstream input(fileName);

    if(!input){ cout<<"Cannot open "<<fileName<<endl; return 0; }

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

    uf = new UnionFind(nReads);

    const size_t numThreads = 3;
    std::thread threads[numThreads];

    for(size_t t=0; t<numThreads; t++)
        threads[t] = std::thread(processReads, nCycles*(t + iteration), nCycles*(t + iteration + 1));

    for(size_t t=0; t<numThreads; t++)
        threads[t].join();

    cout<<"Found "<<uf->nClusters()<<" clusters"<<endl;
    const map<int, list<int> > &clusters = uf->clusters();

    ofstream output("output.csv");
    if( !output ){ cout<<"Cannot open "<<"output.csv"<<endl; return 0; }
    for(map<int, list<int> >::const_iterator iter = clusters.begin(); iter != clusters.end(); iter++){
        int seed = iter->first;
        size_t nClusters = iter->second.size();
        if( nClusters == 0 ) cerr<<" Error: empty cluster for "<<seed<<"!"<<endl;
        output<<seed-1<<","<<nClusters<<endl;
//        if(seed==1)
//        for(list<int>::const_iterator node = iter->second.begin(); node != iter->second.end(); node++)
//            output<<seed-1<<","<<*node-1<<endl;
    }
    output.close();

    for(map<int, list<int> >::const_iterator iter = clusters.begin(); iter != clusters.end(); iter++){
        int seed = iter->first;
        if( iter->second.size() == 0 ) cerr<<" Error: empty cluster for "<<seed<<"!"<<endl;

        if( iter->second.size() < 1000 ) continue;

        stringstream fname;
        fname<<"output"<<(seed-1)<<".fastq";
        ofstream output(fname.str());
        if( !output ){ cout<<"Cannot open "<<fname.str()<<endl; return 0; }

        for(list<int>::const_iterator node = iter->second.begin(); node != iter->second.end(); node++){
            output<<identifier[*node-1]<<endl;
            output<<sequence  [*node-1]<<endl;
            output<<"+"<<endl;
            output<<quality   [*node-1]<<endl;
        }

        output.close();
    }

    return 0;
}
