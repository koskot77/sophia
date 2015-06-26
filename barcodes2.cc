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

#include <thread>
#include <future>
#include <chrono>

#include "./toolbox.h"

using namespace std;
// Comile with:
// g++ -Wl,--no-as-needed -g -Wall -std=c++11 -o barcodes2 barcodes2.cc -lpthread

////////////////////// routine for reading fastq //////////////////////

// keep reads from the input file in global environment
size_t nReads = 0;
#define NREC (10000)
vector<string> identifier(NREC);
vector<string> sequence  (NREC);
vector<string> quality   (NREC);

int readFile(const char *fileName){

    ifstream input(fileName);

    if(!input){ cout<<"Cannot open "<<fileName<<endl; return -1; }

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

    return 0;
}
#undef NREC
///////////////////////////////////////////////////////////////////////




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

    UnionFind(const list<int> &nodes){
        for(list<int>::const_iterator i=nodes.begin(); i!=nodes.end(); i++){
             node2cluster [*i] = *i;
             cluster2nodes[*i].push_back(*i);
        }
    }

    UnionFind(int maxNodes){
        for(int i=1; i<=maxNodes; i++){
             node2cluster [i] = i;
             cluster2nodes[i].push_back(i);
        }
    }

    UnionFind(int begin, int end){
        if( begin >= end || begin <= 0 ) exit(0);
        for(int i=begin; i<end; i++){
             node2cluster [i] = i;
             cluster2nodes[i].push_back(i);
        }
    }
};

// subsets of clusters to be combined in superclusters (one per thread)
UnionFind* uf[100];
///////////////////////////////////////////////////////////////////////




/////////////// grouping reads with similar beginnings ////////////////

bool isLUTgood(const LookUpTable &lut, const char matchPattern[MAX_ADAPTORS][MAX_LENGTH]){
    bool isGood = true;

    // check quality of the hash table
    size_t maxCollisions = 0;
    lut.countCollisions( maxCollisions );
    if( maxCollisions >= MAX_COLLISIONS ){
        for(size_t k=0; k<MAX_ADAPTORS; k++){
            size_t length = strlen( matchPattern[k] );
            if( length > 0 ){
                unsigned short errorPos = 0;
                unsigned long long coreCode = sequence2number(matchPattern[k],strlen(matchPattern[k]),errorPos);
                printf("%ld: %s, bucket: %lld\n",k,matchPattern[k],coreCode % BUCKETS);
            }
        }
        isGood = false;
    }

    return isGood;
}

// check first 'barcodeWidth' positions and request identity of 'viewWidth' consecutive symbols
const size_t barcodeWidth = 12;
const size_t viewWidth    = 10;
/*
size_t chop2patterns(const char *seq, char matchPattern[MAX_ADAPTORS][MAX_LENGTH]){
    // do nothing for short sequences
    if( strlen(seq) < barcodeWidth ) return 0;
    // clear
    bzero(matchPattern, sizeof(matchPattern));
    // slide along the beginning of the sequence 
    size_t  nPatterns = 0;
    for(  ; nPatterns < barcodeWidth-viewWidth; nPatterns++){
        if( nPatterns>=MAX_ADAPTORS ){
            cerr<<"Long record ... exiting"<<endl;
            return 0;
        }
        strncpy(matchPattern[nPatterns], seq + nPatterns, viewWidth);
        matchPattern[nPatterns][viewWidth] = '\0';
    }
    return nPatterns;
}
*/

// check all reads after the 'read' and combine those with similar beginnings
bool groupMatchesFor(size_t read, size_t block, size_t till=nReads){
    // do nothing if read doesn't exist
    if( read >= nReads ) return true;

    // reference sequence
    const char *refSeq = sequence[read].c_str();

    // do nothing for short sequences
    if( strlen(refSeq) < barcodeWidth ) return true;

    // break the read into patterns of viewWidth size
    char matchPattern[MAX_ADAPTORS][MAX_LENGTH];
    bzero(matchPattern, sizeof(matchPattern));

    for(size_t i=0; i<barcodeWidth-viewWidth; i++){
        if( i>=MAX_ADAPTORS ){ cout<<"Long record ... exiting"<<endl; return false; }
        strncpy(matchPattern[i],refSeq+i,viewWidth);
        matchPattern[i][viewWidth] = '\0';
    }

//    size_t nPatterns = chop2patterns(refSeq, matchPattern);
//    // do nothing for short sequences (not an error)
//    if( nPatterns==0 ) return true;

    // build the hash from the patterns
    LookUpTable lookUp(matchPattern);
    // check if it is good
    if( !isLUTgood(lookUp,matchPattern) ) return false; // do not tolerate broken LUT

    // look for other reads matching any of the patterns of the reference record
    for(size_t read2=read+1; read2<till && read2<nReads; read2++){

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
                int cluster1 = uf[block]->findCluster(read+1);
                int cluster2 = uf[block]->findCluster(read2+1);
                if( cluster1 != cluster2 ) uf[block]->joinClusters(cluster1,cluster2);
            }
        }
    }
    return true;
}

bool processReads(size_t begin, size_t end, size_t block){

     if( begin>=end ) return false;

     if( uf[block] ){
         cout<<"More than one call to process block #"<<block<<endl;
         return false;
     }

     uf[block] = new UnionFind(begin+1, end+1);

     for(size_t read=begin; read<end; read++)
         if( !groupMatchesFor(read,block,end) ) return false;

     return true;
}
///////////////////////////////////////////////////////////////////////


/////////// merge matching clusters from different groups /////////////

// build a search machinery
map<int,LookUpTable> buildSearchHelper(const UnionFind &uf){
    // clusters of the first group
    const map<int, list<int> > &clusters = uf.clusters();

    // assign every cluster with a search helper
    map<int,LookUpTable> retval;

    // initialize the search helpers
    for(map<int, list<int> >::const_iterator iter = clusters.begin(); iter != clusters.end(); iter++){
        int seed = iter->first;

        // do nothing for short sequences
        const char *seq = sequence[seed-1].c_str();
        if( strlen(seq) < barcodeWidth ) continue;

        size_t nClusters = iter->second.size();
        if( nClusters == 0 ) cerr<<" Error: empty cluster for "<<seed<<"!"<<endl;
        // construct a set of search patterns
        set<string> barcodeCandidates;
        // start with the parent read of the cluster
        const char *seedSeq = sequence[seed-1].c_str();
        // split the beginning into search views
        for(size_t i=0; i<barcodeWidth-viewWidth; i++){
            if( i>=MAX_ADAPTORS ){ cout<<"Long record ... exiting"<<endl; exit(0); }
            char tmp[viewWidth+1];
            strncpy(tmp,seedSeq+i,viewWidth);
            tmp[viewWidth] = '\0';
            barcodeCandidates.insert(tmp);
        }
        // do the same for every read from the cluster
        for(list<int>::const_iterator node = iter->second.begin(); node != iter->second.end(); node++){
            const char *seq = sequence[*node-1].c_str();
            // split the beginning into search views
            for(size_t i=0; i<barcodeWidth-viewWidth; i++){
                if( i>=MAX_ADAPTORS ){ cout<<"Long record ... exiting"<<endl; exit(0); }
                char tmp[viewWidth+1];
                strncpy(tmp,seq+i,viewWidth);
                tmp[viewWidth] = '\0';
                barcodeCandidates.insert(tmp);
            }
        }
        if( barcodeCandidates.size() >= MAX_ADAPTORS ){ cout<<"Too many search patterns: "<<barcodeCandidates.size()<<endl; exit(0); } 
        // convert the barcodeCandidates into the format accepted by the LookUpTable 
        char matchPattern[MAX_ADAPTORS][MAX_LENGTH];
        bzero(matchPattern, sizeof(matchPattern));
        size_t p=0;
        for(set<string>::const_iterator cand = barcodeCandidates.begin(); cand != barcodeCandidates.end(); cand++,p++){
            strncpy(matchPattern[p],cand->c_str(),viewWidth);
            matchPattern[p][viewWidth] = '\0';
        }
        pair< map<int,LookUpTable>::iterator, bool> res = retval.insert(pair<int,LookUpTable>(seed,LookUpTable(matchPattern)));
        if( !res.second ){ cout<<"Long record ... exiting"<<endl; exit(0); }

        // check quality of the hash table
        size_t maxCollisions = 0;
        res.first->second.countCollisions( maxCollisions );
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
    }

    return retval;
}
/*
UnionFind mergeGroups(UnionFind &group1, UnionFind &group2){

    map<int,LookUpTable> searchHelper = buildSearchHelper(group1);

}
*/
///////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]){

    if( readFile("./data.txt") ) return 0;
    else cout<<"Reads: "<<nReads<<endl;

    const size_t readsInBlock = 20000;
    const size_t nBlocks = 19;

    const size_t maxNumThreads = 4;
    std::future<bool> results [ maxNumThreads ];

    for(size_t block=0; block<nBlocks; block++){

        // define block of records to process in this iteration
        size_t begin = readsInBlock*(block + 0);
        size_t end   = readsInBlock*(block + 1);
        if( end > nReads ) end = nReads;
        cout<<"Block: "<<block<<" ["<<begin<<" - "<<end<<")"<<endl;

        // identify a free thread
        size_t freeThread = 0;
        for( ;  results[freeThread].valid() &&
                results[freeThread].wait_for(std::chrono::milliseconds(100)) != std::future_status::ready ;
                freeThread++ )
            if( freeThread == maxNumThreads-1 ) freeThread = 0;

        // submit
        results[freeThread] = std::async(std::launch::async, processReads, begin, end, block);
    }

    // wait untill all threads finish
    for(size_t thr=0; thr<maxNumThreads; thr++)
        if( results[thr].valid() ) results[thr].wait();

    // report number of clusters in each block
    for(size_t block=0; block<nBlocks; block++)
        cout<<"Found "<<uf[block]->nClusters()<<" clusters in block #"<<block<<endl;

    // group clusters together 
    list<int> parentNodes;
    for(size_t block=0; block<nBlocks; block++){
        const map<int, list<int> > &clusters = uf[block]->clusters();
        for(map<int, list<int> >::const_iterator clust = clusters.begin(); clust != clusters.end(); clust++)
            parentNodes.push_back(clust->first);
    }

    UnionFind superClusters(parentNodes);

    map<int,LookUpTable> searchHelper[nBlocks];

    for(size_t block=0; block<nBlocks; block++){
        searchHelper[block] = buildSearchHelper(*(uf[block]));
        cout<<"SearchHelper for "<<block<<" is generated"<<endl;
    }

    for(size_t block=0; block+1<nBlocks; block++){
        for(map<int,LookUpTable>::const_iterator iter = searchHelper[block].begin(); iter != searchHelper[block].end(); iter++){

            for(size_t block2=block+1; block2<nBlocks; block2++){

                for(map<int,LookUpTable>::const_iterator iter2 = searchHelper[block2].begin(); iter2 != searchHelper[block2].end(); iter2++){

                    if( iter->second.match( iter2->second ) ){
                        int cluster1 = superClusters.findCluster(iter ->first);
                        int cluster2 = superClusters.findCluster(iter2->first);
                        if( cluster1 != cluster2 ) superClusters.joinClusters(cluster1,cluster2);
                    }

                }
            }
//            cout<<iter->first<<" done "<<endl;
        }
        cout<<" Block "<<block<<" done"<<endl;
    }

    cout<<"Found "<<superClusters.nClusters()<<" superclusters "<<endl;


    ofstream output("output.csv");
    if( !output ){ cout<<"Cannot open "<<"output.csv"<<endl; return 0; }
    const map<int, list<int> > &q = superClusters.clusters();
    for(map<int, list<int> >::const_iterator iter = q.begin(); iter != q.end(); iter++){
        int seed = iter->first;
        size_t nClusters = 0;
        for(list<int>::const_iterator clust = iter->second.begin(); clust != iter->second.end(); clust++){
            size_t block = (*clust-1)/readsInBlock;
            const map<int, list<int> > &w = uf[block]->clusters();
            map<int, list<int> >::const_iterator e = w.find(*clust);
            if( e != w.end() )
                nClusters += e->second.size();
            else { cerr<<"Problem"<<endl; exit(0); }
        }
        if( nClusters == 0 ) cerr<<" Error: empty cluster for "<<seed<<"!"<<endl;
        output<<seed-1<<","<<nClusters<<endl;
//        if(seed==1)
//        for(list<int>::const_iterator node = iter->second.begin(); node != iter->second.end(); node++)
//            output<<seed-1<<","<<*node-1<<endl;
    }
    output.close();

/*
    ofstream ungrouped("ungrouped.fastq");
    if( !ungrouped ){ cout<<"Cannot open ungrouped.fastq"<<endl; return 0; }

    for(map<int, list<int> >::const_iterator iter = clusters.begin(); iter != clusters.end(); iter++){
        int seed = iter->first;
        if( iter->second.size() == 0 ) cerr<<" Error: empty cluster for "<<seed<<"!"<<endl;

        if( iter->second.size() < 100 ){
            for(list<int>::const_iterator node = iter->second.begin(); node != iter->second.end(); node++){
                ungrouped<<identifier[*node-1]<<endl;
                ungrouped<<sequence  [*node-1]<<endl;
                ungrouped<<"+"<<endl;
                ungrouped<<quality   [*node-1]<<endl;
            }
            continue;
        }

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

    ungrouped.close();
*/
    return 0;
}
