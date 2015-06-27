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

#include <getopt.h>

#include "./toolbox.h"

using namespace std;
// Compile with:
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

// blocks of clusters
#define MAX_BLOCKS (100)
UnionFind* uf[MAX_BLOCKS];
///////////////////////////////////////////////////////////////////////




/////////////// grouping reads with similar beginnings ////////////////

// check first 'barcodeWidth' positions and request identity of 'viewWidth' consecutive symbols
size_t barcodeWidth;// = 12;
size_t viewWidth;//    = 10;

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
        strncpy( matchPattern[ i ], refSeq+i, viewWidth );
        matchPattern[ i ][ viewWidth ] = '\0';
    }

    // build the hash from the patterns
    LookUpTable lookUp(matchPattern);
    // check if it is good
    size_t maxCollisions = 0;
    lookUp.countCollisions( maxCollisions );
    if( maxCollisions >= MAX_COLLISIONS ) return false; // do not tolerate broken LUT

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

        for(size_t pos=0; pos<barcodeWidth-viewWidth; pos++){
            unsigned long long view = numSeq.view(pos,viewWidth);
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
// parent function of the thread
bool processReads(size_t begin, size_t end, size_t block){
     // safety
     if( begin>=end ) return false;
     // safety once again
     if( uf[block] ){
         cout<<"More than one call to process block #"<<block<<endl;
         return false;
     }
     // initialize data structure for clustering
     uf[block] = new UnionFind(begin+1, end+1); //remember to start numbering from 1
     // build clusters
     for(size_t read=begin; read<end; read++)
         if( !groupMatchesFor(read,block,end) ) return false;

     return true;
}
///////////////////////////////////////////////////////////////////////


/////////////// merge clusters from different blocks //////////////////

// for each cluster leader of the block build a look up table 
map<int,LookUpTable> buildSearchHelper(const UnionFind &uf){
    // all clusters for this block
    const map<int, list<int> > &clusters = uf.clusters();

    // association of cluster leader and lut built of all of the search patterns
    map<int,LookUpTable> retval;

    // for each cluster
    for(auto &clust : clusters){

        // cluster leader
        int leader = clust.first - 1; // remember that UnionFind numeration starts from 1
        const list<int> &children = clust.second;

        // safety
        if( sequence[leader].length() < barcodeWidth ) continue;
        if( children.size() == 0 ){
            cerr<<" Error: empty cluster for "<<leader<<"!"<<endl;
            exit(0);
        }

        // construct a set of search patterns
        set<string> patterns;
        // start with the leader 
        for(size_t pos=0; pos<barcodeWidth-viewWidth; pos++){
            if( pos>=MAX_ADAPTORS ){
                cerr<<"Long record ... exiting"<<endl;
                exit(0);
            }
            patterns.insert( sequence[leader].substr(pos,viewWidth) );
        }

        // do the same for every child
        for(auto &child : children){
            // split sequence beginning into search patterns
            for(size_t pos=0; pos<barcodeWidth-viewWidth; pos++){
                if( pos>=MAX_ADAPTORS ){
                    cerr<<"Long record ... exiting"<<endl;
                    exit(0);
                }
                patterns.insert( sequence[child-1].substr(pos,viewWidth) );
            }
        }
        // safety
        if( patterns.size() >= MAX_ADAPTORS ){
            cerr<<"Too many search patterns: "<<patterns.size()<<endl;
            exit(0);
        } 

        // convert patterns into the format accepted by the LookUpTable 
        char matchPattern[MAX_ADAPTORS][MAX_LENGTH];
        bzero(matchPattern, sizeof(matchPattern));
        size_t p=0;
        for(auto &pat : patterns){
            strncpy(matchPattern[p],pat.c_str(),viewWidth);
            matchPattern[p][viewWidth] = '\0';
            p++;
        }

        pair< map<int,LookUpTable>::iterator, bool> res = retval.insert( pair<int,LookUpTable>( leader+1, LookUpTable(matchPattern) ) );
        if( !res.second ){ cerr<<"Duplicate leader ... exiting"<<endl; exit(0); }

        // check quality of the hash table
        size_t maxCollisions = 0;
        res.first->second.countCollisions( maxCollisions );
        if( maxCollisions >= MAX_COLLISIONS ){
//            for(size_t k=0; k<MAX_ADAPTORS; k++){
//                size_t length = strlen( matchPattern[k] );
//                if( length > 0 ){
//                    unsigned short errorPos = 0;
//                    unsigned long long coreCode = sequence2number(matchPattern[k],strlen(matchPattern[k]),errorPos);
//                    printf("%ld: %s, bucket: %lld\n",k,matchPattern[k],coreCode % BUCKETS);
//                }
//            }
            exit(0); // do not tolerate broken LUT
        }
    }

    return retval;
}

// mutex for the data stracture preventing data race
std::mutex mtx;
UnionFind *joinedGroup;

bool mergeGroups(const UnionFind &group1, const UnionFind &group2){

    // build helpers: relatively fast, but costs memory -> do it on demand and destroy after using
    map<int,LookUpTable> searchHelper1 = buildSearchHelper(group1);
    map<int,LookUpTable> searchHelper2 = buildSearchHelper(group2);

    // for every cluster of the first group 
    for(auto &cluster1 : searchHelper1){

        size_t leader1           = cluster1.first;
        const  LookUpTable &lut1 = cluster1.second;

        // compare it with every cluster of the second group
        for(auto &cluster2 : searchHelper2){

            size_t leader2           = cluster2.first;
            const  LookUpTable &lut2 = cluster2.second;

            bool match = false;
            if( lut1.size() > lut2.size() )
                match = lut2.match( lut1 );
            else
                match = lut1.match( lut2 );

            if( match ){
                mtx.lock();
                int cluster1 = joinedGroup->findCluster(leader1);
                int cluster2 = joinedGroup->findCluster(leader2);
                if( cluster1 != cluster2 ) joinedGroup->joinClusters(cluster1,cluster2);
                mtx.unlock();
            }
        }
    }

    return true;
}
///////////////////////////////////////////////////////////////////////




int main(int argc, char *argv[]){

    // parse the options
    static struct option options[] = {
       {"help",         0, 0, 'h'},
       {"input",        1, 0, 'i'},
       {"cores",        1, 0, 'c'},
       {"nreads",       1, 0, 'n'},
       {"threshold",    1, 0, 't'},
       {"width",        1, 0, 'w'},
       {"length",       1, 0, 'l'},
       {0, 0, 0, 0}
    };

    cout<<"Running: ";
    for(int arg=0; arg<argc; arg++) cout<<argv[arg]<<" ";
    cout<<endl;

    // defaults
    char  *fastqFile    = 0;
    size_t readsInBlock = 1000;
    size_t nCores       = 1;
    size_t threshold    = 1000;

    // global defaults
    barcodeWidth = 12; 
    viewWidth    = 10; 

    while( 1 ){
       int index=0;
       int c = getopt_long(argc, argv, "hi:c:n:t:w:l:",options, &index);
       if( c == -1 ) break;
       switch( tolower(c) ) {
           case 'h':
               cout<<"Usage:"<<endl;
               cout<<"-h     ,   --help              show this message"<<endl;
               cout<<"-i     ,   --input             FASTQ input file"<<endl;
               cout<<"-c     ,   --cores             Number of CPU cores"<<endl;
               cout<<"-n     ,   --nreads            Number of reads to process in one block"<<endl;
               cout<<"-t     ,   --threshold         Minimal cluster size to write in a separate file"<<endl;
               cout<<"-w     ,   --width             Number of consecutive matches in a pattern"<<endl;
               cout<<"-l     ,   --length            Search window in the beginning of a sequence"<<endl;
           break;
           case 'i':
               fastqFile = optarg;
           break;
           case 'c':
               nCores = strtoul(optarg,NULL,0);
           break;
           case 'n':
               readsInBlock = strtoul(optarg,NULL,0);
           break;
           case 't':
               threshold = strtoul(optarg,NULL,0);
           break;
           case 'w':
               barcodeWidth = strtoul(optarg,NULL,0);
           break;
           case 'l':
               viewWidth = strtoul(optarg,NULL,0);
           break;
           default : cout<<"Type -h for help"<<endl; return 0;
       }
    }

    if( !fastqFile ) return 0;

    if( readFile(fastqFile) ) return 0;
    else cout<<"Reads: "<<nReads<<endl;

    const size_t nBlocks = nReads/readsInBlock + 1;
    if( nBlocks >= MAX_BLOCKS ){
        cerr<<"Too many blocks requested (>"<<MAX_BLOCKS<<")"<<endl;
        return 0;
    }

    const size_t maxNumThreads = nCores;
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

    // wait until all threads finish
    for(size_t thr=0; thr<maxNumThreads; thr++)
        if( results[thr].valid() ) results[thr].wait();

    // report number of clusters in each block
    for(size_t block=0; block<nBlocks; block++)
        cout<<"Found "<<uf[block]->nClusters()<<" clusters in block #"<<block<<endl;

    // group clusters together 
    list<int> leaders;
    for(size_t block=0; block<nBlocks; block++){
        const map<int, list<int> > &clusters = uf[block]->clusters();
        for(auto clust : clusters)
            leaders.push_back(clust.first);
    }

    joinedGroup = new UnionFind(leaders);

    for(size_t block1=0; block1+1<nBlocks; block1++){
        cout<<"Grouping block "<<block1<<endl;
        for(size_t block2=block1+1; block2<nBlocks; block2++){

            // identify a free thread
            size_t freeThread = 0;
            for( ;  results[freeThread].valid() &&
                    results[freeThread].wait_for(std::chrono::milliseconds(100)) != std::future_status::ready ;
                    freeThread++ )
                if( freeThread == maxNumThreads-1 ) freeThread = 0;

            cout<<" with block "<<block2<<endl;

            // submit
            results[freeThread] = std::async(std::launch::async, mergeGroups, *uf[block1], *uf[block2]);
        }
    }

    // wait until all threads finish
    for(size_t thr=0; thr<maxNumThreads; thr++)
        if( results[thr].valid() ) results[thr].wait();

    cout<<"Found "<<joinedGroup->nClusters()<<" superclusters "<<endl;


    ofstream output("clustering.csv");
    if( !output ){ cout<<"Cannot open clustering.csv"<<endl; return 0; }

    ofstream ungrouped("ungrouped.fastq");
    if( !ungrouped ){ cout<<"Cannot open ungrouped.fastq"<<endl; return 0; }

    const map<int, list<int> > &superClusters = joinedGroup->clusters();
    for(auto &sc : superClusters){
        int leader = sc.first-1;
        list<int> children;

        for(auto clust : sc.second){
            size_t block = (clust-1)/readsInBlock;
            map<int, list<int> >::const_iterator c = uf[block]->clusters().find(clust);
            if( c != uf[block]->clusters().end() )
                children.insert( children.end(), c->second.begin(), c->second.end() );
            else {
                cerr<<"Problem"<<endl;
                exit(0);
            }
        }
        if( children.size() == 0 ) cerr<<" Error: empty cluster for "<<leader+1<<"!"<<endl;
        output<<leader<<","<<children.size()<<endl;

        if( children.size() < threshold ){
            for(auto read : children ){
                ungrouped<<identifier[read-1]<<endl;
                ungrouped<<sequence  [read-1]<<endl;
                ungrouped<<"+"<<endl;
                ungrouped<<quality   [read-1]<<endl;
            }
        } else {

            stringstream fname;
            fname<<"output"<<leader<<".fastq";
            ofstream out(fname.str());
            if( !out ){ cerr<<"Cannot open "<<fname.str()<<endl; return 0; }

            for(auto read : children){
                out<<identifier[read-1]<<endl;
                out<<sequence  [read-1]<<endl;
                out<<"+"<<endl;
                out<<quality   [read-1]<<endl;
            }
        }
    }
    output.close();
    ungrouped.close();

    return 0;
}
