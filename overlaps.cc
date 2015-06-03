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

#include "./toolbox.h"
#include "./alignment.cc"

const char *fileName = "./testReads.fastq";

#define NREC (10000)
size_t nReads = 0;

// g++ -g -Wall -std=c++0x -o q overlaps.cc -lpthread

using namespace std;
vector<string> identifier(NREC);
vector<string> sequence  (NREC);
vector<string> quality   (NREC);

/////////////////// multithreaded calculations /////////////////
const size_t viewWidth = 31;

bool writeProtect = false;
map< size_t,set<size_t> > similars;

#include <pthread.h>

#define handle_error_en(en, msg) \
        do { errno = en; perror(msg); exit(EXIT_FAILURE); } while (0)
#define handle_error(msg) \
        do { perror(msg); exit(EXIT_FAILURE); } while (0)


struct thread_arg {
    pthread_t thread_id;        // ID returned by pthread_create()
    int       thread_num;       // Application-defined thread #
    size_t    readBegin;
    size_t    readEnd;
};

static void* findMatches(void *arg){
    struct thread_arg *a = (thread_arg*)arg;

    printf("findMatches thread %d: readBegin=%ld readEnd=%ld\n",a->thread_num, a->readBegin, a->readEnd);

    // first pass: identify overlapping records
    for(size_t read1=a->readBegin; read1<a->readEnd; read1++){
        // init the map
        similars[read1].clear();

        const char *seq1 = sequence[read1].c_str();

        // break the read into patterns of viewWidth size
        char matchPattern[MAX_ADAPTORS][MAX_LENGTH];
        bzero(matchPattern, sizeof(matchPattern));

        if( strlen(seq1) >= viewWidth ){
            for(size_t i=0; i<strlen(seq1)-viewWidth; i++){
                if( i>=MAX_ADAPTORS ){ cout<<"Long record ... exiting"<<endl; return 0; }
                strncpy(matchPattern[i],seq1+i,viewWidth);
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
                return arg;
            }

            // look for other reads matching any of the patterns of the reference record
            for(size_t read2=read1+1; read2<nReads; read2++){
                const char *seq2 = sequence[read2].c_str();

                NumericSequence numSeq(seq2);
                for(size_t i=0; i<strlen(seq2)-viewWidth; i++){
                    unsigned long long view = numSeq.view(i,viewWidth);
                    size_t matchLength = viewWidth;
                    if( lookUp.find(view,matchLength) != MAX_ADAPTORS ) // found a match
                        do {
                            if( !writeProtect ){
                                writeProtect = true;
                                similars[read1].insert(read2);
                                writeProtect = false;
                            }
                        } while( writeProtect );
                }
            }
        }
        cout<<"Read "<<read1<<" done"<<endl;
    }

    return 0;
}


static void* lookCloser(void *arg){
    struct thread_arg *a = (thread_arg*)arg;

    printf("lookCloser thread %d: readBegin=%ld readEnd=%ld\n",a->thread_num, a->readBegin, a->readEnd);

    for(size_t read1=a->readBegin; read1<a->readEnd; read1++){

        stringstream fname;
        fname<<"output"<<read1<<".csv";
        ofstream output( fname.str().c_str() );
        if( !output ){ cout<<"Cannot open "<<fname<<endl; return 0; }
        output<<"read,score,len1,len2,overlap,overlapScore"<<endl;

        const char *seq1 = sequence[read1].c_str();

        set<size_t> &candidates = similars[read1];
        for(set<size_t>::const_iterator it=candidates.begin(); it!=candidates.end(); it++){

            size_t read2 = *it;

            const char *seq2 = sequence[read2].c_str();

            size_t len1 = strlen(seq1);
            size_t len2 = strlen(seq2);

            size_t score[len1+1][len2+1];
            bzero(score,sizeof(score));

            size_t s = alignmentScoreMatrix(seq1, seq2, (size_t**)score);

            char a[ len1+len2 ], b[ len1+len2 ];
            bzero(a,sizeof(a));
            bzero(b,sizeof(b));

            reconstruction(seq1, seq2, (const size_t **)score, a, b);

            size_t overlapBegins = 0, overlapEnds = strlen(a);
            while( a[overlapBegins]=='-' || b[overlapBegins]=='-' ) overlapBegins++;
            while( a[overlapEnds]  =='-' || b[overlapEnds]  =='-' ) overlapEnds--;

            size_t overlap = overlapEnds - overlapBegins;
            size_t overlapScore = 0;
            for(size_t i=overlapBegins; i<overlapEnds; i++) overlapScore += (a[i]==b[i]?0:1);

            output<<read1<<","<<read2<<","<<s<<","<<len1<<","<<len2<<","<<overlap<<"," <<overlapScore<<endl;
        }

        cout<<"Read "<<read1<<" done"<<endl;

        output.close();
    }

    return 0;
}

int main(int argc, char *argv[]){
    // open input file:
    size_t iteration = atol(argv[1]);
    const size_t nCycles = 100;

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

////////////////////////////
    const size_t numThreads = 4;

    pthread_attr_t attr;

    int s = pthread_attr_init(&attr);
    if( s != 0 ) handle_error_en(s, "pthread_attr_init");

    struct thread_arg *arg = (thread_arg*)calloc(numThreads, sizeof(struct thread_arg));
    if( arg == NULL ) handle_error("calloc");

    for(size_t tnum=0; tnum<numThreads; tnum++){
        arg[tnum].thread_num = tnum + 1;
        arg[tnum].readBegin = nCycles*(tnum + iteration);
        arg[tnum].readEnd   = arg[tnum].readBegin + nCycles;

        s = pthread_create(&arg[tnum].thread_id, &attr, &findMatches, &arg[tnum]);
        if( s != 0 ) handle_error_en(s, "pthread_create");
    }

    s = pthread_attr_destroy(&attr);
    if( s != 0 ) handle_error_en(s, "pthread_attr_destroy");

    // Now join with each thread, and display its returned value
    for(size_t tnum=0; tnum<numThreads; tnum++){
        void *res;
        s = pthread_join(arg[tnum].thread_id, &res);
        if( s != 0 ) handle_error_en(s, "pthread_join");
        if( res != 0 ){ printf("Bad quality hash table ... exiting\n"); return 0; }

        printf("Joined with thread %d\n", arg[tnum].thread_num);
    }



    ofstream output("output.csv");
    if( !output ){ cout<<"Cannot open "<<"output.csv"<<endl; return 0; }
    output<<"read,nMatches"<<endl;
    for(size_t read1=0; read1<nCycles*numThreads; read1++)
        output<<read1<<","<<similars[read1].size()<<endl;
    output.close();



    s = pthread_attr_init(&attr);
    if( s != 0 ) handle_error_en(s, "pthread_attr_init");

    for(size_t tnum=0; tnum<numThreads; tnum++){
        arg[tnum].thread_num = tnum + 1;
        arg[tnum].readBegin = nCycles*(tnum + iteration);
        arg[tnum].readEnd   = arg[tnum].readBegin + nCycles;

        s = pthread_create(&arg[tnum].thread_id, &attr, &lookCloser, &arg[tnum]);
        if( s != 0 ) handle_error_en(s, "pthread_create");
    }

    s = pthread_attr_destroy(&attr);
    if( s != 0 ) handle_error_en(s, "pthread_attr_destroy");

    // Now join with each thread, and display its returned value
    for(size_t tnum=0; tnum<numThreads; tnum++){
        void *res;
        s = pthread_join(arg[tnum].thread_id, &res);
        if( s != 0 ) handle_error_en(s, "pthread_join");

        printf("Joined with thread %d\n", arg[tnum].thread_num);
    }



    free(arg);
////////////////////////////

    return 0;
}
