#include "toolbox.h"
#include "alignment.cc"

class DNASequencing {
private:
    string reference[24];
    map<unsigned long long, vector<unsigned long long> > lookUp[24]; // k-mer -> positions
    size_t width;

    const static size_t len = 150;

    size_t align1(size_t chId, unsigned long long refPos, size_t readPos, const string &read, unsigned long long &frist, unsigned long long &last);

public:
    int initTest(int testDifficulty){ return 0; }
    int passReferenceGenome(int chromatidSequenceId, const vector<string> &chromatidSequence);
    int preProcessing(void);
    vector<string> getAlignment(int N, double normA, double normS, const vector<string> &readName, const vector<string> &readSequence);

    DNASequencing(void){ width=30; }
    ~DNASequencing(void){}
};


int DNASequencing::passReferenceGenome(int chromatidSequenceId, const vector<string> &chromatidSequence){
    if( chromatidSequenceId < 0 || chromatidSequenceId >= 24 ) return -1;
    for( auto &line : chromatidSequence )
        reference[ chromatidSequenceId ].append( line.substr(0,line.length()-1) ); // Fucking Windows eol extra symbol!
    return 0;
}

int DNASequencing::preProcessing(void){
   for(size_t chId=0; chId<24; chId++){
        if( reference[chId].length() == 0 ) continue;
        const char *seq = reference[chId].c_str();
        NumericSequence numSeq(seq);
        for(unsigned long long pos=0; pos<reference[chId].length(); pos+=width){
            unsigned long long view = numSeq.view(pos,width);
            lookUp[chId][view].push_back(pos);
        }
        cout<<"chId="<<chId<<" done"<<endl;
    }
}


const char complement[128] = {
        // 'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'
        'N','N','N','N','N','N','N','N','N','N', 'N','N','N','N','N','N','N','N','N','N',
        'N','N','N','N','N','N','N','N','N','N', 'N','N','N','N','N','N','N','N','N','N',
        'N','N','N','N','N','N','N','N','N','N', 'N','N','N','N','N','N','N','N','N','N',
        'N','N','N','N','N','T','N','G','N','N', 'N','C','N','N','N','N','N','N','N','N',
        'N','N','N','N','A','N','N','N','N','N', 'n','n','n','n','n','n','n','t','n','g',
        'n','n','n','c','n','n','n','n','n','n', 'n','n','n','n','n','n','a','n','n','n',
        'n','n','n','n','n','n','n','n'
};

size_t DNASequencing::align1(size_t chId, unsigned long long refPos, size_t readPos, const string &read, unsigned long long &first, unsigned long long &last){

    unsigned long long start  = refPos - readPos - (readPos*misCost)/gapCost; // add contingency for potential indels
    if( start  < 0 ) start = 0;

    int length = len + (readPos*misCost)/gapCost + ((len-readPos-width)*misCost)/gapCost;
//((len - width)*misCost)/gapCost; // add contingency for potential indels
    if( start  + length > reference[chId].length() )
        length = reference[chId].length() - start;

    string ref = reference[chId].substr( start, length );
    string pad = string( (readPos*misCost)/gapCost, '+' ). append(read). append( ((len-readPos-width)*misCost)/gapCost, '+' );

    size_t score[pad.length()+1][ref.length()+1];
    size_t s = alignmentScoreMatrix(pad.c_str(), ref.c_str(), (size_t**)score);

    static char a[2*len], b[2*len];
    reconstruction(pad.c_str(), ref.c_str(), (const size_t **)score, a, b);

    size_t skipFront = 0, skipRear = 0, la = strlen(a);

    while( skipFront<la && (a[ skipFront ]=='-' || a[ skipFront ]=='+') ){
       s -= ( a[ skipFront ] == '-' || b[ skipFront ] == '-' ? gapCost : misCost);
       skipFront++;
    }
    while( skipRear <la && (a[la-skipRear-1]=='-' || a[la-skipRear-1]=='+') ){
       s -= ( a[ la-skipRear-1 ] == '-' || b[ la-skipRear-1 ] == '-' ? gapCost : misCost);
       skipRear++;
    }

    int refInserts = 0, refDeletes = 0;
    for(size_t pos = skipFront; pos < la - skipRear; pos++){
        if( a[pos] == '-' ) refInserts++;
        if( b[pos] == '-' ) refDeletes++;
    }

    first = start + skipFront;
    last  = start + la - skipRear - 1 - refDeletes;

/* 
a[la-skipRear] = '\0';
b[la-skipRear] = '\0';
cout<<a+skipFront<<endl;
cout<<b+skipFront<<endl<<endl;;
*/

    return s;
}

vector<string> DNASequencing::getAlignment(int N, double normA, double normS, const vector<string> &readName, const vector<string> &readSequence){
    vector<string> retval;

    for(int read=0; read<N; read++){

        string reverseCompliment;
        for(int pos = readSequence[read].length()-1; pos >= 0; pos--)
            reverseCompliment += complement[ readSequence[read][pos] ];

        map<unsigned long long,unsigned long long> alreadySeen[24];

//        const size_t len = 150; //readSequence[read].length();
        const char *seq1 = readSequence[read].c_str();
        const char *seq2 = reverseCompliment .c_str();

        NumericSequence num1( seq1 );
        NumericSequence num2( seq2 );

        size_t bestScore = 1000000;
        int    bestCh    = -1;
        bool   revCompl  = false;
        unsigned long long bestBegin = 0;
        unsigned long long bestEnd   = 0;

        for(unsigned long long pos=0; pos<len-width; pos++){
            unsigned long long view1 = num1.view(pos,width);
            unsigned long long view2 = num2.view(pos,width);

            for(size_t chId=19; chId<20; chId++){
                map<unsigned long long, vector<unsigned long long> >::const_iterator hit1 = lookUp[chId].find(view1);
                map<unsigned long long, vector<unsigned long long> >::const_iterator hit2 = lookUp[chId].find(view2);

                if( hit1 != lookUp[chId].end() ){

                    if( hit1->second.size() > 2 ) continue;

                    for( auto &refPos : hit1->second ){

                        std::map<unsigned long long,unsigned long long>::const_iterator candidate = alreadySeen[chId].lower_bound(refPos);
                        if( candidate != alreadySeen[chId].end() &&
                            candidate->second < refPos ) continue;

                        unsigned long long first, last;

                        size_t s = align1(chId, refPos, pos, readSequence[read], first, last);

                        alreadySeen[chId][last] = first;

                        if( bestScore > s ){
                            bestScore = s;
                            bestCh    = chId;
                            bestBegin = first;
                            bestEnd   = last;
                            revCompl  = false;
                        }

                        if( s == 0 ) break;
                    }
                }

                if( hit2 != lookUp[chId].end() ){

                    if( hit2->second.size() > 2 ) continue;

                    for( auto &refPos : hit2->second ){

                        std::map<unsigned long long,unsigned long long>::const_iterator candidate = alreadySeen[chId].lower_bound(refPos);
                        if( candidate != alreadySeen[chId].end() &&
                            candidate->second < refPos ) continue;

                        unsigned long long first, last;

                        size_t s = align1(chId, refPos, pos, reverseCompliment, first, last);

                        alreadySeen[chId][last] = first;

                        if( bestScore > s ){
                            bestScore = s;
                            bestCh    = chId;
                            bestBegin = first;
                            bestEnd   = last;
                            revCompl  = true;
                        }

                        if( s == 0 ) break;
                    }
                }

            }
        }

        if( bestCh >= 0 ){
            retval.push_back( reference[bestCh].substr( bestBegin, bestEnd - bestBegin ) );

            cout<<readName[read].substr(1)<<","<<(bestCh+1)<<","<<(bestBegin+1)<<","<<(bestEnd+1)<<","<<(revCompl?"-":"+")<<","
                << endl << (revCompl?reverseCompliment:readSequence[read] )
                << endl << reference[bestCh].substr( bestBegin, bestEnd - bestBegin + 1 )
                << endl << bestScore << ", " << bestEnd - bestBegin << endl;
        }

    }
    return retval;

}

