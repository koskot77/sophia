#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#define MIN(A,B,C) ( A<B ? ( B<C ? A : ( A<C ? A : C ) ) : ( B<C ? B : C ) )
// penalties
const size_t gapCost =  50/2;  // gap
const size_t misCost = (5+8); // mismatch

// Needleman-Wunsch score, the 2D score matrix should already be allocated for [len(seq1)+1][len(seq2)+1] dimentions
size_t alignmentScoreMatrix(const char *seq1, const char *seq2, size_t **score){
    // initialization
    size_t len1 = strlen(seq1);
    size_t len2 = strlen(seq2);
    size_t (*A)[len2+1] = (size_t (*)[len2+1]) score;

    // boundary conditions
    for(size_t i=0; i<=len1; i++) A[i][0] = i*gapCost;
    for(size_t j=0; j<=len2; j++) A[0][j] = j*gapCost;

    // the score matrix
    for(size_t i=1; i<=len1; i++)
        for(size_t j=1; j<=len2; j++)
            A[i][j] = MIN( 
                          A[i-1][j-1] + ( seq1[i-1]==seq2[j-1] ? 0 : misCost ),
                          A[i-1][j]   + gapCost,
                          A[i][j-1]   + gapCost
                         );

     // total score
     return A[len1][len2];
}
#undef MIN

// reconstruct the alignments, the x and y should be both allocated for [len(seq1)+len(seq2)] size
void reconstruction(const char *seq1, const char *seq2, const size_t **score, char *x, char *y){
    // initialization
    size_t len1 = strlen(seq1);
    size_t len2 = strlen(seq2);
    size_t (*A)[len2+1] = (size_t (*)[len2+1]) score;

    // go back over the score matrix
    size_t k = 0, i = len1, j = len2;
    while( i!=0 && j!=0 ){
        if( A[i][j] == A[i-1][j-1] + ( seq1[i-1]==seq2[j-1] ? 0 : misCost ) ){
            x[k] = seq1[i-1];
            y[k] = seq2[j-1];
            i--;
            j--;
        } else
        if( A[i][j] == A[i-1][j] + gapCost ){
            x[k] = seq1[i-1];
            y[k] = '-';
            i--;
        } else
        if( A[i][j] == A[i][j-1] + gapCost ){
            x[k] = '-';
            y[k] = seq2[j-1];
            j--;
        }
        k++;
    }
    while( i==0 && j>0 ){ x[k] = '-'; y[k] = seq2[j-1]; k++; j--; }
    while( j==0 && i>0 ){ y[k] = '-'; x[k] = seq1[i-1]; k++; i--; }

    // reverse the sequences
    for(size_t i=0; i<k/2; i++){
        char tmp = x[k-1-i];
        x[k-1-i] = x[i];
        x[i]     = tmp;
    }
    for(size_t i=0; i<k/2; i++){
        char tmp = y[k-1-i];
        y[k-1-i] = y[i];
        y[i]     = tmp;
    }

    x[k] = '\0';
    y[k] = '\0';
}

// A collisionless hash-like helper function to convert a sequence of bases (limited to 32 symbols) into an integer number 
//   in case of a problem (e.g. non-interpretable sequence) last arguments returns position of the error (starting at 1)
unsigned long long sequence2number(const char *sequence, unsigned short length, unsigned short &errPos){
    if( length > 32 ) return 0;
    // ascii codes of the symbols:
    //  T:  84, G: 71, A:65, C:67
    //  t: 116, g:103, a:97, c:99
    // construct a look-up table: symbolic code -> digit
    // reserve code '4' to find out if the sequence has anything else besides precisely the four expected bases (i.e. T,G,A,C)
    const static unsigned short ascii2digit[128] = {
        // assign each symbol with an unique number from 0 to 4
        #define T 0
        #define G 1
        #define A 2
        #define C 3
        #define t 0
        #define g 1
        #define a 2
        #define c 3
        4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,A,4,C,4,4, 4,G,4,4,4,4,4,4,4,4,
        4,4,4,4,T,4,4,4,4,4, 4,4,4,4,4,4,4,a,4,c,
        4,4,4,g,4,4,4,4,4,4, 4,4,4,4,4,4,t,4,4,4,
        4,4,4,4,4,4,4,4
        #undef c
        #undef a
        #undef g
        #undef t
        #undef C
        #undef A
        #undef G
        #undef T
    };

    // calculating the number
    unsigned long long retval = 0;
    for(size_t pos=0,order=0; pos<length; pos++,order+=2){
        unsigned short code = ascii2digit[(unsigned short)(sequence[pos])];
        if( code>3 ){       // unknown base?
            errPos = pos+1; // raise an error
            code   = 0;     // pretend it is 'T' and carry on
        }
        retval += (unsigned long long)(code) << order;
    }

    // beware: sequence with 'T' symbol(s) in the end converted to the same number as sequence without these 'T's
    //  hence, always keep track of the sequence's length when use result of this function
    return retval;
}

// reverse function
const char* number2sequence(unsigned long long number, unsigned short length){

    if( length > 32 ) return 0;

    static char buffer[32];
    bzero(buffer,sizeof(buffer));

    // construct the look-up table digit -> symbolic code
    const static char digit2ascii[4] = {'T','G','A','C'};

    // decoding the number into a sequence
    for(size_t order=0,pos=0; pos<length; pos++,order+=2)
        buffer[pos] = digit2ascii[ (number>>order)&0x3 ];

    return buffer;
}

// Helper class to scan over the sequence in numeric representation
class NumericSequence {
private:
    unsigned long long *data;   // numeric representation of the sequence
    size_t size, length;        // size of the array and length of the original symbolic sequence
    const static unsigned symbolsInOneElement; // number of symbols coded by a single element of the numeric array
    unsigned long errorPos;    // first occurrence position of an interpretation problem (if any -> starting from 1)

public:
    unsigned long error(void) const { return errorPos; }

    // random access view; start and width are measured in symbols
    unsigned long long view(size_t start, size_t width) const {
        unsigned long long retval = 0;
        // check boundaries for start and width arguments
        if( start >= length || width > symbolsInOneElement ) return retval;
        // identify element in the array, relevant part of the element, and spill over to the next element
        size_t block = start / symbolsInOneElement;
        size_t index = start % symbolsInOneElement;
        long   nLeft = index + width - symbolsInOneElement;
        // select the codes
        retval = data[block] >> (index*2);
        // see if we need to take the rest from next element
        if( nLeft>0 )
            retval |= (data[block+1]&((0x1LL<<(nLeft*2))-1)) << ((symbolsInOneElement-index)*2);
        else
            retval &= (0x1LL<<(width*2)) - 1;

        return retval;
    }

    // construct numeric sequence from a symbolic sequence
    NumericSequence& operator=(const char *symbolicSequence){
        // reset errors
        unsigned short err = 0;
        errorPos = 0;
        // destroy the previous numeric sequence (if any)
        if( data ) delete [] data;
        // allocate new data array for the numeric sequence
        length = strlen(symbolicSequence);
        size   = length/symbolsInOneElement + 1;
        data   = new unsigned long long [size+1];
        // convert the symbolic sequence into the numeric sequence and store it in the array
        const char *ptr = symbolicSequence;
        for(size_t block = 0; block < size-1; block++){
            data[block] = sequence2number(ptr,symbolsInOneElement,err);
            if( err!=0 && errorPos==0 ) errorPos = err + block * symbolsInOneElement;
            ptr += symbolsInOneElement;
        }
        data[size-1] = sequence2number(ptr, length - (ptr-symbolicSequence), err);
        if( err!=0 && errorPos==0 ) errorPos = err + (size-1) * symbolsInOneElement;

        data[size] = 0;

        return *this;
    }

    // construct numeric sequence from a symbolic sequence
    NumericSequence(const char *symbolicSequence):data(0),size(0),length(0){
        this->operator=(symbolicSequence);
    }
    // copying constructor is a "must have thing" whenever objects owns dynamically allocated data
    NumericSequence(const NumericSequence& src){
        length = src.length;
        size   = src.size;
        data   = new unsigned long long [size];
        memcpy( data, src.data, sizeof(unsigned long long)*size );
    }
    // clean up
    ~NumericSequence(void){ delete [] data; }
};
const unsigned NumericSequence::symbolsInOneElement = sizeof(unsigned long long)*4;

#include<map>
#include<set>

// Finally, implementation of the problem as required by the competition
class DNASequencing {
private:
    string reference[25]; // not sure if 24 chromatids Ids start at 0 or 1, let's assume 1
    map< unsigned long long, set<unsigned int> > lookUp[25]; // k-mer -> location (chromotid,positions)
    size_t width, step;   // k-mer size and k-mer step

    const static size_t len = 150;

    size_t alignFast    (size_t chId, size_t refPos, size_t readPos, const string &read, size_t &frist, size_t &last);
    size_t alignAccurate(size_t chId, size_t refPos, size_t readPos, const string &read, size_t &frist, size_t &last);

public:
    int initTest(int testDifficulty){
        switch( testDifficulty ){
            case 0: step = 1; break;
            case 1: step = width; break; // step = width/3;
            case 2: step = width; break;
            default: break;
        }
        return 0;
    }

    int passReferenceGenome(int chromatidSequenceId, const vector<string> &chromatidSequence);

    int preProcessing(void);

    vector<string> getAlignment(size_t N, double normA, double normS, const vector<string> &readName, const vector<string> &readSequence);

    DNASequencing(void){ width = 30; }
    ~DNASequencing(void){}
};


int DNASequencing::passReferenceGenome(int chromatidSequenceId, const vector<string> &chromatidSequence){
    // Initialization of the 24 reference chromatides
    if( chromatidSequenceId < 0 || chromatidSequenceId > 24 ) return -1;
    for( auto &line : chromatidSequence )
        reference[ chromatidSequenceId ].append( line.substr(0,line.length()-1) ); // Fucking Windows eol extra symbol!
    return 0;
}

int DNASequencing::preProcessing(void){
    // Chop every reference chromatid into k-mers of a certain width for fast look-ups
    for(size_t chId=0; chId<25; chId++){
        if( reference[chId].length() == 0 ) continue;
        const char *seq = reference[chId].c_str();
        NumericSequence numSeq(seq);
//        unsigned short err;
//        unsigned long long view = sequence2number(seq,width,err);
        for(unsigned long long pos=0; pos<reference[chId].length()-width; pos+=step){
//            view = (view >> (step*2)) | ( sequence2number(seq+pos,step,err) << ((width-step)*2) );
            unsigned long long view = numSeq.view(pos,width);
            lookUp[chId][view].insert( pos );
//if( view == 0xb6a4f583ec00f02LL ) cout<<"pos1="<<pos<<endl;
//if( view == 0xb66ec0fcee90a6eLL ) cout<<"pos2="<<pos<<endl;
//if( pos == 28494871 ) cout<<"view1="<<hex<<view<<hex<<endl;
//if( pos == 28496839 ) cout<<"view2="<<hex<<view<<hex<<endl;
        }
        cout<<"chId="<<chId<<" done"<<endl;
    }
    return 0;
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

size_t DNASequencing::alignFast(size_t chId, size_t refPos, size_t readPos, const string &read, size_t &first, size_t &last){
// A fast sloppy version that doesn't aim at very accurate begin-end (first-last) alignment position calculation
    static size_t score[len+1][len+1];
    string ref = reference[chId].substr( refPos-readPos, len );
    size_t s = alignmentScoreMatrix(read.c_str(), ref.c_str(), (size_t**)score);

    first = refPos-readPos;
    last  = len;

    return s;
}

size_t DNASequencing::alignAccurate(size_t chId, size_t refPos, size_t readPos, const string &read, size_t &first, size_t &last){
// A more thoral version that estimates accurate begin-end (first-last) alignment position calculation (seemingly even better then in the validation sets)
    unsigned long long start  = refPos - readPos - (readPos*misCost)/gapCost; // add contingency for potential indels
    if( start  < 0 ) start = 0;

    int length = len + (readPos*misCost)/gapCost + ((len-readPos-width)*misCost)/gapCost; // add contingency for potential indels //((len - width)*misCost)/gapCost; 
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

    return s;
}

vector<string> DNASequencing::getAlignment(size_t N, double normA, double normS, const vector<string> &readName, const vector<string> &readSequence){
    vector<string> retval;

    for(size_t read=0; read<N/2; read++){
        // all reads are paired, always consider them together 

        size_t read1 = 2*read;
        size_t read2 = 2*read + 1;

        string reverseCompliment1;
        for(int pos = readSequence[ read1 ].length()-1; pos >= 0; pos--)
            reverseCompliment1 += complement[ readSequence[ read1 ][ pos ] ];

        string reverseCompliment2;
        for(int pos = readSequence[ read2 ].length()-1; pos >= 0; pos--)
            reverseCompliment2 += complement[ readSequence[ read2 ][ pos ] ];

        // label forward ("F") direction when first sequence match to '+' and second to '-'
        const char *seqF1 = readSequence[read1].c_str();
        const char *seqF2 = reverseCompliment2 .c_str();
        // label reverse ("R") direction when first sequence match to '-' and second to '+'
        const char *seqR1 = reverseCompliment1 .c_str();
        const char *seqR2 = readSequence[read2].c_str();

        NumericSequence numF1( seqF1 );
        NumericSequence numF2( seqF2 );
        NumericSequence numR1( seqR1 );
        NumericSequence numR2( seqR2 );

        // for the correct alignment found for the current k-mer, the next k-mer most likely result in the same alignment
        //  we will skip the costly alignment for the intervals that have already bee aligned earlier
        //  the map below is association if alignment positions: end -> (begin,score)
        map< unsigned int, pair<unsigned int,size_t> > alreadySeenF1[25], alreadySeenF2[25], alreadySeenR1[25], alreadySeenR2[25];

        size_t bestScore = 1000000;
        int    bestCh    = -1;
        bool   revCompl  = false;
        unsigned long long bestBegin1 = 0;
        unsigned long long bestBegin2 = 0;
        unsigned long long bestEnd1   = 0;
        unsigned long long bestEnd2   = 0;

///clock_t start = clock();

        // let's slide along the both reads simultaneously
        for(unsigned long long pos=0; pos<len-width; pos++){

            if( bestScore == 0 ) break; // shortcut for the competition

            // get hash codes for the k-mer
            unsigned long long viewF1 = numF1.view(pos,width);
            unsigned long long viewF2 = numF2.view(pos,width);

            unsigned long long viewR1 = numR1.view(pos,width);
            unsigned long long viewR2 = numR2.view(pos,width);

            // try matching every chromatid
            for(size_t chId=0; chId<25; chId++){

                if( lookUp[chId].size() == 0 ) continue;

                // look for a match with forward direction hypothesis
                map< unsigned long long, set<unsigned int> >::const_iterator hitF1 = lookUp[chId].find(viewF1);
                map< unsigned long long, set<unsigned int> >::const_iterator hitF2 = lookUp[chId].find(viewF2);

                // look for a match with reverse direction hypothesis
                map< unsigned long long, set<unsigned int> >::const_iterator hitR1 = lookUp[chId].find(viewR1);
                map< unsigned long long, set<unsigned int> >::const_iterator hitR2 = lookUp[chId].find(viewR2);

//cout<<"pos = "<<pos<<" viewF1 = "<<hex<<viewF1<<" viewF2 = "<<viewF2<<" viewR1 = "<<viewR1<<" viewR2 = "<<viewR2<<dec<<endl;

                // while belonging to the same DNA fragment, the paired reads cannot be far away
                // check it for the first hypothesis
                if( hitF1 != lookUp[chId].end() && hitF2 != lookUp[chId].end() ){ // both of paired reads fire some k-mers in the same chromatid?

                    // found set(!) of hits belong to the same chromatid, now check if their positions are far apart
                    const set<unsigned int> &biggerList       = ( hitF1->second.size() >  hitF2->second.size() ? hitF1->second : hitF2->second );
                    const set<unsigned int> &smallerList      = ( hitF1->second.size() <= hitF2->second.size() ? hitF1->second : hitF2->second );
                    const string &seqBig                      = ( hitF1->second.size() >  hitF2->second.size() ? readSequence[read1]: reverseCompliment2);
                    const string &seqSmall                    = ( hitF1->second.size() <= hitF2->second.size() ? readSequence[read1]: reverseCompliment2);
                    map< unsigned int, pair<unsigned int,size_t> > &seenBig   = ( hitF1->second.size() > hitF2->second.size() ? alreadySeenF1[chId] : alreadySeenF2[chId] );
                    map< unsigned int, pair<unsigned int,size_t> > &seenSmall = ( hitF1->second.size() < hitF2->second.size() ? alreadySeenF1[chId] : alreadySeenF2[chId] );

                    if( smallerList.size() > 20 ) continue; // suspiciously many hits, probably not very informative k-mer -> skip the whole case

//cout<<" same chomo1, nMatches: "<<smallerList.size()<<" sizeF1: "<<hitF1->second.size()<<" sizeF2: "<<hitF2->second.size()<<endl;

                    for( auto &refPos : smallerList ){

//cout << "refPos1:" << refPos << " wrt: "<< *(biggerList.begin()) << endl;

                        // consider all paired alignments close by within 700 base pairs
                        set<unsigned int>::const_iterator complement = biggerList.upper_bound(int(refPos)-int(700));
                        while( complement != biggerList.end() && int(*complement)-int(refPos) < 700 ){

//cout << "complement1:" << *complement << endl;

                            size_t score1, score2;
                            size_t first1, last1;
                            size_t first2, last2;

                            map< unsigned int, pair<unsigned int,size_t> >::const_iterator candidate = seenSmall.lower_bound(refPos);
                            if( candidate == seenSmall.end() || candidate->second.first > refPos ){
                                score1 = alignAccurate(chId, refPos, pos, seqSmall, first1, last1);
                                seenSmall[last1] = pair<unsigned int,size_t>(first1,score1);
//cout<<"alignment1 score1:"<<score1<<endl;
                            } else
                                score1 = candidate->second.second;

                            map< unsigned int, pair<unsigned int,size_t> >::const_iterator candidate2 = seenBig.lower_bound(*complement);
                            if( candidate2 == seenBig.end() || candidate2->second.first > *complement ){
                                score2 = alignAccurate(chId, *complement, pos, seqBig, first2, last2);
                                seenBig[last2] = pair<unsigned int,size_t>(first2,score2);
//cout<<"alignment1 score2:"<<score2<<endl;
                            } else
                                score2 = candidate2->second.second;

                            if( bestScore > score1 + score2 ){
                                bestScore = score1 + score2;
                                bestCh    = chId;
                                bestBegin1 = first1;
                                bestEnd1   = last1;
                                bestBegin2 = first2;
                                bestEnd2   = last2;
                                revCompl  = false;
                            }

                            complement++;
                        }

//                        if( s == 0 ) break;
                    }
                }

                // check it for the second hypothesis
                if( hitR1 != lookUp[chId].end() && hitR2 != lookUp[chId].end() ){ // both of paired reads fire some k-mers in the same chromatid?

                    // found set(!) of hits belong to the same chromatid, now check if their positions are far apart
                    const set<unsigned int> &biggerList       = ( hitR1->second.size() >  hitR2->second.size() ? hitR1->second : hitR2->second );
                    const set<unsigned int> &smallerList      = ( hitR1->second.size() <= hitR2->second.size() ? hitR1->second : hitR2->second );
                    const string &seqBig                      = ( hitR1->second.size() >  hitR2->second.size() ? reverseCompliment1: readSequence[read2]);
                    const string &seqSmall                    = ( hitR1->second.size() <= hitR2->second.size() ? reverseCompliment1: readSequence[read2]);
                    map< unsigned int, pair<unsigned int,size_t> > &seenBig   = ( hitR1->second.size() > hitR2->second.size() ? alreadySeenR1[chId] : alreadySeenR2[chId] );
                    map< unsigned int, pair<unsigned int,size_t> > &seenSmall = ( hitR1->second.size() < hitR2->second.size() ? alreadySeenR1[chId] : alreadySeenR2[chId] );

//cout<<" same chomo2, nMatches: "<<smallerList.size()<<" sizeR1: "<<hitR1->second.size()<<" sizeR2: "<<hitR2->second.size()<<endl;

                    if( smallerList.size() > 20 ) continue; // suspiciously many hits, probably not very informative k-mer -> skip the whole case

                    for( auto &refPos : smallerList ){

//cout << "refPos2:" << refPos << " wrt: "<< *(biggerList.begin()) << endl;

                        // consider all paired alignments close by within 700 base pairs
                        set<unsigned int>::const_iterator complement = biggerList.upper_bound(int(refPos)-int(700));
                        while( complement != biggerList.end() && int(*complement)-int(refPos) < 700 ){

//cout << "complement2:" << *complement << endl;

                            size_t score1, score2;
                            size_t first1, last1;
                            size_t first2, last2;

                            map< unsigned int, pair<unsigned int,size_t> >::const_iterator candidate = seenSmall.lower_bound(refPos);
                            if( candidate == seenSmall.end() || candidate->second.first > refPos ){
                                score1 = alignAccurate(chId, refPos, pos, seqSmall, first1, last1);
                                seenSmall[last1] = pair<unsigned int,size_t>(first1,score1);
//cout<<"alignment2 score1:"<<score1<<endl;
                            } else
                                score1 = candidate->second.second;

                            map< unsigned int, pair<unsigned int,size_t> >::const_iterator candidate2 = seenBig.lower_bound(*complement);
                            if( candidate2 == seenBig.end() || candidate2->second.first > *complement ){
                                score2 = alignAccurate(chId, *complement, pos, seqBig, first2, last2);
                                seenBig[last2] = pair<unsigned int,size_t>(first2,score2);
//cout<<"alignment2 score2:"<<score2<<endl;
                            } else
                                score2 = candidate2->second.second;

                            if( bestScore > score1 + score2 ){
                                bestScore = score1 + score2;
                                bestCh    = chId;
                                bestBegin1 = first1;
                                bestEnd1   = last1;
                                bestBegin2 = first2;
                                bestEnd2   = last2;
                                revCompl  = true;
                            }

                            complement++;
                        }

//                        if( s == 0 ) break;
                    }
                }


            }
        }

        static char buffer[1024];
        if( bestCh >= 0 ){
if( !revCompl ){ 
            sprintf(buffer,"%s,%d,%lld,%lld,%c,%ld",readName[read1].c_str(),bestCh+1,bestBegin1+1,bestEnd1+1,'+',bestScore);
            retval.push_back( buffer );
///cout<<buffer<<" seq="<<reference[bestCh].substr(bestBegin1,bestEnd1-bestBegin1)<<" vs. "<<readSequence[read1]<<endl;
            sprintf(buffer,"%s,%d,%lld,%lld,%c,%ld",readName[read2].c_str(),bestCh+1,bestBegin2+1,bestEnd2+1,'-', bestScore);
            retval.push_back( buffer );
///cout<<buffer<<" seq="<<reference[bestCh].substr(bestBegin2,bestEnd2-bestBegin2)<<" vs. "<<reverseCompliment2<<endl;
} else {
            sprintf(buffer,"%s,%d,%lld,%lld,%c,%ld",readName[read1].c_str(),bestCh+1,bestBegin1+1,bestEnd1+1,'-',bestScore);
            retval.push_back( buffer );
///cout<<buffer<<" seq="<<reference[bestCh].substr(bestBegin1,bestEnd1-bestBegin1)<<" vs. "<<reverseCompliment1<<endl;
            sprintf(buffer,"%s,%d,%lld,%lld,%c,%ld",readName[read2].c_str(),bestCh+1,bestBegin2+1,bestEnd2+1,'+', bestScore);
            retval.push_back( buffer );
///cout<<buffer<<" seq="<<reference[bestCh].substr(bestBegin2,bestEnd2-bestBegin2)<<" vs. "<<readSequence[read2]<<endl;

}
        } else {
            sprintf(buffer,"%s,%d,%lld,%lld,%c,%ld",readName[read1].c_str(),bestCh+1,bestBegin1+1,bestEnd1+1,(revCompl?'-':'+'), bestScore);
            retval.push_back( buffer );
///cout<<buffer<<" seq=NULL"<<endl;
            sprintf(buffer,"%s,%d,%lld,%lld,%c,%ld",readName[read2].c_str(),bestCh+1,bestBegin2+1,bestEnd2+1,(revCompl?'-':'+'), bestScore);
            retval.push_back( buffer );
///cout<<buffer<<" seq=NULL"<<endl;
        }

///clock_t end = clock() - start;
///cout << "Time spent: " << (double)end / ((double)CLOCKS_PER_SEC) << " s" << endl;;

    }
    return retval;

}

