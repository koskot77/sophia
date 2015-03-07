#ifndef TOOLBOX_H
#define TOOLBOX_H
#include <string.h>  // bzero,strlen 

// Maximum length of the adaptor name and the coding sequence
#define MAX_ADAPTORS (128)
#define MAX_LENGTH   (128)

// A collisionless hash-like helper function to convert a sequence of bases (limited to 32 symbols) into an integer number 
//   in case of a problem (e.g. non-interpretable sequence) last arguments returns position of the error (starting from 1)
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

    // calculating the code number
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



// Let us build a classical hash function around '%' operator and construct a look-up table
#define BUCKETS (104743)      // a moderate size prime number
#define MAX_COLLISIONS (10)   // allow up to 10 collisions
class LookUpTable {
private:
    unsigned long long values [MAX_ADAPTORS+1];  // original values (converted keys); indexing starts from 1 
    size_t             lengths[MAX_ADAPTORS+1];  // we should also carry around the length of original sequence in case it was ending with 'T's
    size_t **table; //[BUCKETS][MAX_COLLISIONS]; // look-up table; stores 0 for empty buckets or an index of a potential match for a given hash value
    size_t  *nCollisions; //[BUCKETS];           // number of collisions (should mostly be 0)

public:
    // performance
    size_t countCollisions(void) const {
        size_t retval = 0;
        for(size_t bucket=0; bucket<BUCKETS; bucket++) retval += nCollisions[bucket];
        return retval;
    }

    // quick search function returns MAX_ADAPTORS if requested number is not among the known keys, otherwise key index (start from 0)
    size_t find(unsigned long long number, size_t &length) const {
        // following two lines is as far as we get most of the times
        size_t bucket = number % BUCKETS;
        if( table[bucket][0] == 0 ) return MAX_ADAPTORS;

        // or we may have found a potential match; check if it indeed the case
        size_t matchIndex = MAX_ADAPTORS;

        // look for complete or partial match
        size_t bestMatchLength = 0;
        for(size_t col=0; col <= nCollisions[bucket] && matchIndex == MAX_ADAPTORS; col++){
            size_t             ind = table[bucket][col];
            unsigned long long val = values [ind];
            size_t             len = lengths[ind];
            if( (number & ((0x1LL<<(len*2))-1)) == val && bestMatchLength < len ){
                bestMatchLength = len;
                matchIndex = ind;
            }
        }

        // if match is found make the index starting from 0
        if( matchIndex != MAX_ADAPTORS ) matchIndex--; 

        length = bestMatchLength;

        return matchIndex;
    }

    // construct the table from the set of keys (key's index serves the value)
    LookUpTable(const char keys[MAX_ADAPTORS][MAX_LENGTH]){
        // keys may have arbitrary lengths, but the quick search function operates on the fixed length numbers
        //  this is achieved with padding of the shorter keys to the maximal key length
        //  in the process the shorter keys are multiplexed to all possible combinations of padding symbols

        // initialize static members
        bzero(values,  sizeof(values));
        bzero(lengths, sizeof(lengths));

        // find longest key
        size_t maxKeyLength = 0;
        for(size_t k=0; k<MAX_ADAPTORS; k++){
            size_t length = strlen( keys[k] );
            if( length > maxKeyLength ) maxKeyLength = length;
        }

        // allocate and initialize dynamic tables
        nCollisions = new size_t [BUCKETS];
        bzero(nCollisions, sizeof(size_t)*BUCKETS);

        table = new size_t* [BUCKETS];
        for(size_t k=0; k<BUCKETS; k++){
            table[k] = new size_t [MAX_COLLISIONS];
            bzero(table[k], sizeof(size_t)*MAX_COLLISIONS);
        }

        // build hash values for the padded keys
        for(size_t k=0; k<MAX_ADAPTORS; k++){
            size_t length = strlen( keys[k] );
            if( length > 0 ){
                unsigned short errorPos = 0;
                unsigned long long coreCode = sequence2number(keys[k],length,errorPos);
                // skip any non-interpretable key
                if( errorPos ) continue;
                // remember the sequence in numeric form
                lengths[k+1] = length;
                values [k+1] = coreCode; 

                // now create all new padding combinations and add those to the table
                for(size_t padCode=0; padCode<(0x1<<((maxKeyLength - length)*2)); padCode++){
                    unsigned long long code = coreCode | (padCode<<(length*2));
                    size_t bucket = code % BUCKETS;
                    size_t collision = 0;
                    while( table[bucket][collision] ) collision++;
                    table[bucket][collision] = k+1;
                    nCollisions[bucket] = collision;
                }
            }
        }
    }

    ~LookUpTable(void){
        for(size_t k=0; k<MAX_ADAPTORS; k++)
            delete [] table[k];
        delete [] table;
        delete [] nCollisions;
    }
};

#undef MAX_COLLISIONS
#undef BUCKETS

//#undef MAX_ADAPTORS
//#undef MAX_LENGTH

#endif
