#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// gcc -D__GENOME__ -g -I. -c splitter.cc && gcc -g -o test splitter.o sff2fastq/sff.genome.o -lstdc++
// time for i in `cat ionXpress_barcode.txt | awk '{print $2}'`; do echo $i; grep $i data.txt | wc; done

// Input file names
const char *inputFileName   = "input.sff";
const char *adaptorFileName = "ionXpress_barcode.txt";

// All the "fixed length" arrays and odd hashes could have been avoided within the STL framework,
//  but I don't know if the libstdc++ is installed on the testing machine so I don't use STL

// Maximum length of the adaptor name and the coding sequence
#define MAX_ADAPTORS (128)
#define MAX_LENGTH   (128)

// Read specified text file with adaptors' names and coding sequences separated by tabs
unsigned int readAdaptors(const char *fileName, char names[MAX_ADAPTORS][MAX_LENGTH], char bases[MAX_ADAPTORS][MAX_LENGTH]){

    FILE *file = NULL;
    if( (file = fopen(fileName,"r")) == NULL ){
        printf("No adaptor file found\n");
        return 0;
    }

    // safety first: make sure we won't get an overflow while scanning the components 
    char format[128];
    sprintf(format,"%%%ds\t%%%ds\n",MAX_ADAPTORS-1,MAX_ADAPTORS-1); 

    unsigned int num = 0;
    while( !feof(file) && num<MAX_ADAPTORS && 
            fscanf(file, format, names[num], bases[num])==2 ) num++;

    fclose(file);

    if( num == MAX_ADAPTORS )
        printf("Warning: only first %d adaptors are to be processed\n",num);

    return num;
}



// A colisionless hash-like helper function to convert a sequence of bases (limited to 32 symbols) into an integer number 
unsigned long long sequence2number(const char *sequence, unsigned short length){
    if( length > 32 ) return 0;
    // ascii codes of the symbols:
    //  T:  84, G: 71, A:65, C:67
    //  t: 116, g:103, a:97, c:99
    // construct the look-up table symbolic code -> a digit
    const static unsigned short ascii2digit[128] = {
        // assign each symbol with an unique number from 0 to 3
        #define T 0
        #define G 1
        #define A 2
        #define C 3
        #define t 0
        #define g 1
        #define a 2
        #define c 3
        0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,A,0,C,0,0, 0,G,0,0,0,0,0,0,0,0,
        0,0,0,0,T,0,0,0,0,0, 0,0,0,0,0,0,0,a,0,c,
        0,0,0,g,0,0,0,0,0,0, 0,0,0,0,0,0,t,0,0,0,
        0,0,0,0,0,0,0,0
        #undef c
        #undef a
        #undef g
        #undef t
        #undef C
        #undef A
        #undef G
        #undef T
    };
    // remark about safety: if ever we encounter an unknown (not one of the four) symbol, it'll be considered symbol 'T' 
    //  although this wouldn't make much sence, we would at least never get an overflow

    // calculating the code number
    unsigned long long retval = 0;
    for(size_t pos=0,order=0; pos<length; pos++,order+=2)
        retval += (unsigned long long)( ascii2digit[(unsigned short)(sequence[pos])] ) << order;

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

public:
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
            retval |= (data[block+1]&(((unsigned long long)(0x1)<<(nLeft*2))-1)) << ((symbolsInOneElement-index)*2);
        else
            retval &= ((unsigned long long)(0x1)<<(width*2)) - 1;

        return retval;
    }

    // construct numeric sequence from a symbolic sequence
    NumericSequence(const char *symbolicSequence){
        // allocate data array for the numeric sequence
        length = strlen(symbolicSequence);
        size   = length/symbolsInOneElement + 1;
        data   = new unsigned long long [size+1];
        // convert the symbolic sequence into the numeric sequence and store it in the array
        const char *ptr = symbolicSequence;
        for(size_t block = 0; block < size-1; block++){
            data[block] = sequence2number(ptr,symbolsInOneElement);
            ptr += symbolsInOneElement;
        }
        data[size-1] = sequence2number(ptr, length - (ptr-symbolicSequence));
        data[size]   = 0;
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
    unsigned long long values [MAX_ADAPTORS+1];  // original values; indexing starts from 1 
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

    // quick search function returns MAX_ADAPTORS if requested number is not among the original values
    size_t find(unsigned long long number, size_t length) const {
        // following two lines is as far as we get most of the times
        size_t bucket = number % BUCKETS;
        if( table[bucket][0] == 0 ) return MAX_ADAPTORS;

        // or we may have found a potential match; check if it indeed the case
        size_t matchIndex = MAX_ADAPTORS;

        for(size_t col=0; col <= nCollisions[bucket] && matchIndex == MAX_ADAPTORS; col++)
            if( number == values [ table[bucket][col] ] &&
                length == lengths[ table[bucket][col] ] )
                matchIndex = table[bucket][col];

        // if match is found make the index starting from 0
        if( matchIndex != MAX_ADAPTORS ) matchIndex--; 

        return matchIndex;
    }

    LookUpTable(const char keys[MAX_ADAPTORS][MAX_LENGTH]){
        // initialize static members
        bzero(values,  sizeof(values));
        bzero(lengths, sizeof(lengths));

        // allocate and initialize dynamic tables
        nCollisions = new size_t [BUCKETS];
        bzero(nCollisions, sizeof(size_t)*BUCKETS);

        table = new size_t* [BUCKETS];
        for(size_t k=0; k<BUCKETS; k++){
            table[k] = new size_t [MAX_COLLISIONS];
            bzero(table[k], sizeof(size_t)*MAX_COLLISIONS);
        }

        // Build hash values for the adaptors
        for(size_t k=0; k<MAX_ADAPTORS; k++){
            size_t length = strlen( keys[k] );
            if( length > 0 ){
                lengths[k+1] = length;
                values [k+1] = sequence2number(keys[k],length);
                size_t bucket = values[k+1] % BUCKETS;
                size_t collision = 0;
                while( table[bucket][collision] ) collision++;
                table[bucket][collision] = k+1;
                nCollisions[bucket] = collision;
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
/*
#include <endian.h>
#include <stdint.h>

int   readSFFcommonHeader(FILE *sff_fp);
char* readSFFrecord(FILE *sff_fp);

int main(int argc, char *argv[]){

    // List of adaptors
    char adaptorNames[MAX_ADAPTORS][MAX_LENGTH];
    char adaptorBases[MAX_ADAPTORS][MAX_LENGTH];
    bzero(adaptorNames,sizeof(adaptorNames));
    bzero(adaptorBases,sizeof(adaptorBases));

    // Read the adaptors
    size_t nAdaptors = readAdaptors(adaptorFileName,adaptorNames,adaptorBases);
    printf("nAdaptors = %ld\n",nAdaptors);

    // Find longest and shortest adaptor lengths
    size_t maxAdaptorLength = 0, minAdaptorLength = 1000000;
    for(size_t k=0; k<MAX_ADAPTORS; k++){
        size_t length = strlen( adaptorBases[k] );
        if( length == 0 ) continue;
        if( length > maxAdaptorLength ) maxAdaptorLength = length;
        if( length < minAdaptorLength ) minAdaptorLength = length;
    }

    // build a table
    LookUpTable hash2index(adaptorBases);
    printf("Found %ld collisions\n",hash2index.countCollisions());

//inputFileName = argv[1];

    // open input file:
    FILE *inputFile = NULL;
    if( (inputFile = fopen(inputFileName,"r")) == NULL ){
        printf("Cannot open %s\n",inputFileName);
        exit(0);
    }

    // make sure we remenber where the common header is in the input file
    long commonHeaderBegins = ftell(inputFile);
    int nReads = readSFFcommonHeader(inputFile);
    long commonHeaderLength = ftell(inputFile) - commonHeaderBegins;

    printf("iterating over %d reads\n",nReads);

    // every time we find an adaptor, remember the record position
    unsigned int adaptorsFound[MAX_ADAPTORS];
    long  *adaptorBegins[MAX_ADAPTORS]; //[nReads];
    long  *adaptorLength[MAX_ADAPTORS]; //[nReads];
    bzero(adaptorsFound, sizeof(adaptorsFound));
    for(size_t a=0; a<MAX_ADAPTORS; a++){
        adaptorBegins[a] = new long [nReads];
        adaptorLength[a] = new long [nReads];
        bzero(adaptorBegins[a], sizeof(long)*nReads);
        bzero(adaptorLength[a], sizeof(long)*nReads);
    }

    for(size_t read=0; read<nReads; read++){

        long recordBegins = ftell(inputFile);
        char *seq = readSFFrecord(inputFile);
        long recordLength = ftell(inputFile) - recordBegins;

        NumericSequence numSeq(seq);

        for(size_t i=0; i<strlen(seq); i++){
            unsigned long long view = numSeq.view(i,maxAdaptorLength);
            for(size_t len = maxAdaptorLength; len>=minAdaptorLength; len--){
                const unsigned long long mask[] = {0x0,
                             0x3,        0xF,        0x3F,        0xFF,        0x3FF,        0xFFF,        0x3FFF,        0xFFFF,
                         0x3FFFF,    0xFFFFF,    0x3FFFFF,    0xFFFFFF,    0x3FFFFFF,    0xFFFFFFF,    0x3FFFFFFF,    0xFFFFFFFF,
                     0x3FFFFFFFF,0xFFFFFFFFF,0x3FFFFFFFFF,0xFFFFFFFFFF,0x3FFFFFFFFFF,0xFFFFFFFFFFF,0x3FFFFFFFFFFF,0xFFFFFFFFFFFF}; 
                view &= mask[len];
                size_t ind = hash2index.find( view, len ); // terrible implementation of the hash!!!

                if( ind != MAX_ADAPTORS ){
                    adaptorBegins[ind][ adaptorsFound[ind] ] = recordBegins;
                    adaptorLength[ind][ adaptorsFound[ind] ] = recordLength;
                    adaptorsFound[ind]++;
                }
            }
        }

        free(seq);
    }


    // Now we need to write back sff files
    for(size_t a=0; a<MAX_ADAPTORS; a++){
        if( adaptorsFound[a] ){ //&& false ){
            // open output file:
            char name[MAX_LENGTH+4];
            sprintf(name,"%s.sff",adaptorNames[a]);
            FILE *outputFile = NULL;
            if( (outputFile = fopen(name,"w")) == NULL ){
                printf("Cannot open %s\n",adaptorNames[a]);
                exit(0);
            }
            // write the common header
            char  rawBuffer[ commonHeaderLength ];
            fseek(inputFile, commonHeaderBegins, SEEK_SET);
            fread(rawBuffer, commonHeaderLength, 1, inputFile);
            // a dirty hack to correct nreads
            uint32_t newNreads = htobe32( adaptorsFound[a] );
            fwrite(rawBuffer, 20, 1, outputFile);
            fwrite(&newNreads, sizeof(uint32_t), 1, outputFile);
            fwrite(rawBuffer+20+sizeof(uint32_t), commonHeaderLength - sizeof(uint32_t) - 20, 1, outputFile);

            for(size_t i=0; i<adaptorsFound[a]; i++){
                char rawBuffer[ adaptorLength[a][i] ];
                fseek (inputFile, adaptorBegins[a][i], SEEK_SET);
                fread (rawBuffer, adaptorLength[a][i], 1, inputFile);
                fwrite(rawBuffer, adaptorLength[a][i], 1, outputFile);
            }
            fclose(outputFile);
        }
        printf("%s: %d\n",adaptorNames[a],adaptorsFound[a]);
    }

    fclose(inputFile);
    return 0;
}

#undef MAX_ADAPTORS
#undef MAX_LENGTH



///////// Inhereted code
#include "sff2fastq/sff.h"
#define VERSION "0.9.0"
#define PRG_NAME "sff2fastq"
int trim_flag  = 1;  // trimming is enabled by default -- like 454 tools

sff_common_header h;
sff_read_header rh;
sff_read_data rd;

int readSFFcommonHeader(FILE *sff_fp){
    read_sff_common_header(sff_fp, &h);
    verify_sff_common_header(PRG_NAME, VERSION, &h);
    return int(h.nreads);
}

// free_sff_common_header(&h);

// read records one-by-one (it doesn't look like there is a random access functionality)
char* readSFFrecord(FILE *sff_fp){
    int left_clip = 0, right_clip = 0, nbases = 0;
    char *name;
    char *bases;
    uint8_t *quality;
        read_sff_read_header(sff_fp, &rh);
        read_sff_read_data(sff_fp, &rd, h.flow_len, rh.nbases);

        // get clipping points
        get_clip_values(rh, trim_flag, &left_clip, &right_clip);
        nbases = right_clip - left_clip;

        // create bases string
        bases = get_read_bases(rd, left_clip, right_clip);

        // create quality array
        quality = get_read_quality_values(rd, left_clip, right_clip);

        // create read name string
        int name_length = (int) rh.name_len + 1; // account for NULL termination
        name = (char *) malloc( name_length * sizeof(char) );
        if (!name) {
            fprintf(stderr, "Out of memory! For read name string!\n");
            exit(1);
        }
        memset(name, '\0', (size_t) name_length);
        strncpy(name, rh.name, (size_t) rh.name_len);
        
        free(name);
///        free(bases);
        free(quality);
        free_sff_read_header(&rh);
        free_sff_read_data(&rd);

return bases;
}
*/
