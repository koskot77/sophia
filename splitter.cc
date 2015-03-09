#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "./toolbox.h"

// Input file names
const char *inputFileName   = "input.sff";
const char *adaptorFileName = "ionXpress_barcode.txt";

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
    sprintf(format,"%%%ds\t%%%ds\n",MAX_LENGTH-1,MAX_LENGTH-1); 

    unsigned int num = 0;
    while( !feof(file) && num<MAX_ADAPTORS && 
            fscanf(file, format, names[num], bases[num])==2 ) num++;

    fclose(file);

    if( num == MAX_ADAPTORS )
        printf("Warning: only first %d adaptors are to be processed\n",num);

    return num;
}


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
    if( maxAdaptorLength - minAdaptorLength > 4 ){
        printf("Adaptors should have length variation within 4 symbols max\n");
        return 0;
    }

    // build a table
    LookUpTable number2index(adaptorBases);
    printf("Found %ld collisions\n",hash2index.countCollisions());

    // open input file:
    FILE *inputFile = NULL;
    if( (inputFile = fopen(inputFileName,"r")) == NULL ){
        printf("Cannot open %s\n",inputFileName);
        exit(0);
    }

    // make sure we remember where the common header is in the input file
    long commonHeaderBegins = ftell(inputFile);
    int nReads = readSFFcommonHeader(inputFile);
    long commonHeaderLength = ftell(inputFile) - commonHeaderBegins;

    printf("iterating over %d reads\n",nReads);

//#include "l.h"

    // every time we find an adaptor, we will remember the record position
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
        if( numSeq.error() ){
            printf("Warning\n");
        }

        if( strlen(seq) < minAdaptorLength ) continue;

        for(size_t i=0; i<strlen(seq) - minAdaptorLength; i++){
            unsigned long long view = numSeq.view(i,maxAdaptorLength);
            size_t len = maxAdaptorLength;
            size_t ind = number2index.find( view, len );
            if( ind != MAX_ADAPTORS ){
                adaptorBegins[ind][ adaptorsFound[ind] ] = recordBegins;
                adaptorLength[ind][ adaptorsFound[ind] ] = recordLength;
                adaptorsFound[ind]++;
            }
        }
        free(seq);
    }

/*
    printf("unsigned int adaptorsFound[MAX_ADAPTORS] = {\n");
    for(size_t a=0; a<MAX_ADAPTORS; a++){
        printf("%d,",adaptorsFound[a]);
        if((a%19)==0) printf("\n");
    }
    printf("}\n");

    printf("long adaptorBegins[MAX_ADAPTORS][1209] = {\n");
    for(size_t a=0; a<MAX_ADAPTORS; a++){
        printf("{");
        for(size_t i=0; i<1209; i++){
            printf("%d",adaptorBegins[a][i]);
            if(i!=1208) printf(",");
            if(((i+1)%80)==0) printf("\n");
        }
        printf("}\n");
        if(a!=MAX_ADAPTORS-1) printf(",\n");
    }
    printf("}\n");

    printf("long adaptorLength[MAX_ADAPTORS][1209] = {\n");
    for(size_t a=0; a<MAX_ADAPTORS; a++){
        printf("{");
        for(size_t i=0; i<1209; i++){
            printf("%d",adaptorLength[a][i]);
            if(i!=1208) printf(",");
            if(((i+1)%80)==0) printf("\n");
        }
        printf("}\n");
        if(a!=MAX_ADAPTORS-1) printf(",\n");
    }
    printf("}\n");

return 0;
*/
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
        if( strlen(adaptorNames[a]) )
            printf("%s: %d\n",adaptorNames[a],adaptorsFound[a]);
    }

    fclose(inputFile);
    return 0;
}

#undef MAX_ADAPTORS
#undef MAX_LENGTH



///////// Inherited code
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

