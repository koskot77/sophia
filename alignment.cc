#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#define MIN(A,B,C) ( A<B ? ( B<C ? A : ( A<C ? A : C ) ) : ( B<C ? B : C ) )
// penalties
const size_t gapCost =  5; // gap
const size_t misCost = 10; // mismatch

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
    char *tmp = new char [len1+len2];
    memcpy(tmp,x,sizeof(char)*(len1+len2));
    for(size_t i=0; i<k; i++)
        x[k-1-i] = tmp[i];
    memcpy(tmp,y,sizeof(char)*(len1+len2));
    for(size_t i=0; i<k; i++)
        y[k-1-i] = tmp[i];
    delete [] tmp;
}
/*
int main(void){

const char *seq1 = "CGTGTCTCTAAGACTCGGCAGCATCTCCATCATGTGGTTTATGCAGCAGATGCAAGGTATTCTGTAAAGGTTCTTGGTATACCTGTTTCGTAACAACATGAGTAGTCTCTTCAGTAATTAGATTAGTTAAAGTGATGTGGTGTTTCTGGCAAACTTGTACACGAGCATCTGAAATTAAATCAAATATTCCATTATCATGAGTTACCTCTAGCACACAGCTCAGAATACTAGTTATCCCACCATGGCATATGTTTACCTACGTAGCAATCCTGCACGTTCTACACGTGTCCTGGAACTATTTAAAGTGATGGAGAACAGTGACGATCGCTAGAGACACG";
const char *seq2 = "CGTGTCTCTAAGACTCGGCAGCATCTCCATCCCAAATGGTCTTCAGAATAATCTAATTACAGTACTGTATCTACCCACTCTCTTTTCAGTGCCTGTTAAGTTGGCAAACTTTGCCATTACCCTTTTTTGCAGAATCCAAACTGATTTCATCCCTGGTTCCTTGAGGGGTGATTTGTAACAATTCTTGATCTCTCACACTATAGGGAAAAGACAGAGTCCTAATAAGAAACACTAGTTACATGTATGCAGAACTGTCAAATGACCAAGATCAAACATTTTAGCTCTTTCGATTACAGAAAGCTGACCAATCTTATTTAGTTAGCGAAAGCTGCTCTCTCCTGGAGAACAGTGACGATCGCTAGAGACACG";

    size_t len1 = strlen(seq1);
    size_t len2 = strlen(seq2);

    char a[ len1+len2 ], b[ len1+len2 ];
    bzero(a,sizeof(a));
    bzero(b,sizeof(b));

    size_t score[len1+1][len2+1];
    bzero(score,sizeof(score));

    size_t s = alignmentScoreMatrix(seq1, seq2, (size_t**)score);

    // reconstruction
    printf("Score: %ld, normalized score: %f\n",s,s/double(len1+len2));

    reconstruction(seq1, seq2, (const size_t **)score, a, b);

////////////////////
    for(size_t i=0; i<strlen(a); i++)
        printf("%c",a[i]);
    printf("\n");
    for(size_t i=0; i<strlen(b); i++)
        printf("%c",b[i]);
    printf("\n");
////////////////////

    return 0;
}
*/
