// main.h
#ifndef MAIN_H
#define MAIN_H

#define NEG_INF INT_MIN
#define MAX_LEN 20000  // Define a maximum sequence length, adjust as needed
#define MIN2(a,b) (a < b ? a : b)
#define MIN3(a,b,c) (a < b ? MIN2(a,c) : MIN2(b,c))

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <time.h>
#include <math.h>
#include <sys/resource.h>

typedef struct {
    const int bases;
    const int ws;   // Window Size
    const int os;   // Overlap Size
    int Cm;         // Cost Match
    int Cx;         // Cost Mismatch
    int Co;         // Cost Open
    int Ci;         // Cost Extend Insertion
    int Cd;         // Cost Extend Deletion
    int alpha;
    int beta;
    int gamma;
} GapAffine_Parameters;

typedef struct {
    char *query;
    char *target;
    int len_query;
    int len_target;
    char *cigar;
    int cigar_len;
    int **M;
    int **I;
    int **D;
} GapAffine_Alignment;

typedef struct {
    int score;
    int bound;
    int memory;
    double elapsed;
} GapAffine_Results;



// Functions
int get_memory_usage();
__uint32_t **create_matrix(int rows, int cols);


#endif