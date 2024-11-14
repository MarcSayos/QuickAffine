// main.h
#ifndef MAIN_H
#define MAIN_H

#define NEG_INF INT_MIN
#define MAX_LEN 20000  // Define a maximum sequence length, adjust as needed
#define MAX2(a,b) (a > b ? a : b)
#define MIN2(a,b) (a < b ? a : b)
#define MIN3(a,b,c) (a < b ? MIN2(a,c) : MIN2(b,c))
#define DEBUG 1 // 0 for real tests, 1 for terminal tests, 3 for python plots

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <time.h>
#include <math.h>
#include <sys/resource.h>

typedef struct GapAffine_Parameters{
    const int bases;
    const int ws;                   // Window Size
    const int os;                   // Overlap Size
    int Cm;                         // Match Cost
    int Cx;                         // Mismatch Cost
    int Co;                         // Open Cost
    int Ci;                         // Extend Insertion Cost
    int Cd;                         // Extend Deletion Cost
    int or_Cm;                      // Original Match Cost
    int or_Cx;                      // Original Mismatch Cost
    int or_Co;                      // Original Open Cost
    int or_Ci;                      // Original Extend Insertion Cost
    int or_Cd;                      // Original Extend Deletion Cost
    int alpha;                      // Cost transformation values
    int beta;                       // Cost transformation values
    int gamma;                      // Cost transformation values
} GapAffine_Parameters;

typedef struct GapAffine_Alignment{
    char *query;                    // Query sequence
    char *target;                   // Target sequence
    int len_query;                  // Query length
    int len_target;                 // Target length
    char *cigar;                    // CIGAR alignment
    int cigar_len;                  // CIGAR length
    uint16_t **M;                   // Matches matrix
    uint16_t **I;                   // Insertions matrix
    uint16_t **D;                   // Deletions matrix
} GapAffine_Alignment;

typedef struct GapAffine_Results{
    int score;                      // Final alignment score
    int original_score;             // Computed score adjusted with original penalties
    int bound;                      // Upper-bound from windowed
    int bound_size;                 // Expected upper-bound size in the matrix
    int memory;                     // Memory utilized in the execution
    double elapsed;                 // Elapsed time during execution
    int computed_cells_score;       // Computed cells during normal GapAffine
    int computed_cells_windowed;    // Computed cells during windowed upper-bound calculation
    int computed_cells_banded;      // Computed cells during exact banded alignment
} GapAffine_Results;



// Functions
int get_memory_usage();
uint16_t **create_matrix(int rows, int cols);
void reset_matrices(GapAffine_Alignment *ga_algn);
void python_plot_print(GapAffine_Alignment *ga_algn);


#endif