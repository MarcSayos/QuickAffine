// main.h
#ifndef MAIN_H
#define MAIN_H

#define NEG_INF INT_MIN
#define MAX_LEN 100000  // Define a maximum sequence length, adjust as needed
#define MAX2(a,b) (a > b ? a : b)
#define MIN2(a,b) (a < b ? a : b)
#define MIN3(a,b,c) (a < b ? MIN2(a,c) : MIN2(b,c))
#define DEBUG 4 // 0 for real tests, 1 for terminal tests, 3 for python plots, 4 for postprocessing
#define ELAPSED_DIV 1000 // 1000 to print seconds, 1 to print milliseconds
#define CELLS_DIV 1000000 // 1000000 to print millions, 1 to print normal
#define PATH_MAX 300

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <time.h>
#include <math.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/resource.h>

typedef struct GapAffine_Parameters{
    int ws;                         // Window Size
    int os;                         // Overlap Size
    char *penalty_set;              // Name (if any) of the penalty set
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
    int swg;                        // Boolean, should it compute SWG
    int windowed;                   // Boolean, should it compute Windowed
    int banded;                     // Boolean, should it compute Banded (only if Windowed was computed previously)
    int parasail;                   // Boolean, should it compute Parasail
    int upper_bound;                // Upper-bound (only if computing Banded and Windowed was computed previously)
    char *dataset_type;             // Real or Simulated
} GapAffine_Parameters;

typedef struct GapAffine_Alignment{
    char *query;                    // Query sequence
    char *target;                   // Target sequence
    int len_query;                  // Query length
    int len_target;                 // Target length
    char *cigar;                    // CIGAR alignment
    int cigar_len;                  // CIGAR length
    char *cigar_position;           // Pointer to the position in the CIGAR allocation
    uint16_t **M;                   // Matches matrix
    uint16_t **I;                   // Insertions matrix
    uint16_t **D;                   // Deletions matrix
} GapAffine_Alignment;

typedef struct GapAffine_Results{
    int score;                      // Final alignment score
    int original_score;             // Computed score adjusted with original penalties
    double memory;                     // Memory utilized in the execution
    double elapsed;                 // Elapsed time during execution
    int cells;                      // Computed cells
} GapAffine_Results;

typedef struct GapAffine_Totals{
    double elapsed_windowed;        // Elapsed time windowed algorithm
    double elapsed_banded;          // Elapsed time banded algorithm
    double elapsed_SWG;             // Elapsed time SWG algorithm
    double elapsed_parasail_scan;   // Elapsed time parasail algorithm scan
    double elapsed_parasail_diag;   // Elapsed time parasail algorithm diag
    int cells_windowed;             // Computed cells windowed algorithm
    int cells_banded;               // Computed cells banded algorithm
    int cells_SWG;                  // Computed cells SWG algorithm
    double memory_windowed;            // Allocated memory windowed algorithm
    double memory_banded;              // Allocated memory banded algorithm
    double memory_SWG;                 // Allocated memory SWG algorithm
    double memory_parasail_scan;       // Allocated memory parasail algorithm scan
    double memory_parasail_diag;       // Allocated memory parasail algorithm diag
    int score_windowed;             // Sum of scores (check for windowed accuracy)
    int score_banded;               // Sum of scores (check for windowed accuracy)
    int score_SWG;                  // Sum of scores (check for windowed accuracy)
    int score_parasail_scan;        // Sum of scores parasail scan
    int score_parasail_diag;        // Sum of scores parasail diag
    int avg_query_length;
    int num_queries;
} GapAffine_Totals;

// Functions
double get_memory_usage();
uint16_t **create_matrix(int rows, int cols);

void costs_transform(GapAffine_Parameters *ga_params);
void read_inputs(int argc, char *argv[], GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Totals *ga_totals, 
                char *name, char *input_file_path, char *output_file_path, char *output_scores_file_path);
void print_postprocessing(FILE *outfile, GapAffine_Parameters *ga_params, GapAffine_Totals *ga_totals, char *name);

#endif