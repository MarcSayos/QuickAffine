// main.h
#ifndef MAIN_H
#define MAIN_H

#define NEG_INF INT_MIN
#define MAX_LEN 100000  // Define a maximum sequence length, adjust as needed
#define MAX2(a,b) (a > b ? a : b)
#define MIN2(a,b) (a < b ? a : b)
#define MIN3(a,b,c) (a < b ? MIN2(a,c) : MIN2(b,c))
#define DEBUG 1 // 0 for real tests, 1 for terminal tests, 3 for python plots, 4 for postprocessing
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
#include "henly.h"

typedef struct GapAffine_Parameters{
    int ws;                         // Window Size
    int os;                         // Overlap Size
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
    int upper_bound;                // Upper-bound (only if computing Banded and Windowed was computed previously)
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
    int memory;                     // Memory utilized in the execution
    double elapsed;                 // Elapsed time during execution
    int cells;                      // Computed cells
} GapAffine_Results;

typedef struct Position{
    int i;                          // i
    int j;                          // j
} Position;

typedef struct {
    Position qt_bottom;             // Bottom right corner of query/target
    Position qt_top;                // Top left corner of query/target
    Position w_bottom;              // Bottom right corner of matrices
    Position last_stop;             // The stoping point of the last window
    int current_matrix;             // 0 for M, 1 for I, 2 for D
} Windowed_positions;

typedef struct GapAffine_Totals{
    double elapsed_windowed;        // Elapsed time windowed algorithm
    double elapsed_banded;          // Elapsed time banded algorithm
    double elapsed_SWG;             // Elapsed time SWG algorithm
    int cells_windowed;             // Computed cells windowed algorithm
    int cells_banded;               // Computed cells banded algorithm
    int cells_SWG;                  // Computed cells SWG algorithm
    int memory_windowed;            // Allocated memory windowed algorithm
    int memory_banded;              // Allocated memory banded algorithm
    int memory_SWG;                 // Allocated memory SWG algorithm
} GapAffine_Totals;

// Functions
int get_memory_usage();
uint16_t **create_matrix(int rows, int cols);
void reset_matrices(GapAffine_Alignment *ga_algn, int previous_length, int cols, int rows);
void python_plot_print(GapAffine_Alignment *ga_algn);

// void read_inputs(int argc, char *argv[], GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Totals *ga_totals, char *name, char *input_file_path, char *output_file_path);
void print_qt_results(FILE *outfile, GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res_swg, GapAffine_Results *ga_res_windowed, GapAffine_Results *ga_res_banded);
void print_total_results(FILE *outfile, GapAffine_Parameters *ga_params, GapAffine_Totals *ga_totals, char *name);

void QuickAffine(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res_windowed, GapAffine_Results *ga_res_banded);

#endif