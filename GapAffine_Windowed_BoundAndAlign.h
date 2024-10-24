#ifndef GAPAFFINE_WINDOWED_BOUNDANDALIGN_H
#define GAPAFFINE_WINDOWED_BOUNDANDALIGN_H

#include "main.h"

// Function prototypes
void update_window_bound(char *Q, char *T, GapAffine_Parameters *ga_params, 
                         int start_i, int start_j, int end_i, int end_j, int is_last_window);
void backtrace_bucle(char *Q, char *T, GapAffine_Parameters *ga_params,
                     int end_i, int end_j, int *i, int *j, char *cigar, int *cigar_len);
int calculate_cigar_score(GapAffine_Parameters *ga_params, char *cigar, int cigar_len);
char* print_cigar_windowed(char *cigar_ops, int cigar_len, int score);
int windowed_gapAffine_bound(char *Q, char *T, GapAffine_Parameters *ga_params, char *cigar, int *cigar_len, int *len_query, int *len_target);
int windowed_gapAffine_align(int bound, char *Q, char *T, GapAffine_Parameters *ga_params,
                    char *cigar, int *cigar_len, int *len_query, int *len_target);
void GapAffine_windowed(char *Q, char *T, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res);


#endif // GAPAFFINE_WINDOWED_BOUNDANDALIGN_H
