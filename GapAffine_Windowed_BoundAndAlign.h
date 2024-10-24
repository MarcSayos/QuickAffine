#ifndef GAPAFFINE_WINDOWED_BOUNDANDALIGN_H
#define GAPAFFINE_WINDOWED_BOUNDANDALIGN_H

#include "main.h"

// Function prototypes
void update_window_bound(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, int start_i, int start_j, int end_i, int end_j, int is_last_window);
void backtrace_bucle(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, int end_i, int end_j, int *i, int *j);
void calculate_cigar_score(GapAffine_Parameters *ga_params, GapAffine_Alignment *ga_algn, GapAffine_Results *ga_res);
void calculate_cigar_bound(GapAffine_Parameters *ga_params, GapAffine_Alignment *ga_algn, GapAffine_Results *ga_res);
char* print_cigar_windowed(GapAffine_Alignment *ga_algn, int score);
void windowed_gapAffine_bound(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res);
void windowed_gapAffine_align(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res);
void GapAffine_windowed(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res);


#endif // GAPAFFINE_WINDOWED_BOUNDANDALIGN_H
