#ifndef GAPAFFINE_WINDOWED_BOUNDANDALIGN_H
#define GAPAFFINE_WINDOWED_BOUNDANDALIGN_H

#include "main.h"

// Function prototypes
Position BOUND_local_window_backtrace(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res, Windowed_positions win_pos, int is_last_window);
void BOUND_compute_window(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res, Windowed_positions win_pos, int is_last_window);
void BOUND_windowed_gapAffine(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res);

void ALIGN_backtrace_bucle(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, int end_i, int end_j, int *i, int *j);
void ALIGN_calculate_cigar_score(GapAffine_Parameters *ga_params, GapAffine_Alignment *ga_algn, GapAffine_Results *ga_res);
char* ALIGN_print_cigar(GapAffine_Alignment *ga_algn);
void ALIGN_reverse_cigar(GapAffine_Alignment *ga_algn);
void ALIGN_windowed_gapAffine(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res);

void GapAffine_windowed(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res);


#endif // GAPAFFINE_WINDOWED_BOUNDANDALIGN_H
