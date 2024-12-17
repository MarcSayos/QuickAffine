#ifndef GAPAFFINE_WINDOWED_H
#define GAPAFFINE_WINDOWED_H

#include "main.h"

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
    int is_last_window;             // Bool checking whether it's the last window being computed
} Windowed_positions;


// Function prototypes
void windowed_local_backtrace(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res, Windowed_positions *win_pos);
void windowed_compute_window(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res, Windowed_positions *win_pos);
void windowed_update_positions(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res, Windowed_positions *win_pos);
void windowed_GapAffine(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res);


#endif // GAPAFFINE_WINDOWED_H