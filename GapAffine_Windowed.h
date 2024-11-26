#ifndef GAPAFFINE_WINDOWED_H
#define GAPAFFINE_WINDOWED_H

#include "main.h"

// Function prototypes
void windowed_local_backtrace(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res, Windowed_positions *win_pos);
void windowed_compute_window(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res, Windowed_positions win_pos);
void windowed_update_positions(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res, Windowed_positions *win_pos);
void windowed_GapAffine(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res);


#endif // GAPAFFINE_WINDOWED_H