#ifndef GAPAFFINE_BANDED_H
#define GAPAFFINE_BANDED_H

#include "main.h"

// void banded_backtrace(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, int end_i, int end_j, int *i, int *j);
void banded_backtrace_v2(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res);
void banded_calculate_cigar_score(GapAffine_Parameters *ga_params, GapAffine_Alignment *ga_algn, GapAffine_Results *ga_res);
char* banded_print_cigar(GapAffine_Alignment *ga_algn);
void banded_reverse_cigar(GapAffine_Alignment *ga_algn);
void banded_GapAffine(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res);


#endif // GAPAFFINE_BANDED_H