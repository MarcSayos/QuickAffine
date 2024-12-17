#include "GapAffine_Banded.h"

void banded_backtrace(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res) {
    int i = ga_algn->len_query;
    int j = ga_algn->len_target;

    int current_matrix =    (ga_algn->M[i][j] <= ga_algn->I[i][j] && ga_algn->M[i][j] <= ga_algn->D[i][j]) ? 'M' : 
                            (ga_algn->I[i][j] <= ga_algn->M[i][j] && ga_algn->I[i][j] <= ga_algn->D[i][j]) ? 'I' : 'D';
    while (i > 0 || j > 0) {
        if (i == 0) {
            if (current_matrix != 'I') {
                ga_res->score += ga_params->Co; 
                ga_res->original_score += ga_params->or_Co;
                current_matrix = 'I';
            }
            ga_res->score +=  ga_params->Ci; 
            ga_res->original_score +=  ga_params->or_Ci;
            (j)--;
        } else if (j == 0) {
            if (current_matrix != 'D') {
                ga_res->score += ga_params->Co; 
                ga_res->original_score += ga_params->or_Co;
                current_matrix = 'D';
            }
            ga_res->score +=  ga_params->Cd; 
            ga_res->original_score +=  ga_params->or_Cd;
            (i)--;
        } else if (current_matrix == 'M') {
            int C = (ga_algn->query[i-1] == ga_algn->target[j-1]) ? ga_params->Cm : ga_params->Cx;
            int or_C = (ga_algn->query[i-1] == ga_algn->target[j-1]) ? ga_params->or_Cm : ga_params->or_Cx;
            if      (ga_algn->M[i][j] == ga_algn->M[i-1][j-1] + C) current_matrix = 'M'; 
            else if (ga_algn->M[i][j] == ga_algn->I[i-1][j-1] + C) current_matrix = 'I';
            else if (ga_algn->M[i][j] == ga_algn->D[i-1][j-1] + C) current_matrix = 'D';
            ga_res->score += C; 
            ga_res->original_score += or_C;
            i--;
            j--;
        } else if (current_matrix == 'I') {
            if (ga_algn->I[i][j] == ga_algn->M[i-1][j] + ga_params->Co + ga_params->Ci) {
                current_matrix = 'M';
                ga_res->score += ga_params->Co;
                ga_res->original_score += ga_params->or_Co;
            } else if (ga_algn->I[i][j] == ga_algn->I[i-1][j] + ga_params->Ci) {}
            ga_res->score += ga_params->Ci;
            ga_res->original_score += ga_params->or_Ci;
            i--;
        } else if (current_matrix == 'D') {
            if (ga_algn->D[i][j] == ga_algn->M[i][j-1] + ga_params->Co + ga_params->Cd) {
                current_matrix = 'M';
                ga_res->score += ga_params->Co;
                ga_res->original_score += ga_params->or_Co;
            } else if (ga_algn->D[i][j] == ga_algn->D[i][j-1] + ga_params->Cd) {}
            ga_res->score += ga_params->Cd;
            ga_res->original_score += ga_params->or_Cd;
            j--;
        }
    }
}


// Function to calculate the CIGAR score
void banded_calculate_cigar_score(GapAffine_Parameters *ga_params, GapAffine_Alignment *ga_algn, GapAffine_Results *ga_res) {

    ga_res->score = 0, ga_res->original_score = 0;
    for (int i = 0; i < ga_algn->cigar_len; i++) {
        switch (ga_algn->cigar[i]) {
            case 'M':
                ga_res->score += ga_params->Cm;
                ga_res->original_score += ga_params->or_Cm;
                break;
            case 'X':
                ga_res->score += ga_params->Cx;
                ga_res->original_score += ga_params->or_Cx;
                break;
            case 'I':
                ga_res->score += (i == 0 || ga_algn->cigar[i-1] != 'I') ? (ga_params->Ci + ga_params->Co) : ga_params->Ci;
                ga_res->original_score += (i == 0 || ga_algn->cigar[i-1] != 'I') ? (ga_params->or_Ci + ga_params->or_Co) : ga_params->or_Ci;
                break;
            case 'D':
                ga_res->score += (i == 0 || ga_algn->cigar[i-1] != 'D') ? (ga_params->Cd + ga_params->Co) : ga_params->Cd;
                ga_res->original_score += (i == 0 || ga_algn->cigar[i-1] != 'D') ? (ga_params->or_Cd + ga_params->or_Co) : ga_params->or_Cd;
                break;
        }
    }
}

void banded_initialize_matrices(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res) {
    // Create matrices
    ga_algn->M = create_matrix(ga_algn->len_query + 1, ga_algn->len_target + 1);
    ga_algn->I = create_matrix(ga_algn->len_query + 1, ga_algn->len_target + 1);
    ga_algn->D = create_matrix(ga_algn->len_query + 1, ga_algn->len_target + 1);

    // Initialize values
    ga_algn->M[0][0] = 0;  
    ga_algn->I[0][0] = 0;
    ga_algn->D[0][0] = 0;
    for (int i = 1; i <= ga_algn->len_query; i++) {
        if (ga_params->Co + ga_params->Cd * i > ga_params->upper_bound) break;
        ga_res->cells += 1;
        ga_algn->M[i][0] = ga_params->Co + ga_params->Cd * i;
        ga_algn->I[i][0] = __UINT16_MAX__;
        ga_algn->D[i][0] = ga_params->Co + ga_params->Cd * i;
    }
    for (int j = 1; j <= ga_algn->len_target; j++) {
        if (ga_params->Co + ga_params->Ci * j > ga_params->upper_bound) break;
        ga_res->cells += 1;
        ga_algn->M[0][j] = ga_params->Co + ga_params->Ci * j;
        ga_algn->I[0][j] = ga_params->Co + ga_params->Ci * j;
        ga_algn->D[0][j] = __UINT16_MAX__;
    }
}

void banded_GapAffine(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res) {
    

    ga_algn->cigar = (char *)malloc((ga_algn->len_query + ga_algn->len_target) * sizeof(char) * 2);
    if (ga_algn->cigar == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
    }
    ga_algn->cigar_len = 0;
    ga_res->cells = 0;

    banded_initialize_matrices(ga_algn, ga_params, ga_res);
    
    int start_time = clock();
    ga_res->memory = get_memory_usage();
    
    int bottom_j = 1;
    // Align Score Calculation
    for (int i = 1; i <= ga_algn->len_query; i++) {
        for (int j = bottom_j; j <= ga_algn->len_target; j++) {
            int Cost = (ga_algn->query[i - 1] == ga_algn->target[j - 1]) ? ga_params->Cm : ga_params->Cx;

            ga_algn->M[i][j] = MIN2(MIN3(ga_algn->M[i - 1][j - 1], ga_algn->I[i - 1][j - 1], ga_algn->D[i - 1][j - 1]) + Cost, __UINT16_MAX__);
            ga_algn->I[i][j] = MIN3(ga_algn->I[i-1][j] + ga_params->Ci, ga_algn->M[i-1][j] + ga_params->Co + ga_params->Ci, __UINT16_MAX__);
            ga_algn->D[i][j] = MIN3(ga_algn->D[i][j-1] + ga_params->Cd, ga_algn->M[i][j-1] + ga_params->Co + ga_params->Cd, __UINT16_MAX__);

            ga_res->cells += 1;

            // Threshold check
            if (ga_algn->len_query > ga_params->upper_bound && ga_algn->M[i][j] > ga_params->upper_bound && ga_algn->I[i][j] > ga_params->upper_bound && ga_algn->D[i][j] > ga_params->upper_bound) {
                if ((i-1 == 0) || ga_algn->M[i-1][j] == __UINT16_MAX__) { // Top
                    ga_algn->M[i][j] = ga_algn->I[i][j] = ga_algn->D[i][j] = __UINT16_MAX__;
                    break;
                } else if ((j-1 == 0) || ga_algn->M[i][j-1] == __UINT16_MAX__) {
                    bottom_j++;
                    ga_algn->M[i][j] = ga_algn->I[i][j] = ga_algn->D[i][j] = __UINT16_MAX__;
                }
            }
        }
    }
    ga_res->score = MIN3(ga_algn->M[ga_algn->len_query][ga_algn->len_target],
                         ga_algn->I[ga_algn->len_query][ga_algn->len_target],
                         ga_algn->D[ga_algn->len_query][ga_algn->len_target]);
    // banded_backtrace(ga_algn, ga_params, ga_res);
    
    int end_time = clock();
    ga_res->elapsed = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000.0;
    ga_res->memory = get_memory_usage() - ga_res->memory;
}