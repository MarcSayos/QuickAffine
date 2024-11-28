#include "GapAffine_Windowed.h"

// Function to compute the window bounds
void windowed_compute_window(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res, Windowed_positions win_pos) {
    ga_algn->M[0][0] = 0;
    ga_algn->I[0][0] = 0;
    ga_algn->D[0][0] = 0;

    int is_last_window = 0;
    if (win_pos.qt_top.i == 0 && win_pos.qt_top.j == 0) is_last_window = 1;

    for(int i = 1; i <= win_pos.w_bottom.i; i++){
        ga_algn->M[i][0] = (is_last_window ? ga_params->Co + ga_params->Cd * i : 0);
        ga_algn->I[i][0] = UINT16_MAX;
        ga_algn->D[i][0] = (is_last_window ? ga_params->Co + ga_params->Cd * i : 0);
    }
    for (int j = 1; j <= win_pos.w_bottom.j; j++) {
        ga_algn->M[0][j] = (is_last_window ? ga_params->Co + ga_params->Ci * j : 0);
        ga_algn->I[0][j] = (is_last_window ? ga_params->Co + ga_params->Ci * j : 0);
        ga_algn->D[0][j] = UINT16_MAX;
    }

    for (int i = 1; i <= win_pos.w_bottom.i; i++) {
        for (int j = 1; j <= win_pos.w_bottom.j; j++) {
            int C = (ga_algn->query[win_pos.qt_top.i + i - 1] == ga_algn->target[win_pos.qt_top.j + j - 1]) ? ga_params->Cm : ga_params->Cx;
            ga_algn->M[i][j] = MIN3(ga_algn->M[i-1][j-1], ga_algn->I[i-1][j-1], ga_algn->D[i-1][j-1]) + C;
            ga_algn->I[i][j] = MIN2(ga_algn->I[i-1][j], ga_algn->M[i-1][j] + ga_params->Co) + ga_params->Ci;
            ga_algn->D[i][j] = MIN2(ga_algn->D[i][j-1], ga_algn->M[i][j-1] + ga_params->Co) + ga_params->Cd;
            ga_res->cells += 1;
        }
    }
}

void windowed_local_backtrace(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res, Windowed_positions *win_pos) {
    int i = win_pos->w_bottom.i;
    int j = win_pos->w_bottom.j;
    int os = ga_params->os;
    int stop_i_at_overlap = (win_pos->qt_top.i == 0 ? 0 : 1);
    int stop_j_at_overlap = (win_pos->qt_top.j == 0 ? 0 : 1);
    if (win_pos->current_matrix == -1) {
        win_pos->current_matrix =   (ga_algn->M[i][j] <= ga_algn->I[i][j] && ga_algn->M[i][j] <= ga_algn->D[i][j]) ? 0 : 
                                    (ga_algn->I[i][j] <= ga_algn->M[i][j] && ga_algn->I[i][j] <= ga_algn->D[i][j]) ? 1 : 2;
    }
    int start_matrix = win_pos->current_matrix;

    while (
        (i > os && j > os) ||
        (i > os && j <= os && !stop_j_at_overlap) ||
        (j > os && i <= os && !stop_i_at_overlap) ||
        ((i <= os && !stop_i_at_overlap) && (j <= os && !stop_j_at_overlap))
    ) {
        if (i == 0 || j == 0) {
            i = j = 0;
            break;
        }
        
        if (win_pos->current_matrix == 0) {
            int C = (ga_algn->query[win_pos->qt_top.i + i - 1] == ga_algn->target[win_pos->qt_top.j + j - 1]) ? ga_params->Cm : ga_params->Cx;
            if (ga_algn->M[i][j] == ga_algn->I[i-1][j-1] + C) win_pos->current_matrix = 1;
            else if (ga_algn->M[i][j] == ga_algn->D[i-1][j-1] + C) win_pos->current_matrix = 2;
            i--;
            j--;
        } else if (win_pos->current_matrix == 1) {
            if (ga_algn->I[i][j] == ga_algn->M[i-1][j] + ga_params->Co + ga_params->Ci) win_pos->current_matrix = 0;
            i--;
        } else if (win_pos->current_matrix == 2) {
            if (ga_algn->D[i][j] == ga_algn->M[i][j-1] + ga_params->Co + ga_params->Cd) win_pos->current_matrix = 0;
            j--;
        }
    }
    
    if (start_matrix == 0) ga_res->score += ga_algn->M[win_pos->w_bottom.i][win_pos->w_bottom.j];
    else if (start_matrix == 1) ga_res->score += ga_algn->I[win_pos->w_bottom.i][win_pos->w_bottom.j];
    else ga_res->score += ga_algn->D[win_pos->w_bottom.i][win_pos->w_bottom.j];
    win_pos->last_stop.i = i;
    win_pos->last_stop.j = j;
}

void windowed_update_positions(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res, Windowed_positions *win_pos) {
    
    int shift_i = win_pos->w_bottom.i - win_pos->last_stop.i;
    int shift_j = win_pos->w_bottom.j - win_pos->last_stop.j;

    // Update Windows
    win_pos->qt_bottom.i -= shift_i;
    win_pos->qt_bottom.j -= shift_j;
    win_pos->qt_top.i = MAX2(0, win_pos->qt_top.i - shift_i);
    win_pos->qt_top.j = MAX2(0, win_pos->qt_top.j - shift_j);
    win_pos->w_bottom.i = MIN2(ga_params->ws, win_pos->qt_bottom.i);
    win_pos->w_bottom.j = MIN2(ga_params->ws, win_pos->qt_bottom.j);
}


void windowed_GapAffine(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res) {
    int start_time = clock();
    ga_res->memory = get_memory_usage();
    
    // Initialize matrices
    ga_algn->M = create_matrix(MIN2(ga_params->ws, ga_algn->len_query) + 1, MIN2(ga_params->ws, ga_algn->len_target) + 1);
    ga_algn->I = create_matrix(MIN2(ga_params->ws, ga_algn->len_query) + 1, MIN2(ga_params->ws, ga_algn->len_target) + 1);
    ga_algn->D = create_matrix(MIN2(ga_params->ws, ga_algn->len_query) + 1, MIN2(ga_params->ws, ga_algn->len_target) + 1);

    ga_algn->cigar = (char *)malloc((ga_algn->len_query + ga_algn->len_target) * sizeof(char) * 2);
    if (ga_algn->cigar == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
    }
    ga_algn->cigar_len = 0;
    Windowed_positions win_pos = {
        {ga_algn->len_query, ga_algn->len_target},                                                    // qt_bottom
        {MAX2(0, ga_algn->len_query - ga_params->ws), MAX2(0, ga_algn->len_target - ga_params->ws)},  // qt_top
        {MIN2(ga_params->ws, ga_algn->len_query), MIN2(ga_params->ws, ga_algn->len_target)},          // w_bottom
        {MIN2(ga_params->ws, ga_algn->len_query), MIN2(ga_params->ws, ga_algn->len_target)},          // last_stop
        -1                                                                                            // last matrix, -1 for undef
    };
    while (win_pos.last_stop.i > 0 && win_pos.last_stop.j > 0) {
        
        windowed_compute_window(ga_algn, ga_params, ga_res, win_pos);
        
        windowed_local_backtrace(ga_algn, ga_params, ga_res, &win_pos);

        windowed_update_positions(ga_algn, ga_params, ga_res, &win_pos);
    }

    if (DEBUG == 1) {
        printf(" Windowed Gap-Affine:       %5d\n", ga_res->score);
    }

    free(ga_algn->cigar);

    int end_time = clock();
    ga_res->elapsed = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000.0;
    ga_res->memory = get_memory_usage() - ga_res->memory;
    
    // Free matrices
    for (int i = 0; i <= MIN2(ga_params->ws, ga_algn->len_query); i++) {
        free(ga_algn->M[i]);
        free(ga_algn->I[i]);
        free(ga_algn->D[i]);
    }
    free(ga_algn->M);
    free(ga_algn->I);
    free(ga_algn->D);
}
