#include "GapAffine_Windowed_BoundAndAlign.h"

// Function to update the window bounds
void update_window_bound(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res,
                        int start_i, int start_j, int end_i, int end_j, int is_last_window) {

    ga_algn->M[start_i][start_j] = 0;
    ga_algn->I[start_i][start_j] = 0;
    ga_algn->D[start_i][start_j] = 0;

    for (int i = start_i; i <= end_i; i++) {
        if (i != start_i) {
            ga_algn->M[i][start_j] = (is_last_window ? (ga_params->Co + ga_params->Cd * i) : 0);
            ga_algn->I[i][start_j] = __UINT16_MAX__;
            ga_algn->D[i][start_j] = (is_last_window ? (ga_params->Co + ga_params->Cd * i) : 0);
        }

        for (int j = start_j+1; j <= end_j; j++) {
            if (i == start_i) {
                ga_algn->M[i][j] = (is_last_window ? (ga_params->Co + ga_params->Ci * j) : 0);
                ga_algn->I[i][j] = (is_last_window ? (ga_params->Co + ga_params->Ci * j) : 0);
                ga_algn->D[i][j] = __UINT16_MAX__;
            } else {
                int C = (ga_algn->query[i-1] == ga_algn->target[j-1] ? ga_params->Cm : ga_params->Cx);
                ga_algn->M[i][j] = MIN3(ga_algn->M[i-1][j-1], ga_algn->I[i-1][j-1], ga_algn->D[i-1][j-1]) + C;
                ga_algn->I[i][j] = MIN2(ga_algn->I[i-1][j], ga_algn->M[i-1][j] + ga_params->Co) + ga_params->Ci;
                ga_algn->D[i][j] = MIN2(ga_algn->D[i][j-1], ga_algn->M[i][j-1] + ga_params->Co) + ga_params->Cd;
            }
            ga_res->computed_cells_windowed += 1;
        }
    }
}


// Backtrace function to create CIGAR string
void backtrace_bucle(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params,
                     int end_i, int end_j, int *i, int *j) {
    int current_table;
    int current_value;
    if (MIN3(ga_algn->M[*i][*j], ga_algn->I[*i][*j], ga_algn->D[*i][*j]) == ga_algn->M[*i][*j]) {
        current_value = ga_algn->M[*i][*j];
        current_table = 'M';
    } else if (MIN2(ga_algn->I[*i][*j],ga_algn->D[*i][*j]) == ga_algn->I[*i][*j]) {
        current_value = ga_algn->I[*i][*j];
        current_table = 'I';
    } else {
        current_value = ga_algn->D[*i][*j];
        current_table = 'D';
    }

    while ((*i > end_i && *j > end_j) || (*i == 0 && *j > 0) || (*j == 0 && *i > 0)) {
        if (*i == 0) {
            ga_algn->cigar[(ga_algn->cigar_len)++] = 'I';
            (*j)--;
        } else if (*j == 0) {
            ga_algn->cigar[(ga_algn->cigar_len)++] = 'D';
            (*i)--;
        } else if (current_table == 'M') {
            double C = (ga_algn->query[*i-1] == ga_algn->target[*j-1]) ? ga_params->Cm : ga_params->Cx;
            if (current_value == ga_algn->M[*i-1][*j-1] + C) {
                current_value = ga_algn->M[*i-1][*j-1];
            } else if (current_value == ga_algn->I[*i-1][*j-1] + C) {
                current_table = 'I';
                current_value = ga_algn->I[*i-1][*j-1];
            } else if (current_value == ga_algn->D[*i-1][*j-1] + C) {
                current_table = 'D';
                current_value = ga_algn->D[*i-1][*j-1];
            }
            ga_algn->cigar[(ga_algn->cigar_len)++] = ((ga_algn->query[*i-1] == ga_algn->target[*j-1]) ? 'M' : 'X');
            (*i)--;
            (*j)--;
        } else if (current_table == 'I') {
            if (current_value == ga_algn->M[*i-1][*j] + ga_params->Co + ga_params->Ci) {
                current_value = ga_algn->M[*i-1][*j];
                current_table = 'M';
            } else if (current_value == ga_algn->I[*i-1][*j] + ga_params->Ci) {
                current_value = ga_algn->I[*i-1][*j];
                current_table = 'I';
            }
            ga_algn->cigar[(ga_algn->cigar_len)++] = 'D';
            (*i)--;
        } else if (current_table == 'D') {
            if (current_value == ga_algn->M[*i][*j-1] + ga_params->Co + ga_params->Cd) {
                current_value = ga_algn->M[*i][*j-1];
                current_table = 'M';
            } else if (current_value == ga_algn->D[*i][*j-1] + ga_params->Cd) {
                current_value = ga_algn->D[*i][*j-1];
                current_table = 'D';
            }
            ga_algn->cigar[(ga_algn->cigar_len)++] = 'I';
            (*j)--;
        }
    }
}
// Function to calculate the CIGAR bound
void calculate_cigar_bound(GapAffine_Parameters *ga_params, GapAffine_Alignment *ga_algn, GapAffine_Results *ga_res) {
    ga_res->bound = 0;
    for (int i = 0; i < ga_algn->cigar_len; i++) {
        switch (ga_algn->cigar[i]) {
            case 'M': ga_res->bound += ga_params->Cm; break;
            case 'X': ga_res->bound += ga_params->Cx; break;
            case 'I': ga_res->bound += (i == 0 || ga_algn->cigar[i-1] != 'I') ? (ga_params->Ci + ga_params->Co) : ga_params->Ci; break;
            case 'D': ga_res->bound += (i == 0 || ga_algn->cigar[i-1] != 'D') ? (ga_params->Cd + ga_params->Co) : ga_params->Cd; break;
        }
    }
}

// Function to calculate the CIGAR score
void calculate_cigar_score(GapAffine_Parameters *ga_params, GapAffine_Alignment *ga_algn, GapAffine_Results *ga_res) {

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

// Function to build and print the CIGAR string from operations
char* print_cigar_windowed(GapAffine_Alignment *ga_algn) {
    // Dynamically allocate memory for formatted_cigar
    char *formatted_cigar = (char *)malloc(2 * ga_algn->cigar_len * sizeof(char));
    if (formatted_cigar == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
    }
    int pos = 0;

    char current_op = ga_algn->cigar[0];
    int count = 0;
    for (int i = 0; i < ga_algn->cigar_len; i++) {
        if (ga_algn->cigar[i] == current_op) {
            count++;
        } else {
            pos += sprintf(&formatted_cigar[pos], "%d%c", count, current_op);
            current_op = ga_algn->cigar[i];
            count = 1;
        }
    }
    pos += sprintf(&formatted_cigar[pos], "%d%c", count, current_op);
    formatted_cigar[pos] = '\0';
    
    return formatted_cigar;
}

void reverse_cigar(GapAffine_Alignment *ga_algn) {
    for (int i = 0; i < ga_algn->cigar_len / 2; i++) {
        char temp = ga_algn->cigar[i];
        ga_algn->cigar[i] = ga_algn->cigar[ga_algn->cigar_len - i - 1];
        ga_algn->cigar[ga_algn->cigar_len - i - 1] = temp;
    }

}

// Main function to handle windowed gap-affine bound
void windowed_gapAffine_bound(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res) {
    int current_i = ga_algn->len_query;
    int current_j = ga_algn->len_target;
    int is_last_window = 0;

    int start_i = (current_i - ga_params->ws > 0) ? current_i - ga_params->ws : 0;
    int start_j = (current_j - ga_params->ws > 0) ? current_j - ga_params->ws : 0;
    int end_i = current_i;
    int end_j = current_j;
    while (start_i >= 0 && start_j >= 0 && !is_last_window) {
        is_last_window = (start_i == 0 && start_j == 0);
        
        update_window_bound(ga_algn, ga_params, ga_res, start_i, start_j, end_i, end_j, is_last_window);
        int temp_end_i, temp_end_j;
        if (is_last_window) {
            temp_end_i = (start_i > 0) ? start_i : 0;
            temp_end_j = (start_j > 0) ? start_j : 0;
        } else {
            temp_end_i = (end_i - (ga_params->ws - ga_params->os) > 0) ? end_i - (ga_params->ws - ga_params->os) : 0;
            temp_end_j = (end_j - (ga_params->ws - ga_params->os) > 0) ? end_j - (ga_params->ws - ga_params->os) : 0;
        }
        backtrace_bucle(ga_algn, ga_params, temp_end_i, temp_end_j, &current_i, &current_j);

        start_i = (current_i - ga_params->ws > 0) ? current_i - ga_params->ws : 0;
        start_j = (current_j - ga_params->ws > 0) ? current_j - ga_params->ws : 0;
        end_i = current_i;
        end_j = current_j;
    }

    for (int i = 0; i < ga_algn->cigar_len / 2; i++) {
        char temp = ga_algn->cigar[i];
        ga_algn->cigar[i] = ga_algn->cigar[ga_algn->cigar_len - i - 1];
        ga_algn->cigar[ga_algn->cigar_len - i - 1] = temp;
    }
    calculate_cigar_bound(ga_params, ga_algn, ga_res);
    if (DEBUG == 1) {
        char* formatted_cigar = print_cigar_windowed(ga_algn);

        printf("Windowed Gap-Affine upper_bound:   %d    %s\n", ga_res->bound, formatted_cigar);
    }
}

void windowed_gapAffine_align(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res) {
    ga_res->computed_cells_banded = 0;
    // Initialize the matrices M, I, and D
    ga_algn->M[0][0] = 0;  
    ga_algn->I[0][0] = 0;
    ga_algn->D[0][0] = 0;
    for (int i = 1; i <= ga_algn->len_query; i++) {
            ga_algn->M[i][0] = ga_params->Co + ga_params->Cd * i;
            ga_algn->I[i][0] = __UINT16_MAX__;
            ga_algn->D[i][0] = ga_params->Co + ga_params->Cd * i;
    }
    for (int j = 1; j <= ga_algn->len_target; j++) {
            ga_algn->M[0][j] = ga_params->Co + ga_params->Ci * j;
            ga_algn->I[0][j] = ga_params->Co + ga_params->Ci * j;
            ga_algn->D[0][j] = __UINT16_MAX__;
    }
    int bottom_j = 1;
    // Align Score Calculation
    for (int i = 1; i <= ga_algn->len_query; i++) {
        for (int j = bottom_j; j <= ga_algn->len_target; j++) {
            int Cost = (ga_algn->query[i - 1] == ga_algn->target[j - 1]) ? ga_params->Cm : ga_params->Cx;

            ga_algn->M[i][j] = MIN3(ga_algn->M[i - 1][j - 1], ga_algn->I[i - 1][j - 1], ga_algn->D[i - 1][j - 1]) + Cost;
            ga_algn->I[i][j] = MIN2(ga_algn->I[i - 1][j] + ga_params->Ci, ga_algn->M[i - 1][j] + ga_params->Co + ga_params->Ci);
            ga_algn->D[i][j] = MIN2(ga_algn->D[i][j - 1] + ga_params->Cd, ga_algn->M[i][j - 1] + ga_params->Co + ga_params->Cd);

            ga_res->computed_cells_banded += 1;

            // Threshold check
            if (ga_algn->M[i][j] > ga_res->bound && ga_algn->I[i][j] > ga_res->bound && ga_algn->D[i][j] > ga_res->bound) {
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
    backtrace_bucle(ga_algn, ga_params, 0, 0, &ga_algn->len_query, &ga_algn->len_target);
    // ga_algn->len_query = (int) strlen(ga_algn->query);
    // ga_algn->len_target = (int) strlen(ga_algn->target);
    ga_algn->cigar_len = strlen(ga_algn->cigar);

    reverse_cigar(ga_algn);

    calculate_cigar_score(ga_params, ga_algn, ga_res);
    if (DEBUG == 1) {
        char* formatted_cigar = print_cigar_windowed(ga_algn);

        printf("Windowed Gap-Affine align:      %d    %s\n", ga_res->score, formatted_cigar);
        printf("Windowed Gap-Affine original:   %d    %s\n", ga_res->original_score, formatted_cigar);
    }
    // char* formatted_cigar = print_cigar_windowed(ga_algn, ga_res->score);
    // printf("Windowed Gap-Affine align:         %d    %s\n", ga_res->score, formatted_cigar);
}


// void GapAffine_windowed(const char *Q, const char *T, const int *ws, const int *os, int *Cm, int *Cx, int *Co, int *Ci, int *Cd, 
//                 int *score, int *bound, int *memory, double *elapsed, int alpha, int beta, int gamma) {
void GapAffine_windowed(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res) {
    
    int start_time = clock();
    ga_res->memory = get_memory_usage();
    // Initialize matrices
    
    ga_algn->M = create_matrix(ga_algn->len_query + 1, ga_algn->len_target + 1);
    ga_algn->I = create_matrix(ga_algn->len_query + 1, ga_algn->len_target + 1);
    ga_algn->D = create_matrix(ga_algn->len_query + 1, ga_algn->len_target + 1);

    ga_algn->cigar = (char *)malloc((ga_algn->len_query + ga_algn->len_target) * sizeof(char) * 2);
    if (ga_algn->cigar == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
    }
    ga_algn->cigar_len = 0;
    
    windowed_gapAffine_bound(ga_algn, ga_params, ga_res);
    free(ga_algn->cigar);
    ga_algn->cigar = (char *)malloc((ga_algn->len_query + ga_algn->len_target) * sizeof(char) * 2);
    if (ga_algn->cigar == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
    }
    ga_algn->cigar_len = 0;
    
    windowed_gapAffine_align(ga_algn, ga_params, ga_res);
    int end_time = clock();
    ga_res->elapsed = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000.0;
    ga_res->memory = get_memory_usage() - ga_res->memory;
}
