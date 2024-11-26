#include "GapAffine_Banded.h"

// Backtrace function to create CIGAR string
void banded_backtrace(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params,
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

// Function to build and print the CIGAR string from operations
char* banded_print_cigar(GapAffine_Alignment *ga_algn) {
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

void banded_reverse_cigar(GapAffine_Alignment *ga_algn) {
    for (int i = 0; i < ga_algn->cigar_len / 2; i++) {
        char temp = ga_algn->cigar[i];
        ga_algn->cigar[i] = ga_algn->cigar[ga_algn->cigar_len - i - 1];
        ga_algn->cigar[ga_algn->cigar_len - i - 1] = temp;
    }

}


void banded_GapAffine(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res) {

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
    ga_res->cells = 0;
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

            ga_res->cells += 1;

            // Threshold check
            if (ga_algn->M[i][j] > ga_params->upper_bound && ga_algn->I[i][j] > ga_params->upper_bound && ga_algn->D[i][j] > ga_params->upper_bound) {
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
    banded_backtrace(ga_algn, ga_params, 0, 0, &ga_algn->len_query, &ga_algn->len_target);
    ga_algn->cigar_len = strlen(ga_algn->cigar);

    banded_reverse_cigar(ga_algn);

    banded_calculate_cigar_score(ga_params, ga_algn, ga_res);
    if (DEBUG == 1) {
        char* formatted_cigar = banded_print_cigar(ga_algn);

        printf("Windowed Gap-Affine align:         %5d    %s\n", ga_res->score, formatted_cigar);
        printf("Windowed Gap-Affine original pens: %5d    %s\n", ga_res->original_score, formatted_cigar);
    }

    int end_time = clock();
    ga_res->elapsed = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000.0;
    ga_res->memory = get_memory_usage() - ga_res->memory;
}
