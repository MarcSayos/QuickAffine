
#include "main.h"
    
// Matrix Initialization
double **M, **I, **D;

int debug = 0;

// Function to update the window bounds
void update_window_bound(const char *Q, const char *T, GapAffine_Parameters *ga_params,
                        int start_i, int start_j, int end_i, int end_j, int is_last_window) {
    for (int i = start_i; i <= end_i; i++) {
        for (int j = start_j; j <= end_j; j++) {
            if (i == start_i && j == start_j) {
                M[i][j] = 0;
                I[i][j] = 0;
                D[i][j] = 0;
            } else if (i == start_i && j != start_j) {
                M[i][j] = is_last_window ? (ga_params->Co + ga_params->Ci * j) : 0;
                I[i][j] = is_last_window ? (ga_params->Co + ga_params->Ci * j) : 0;
                D[i][j] = DBL_MAX;
            } else if (j == start_j && i != start_i) {
                M[i][j] = is_last_window ? (ga_params->Co + ga_params->Cd * i) : 0;
                I[i][j] = DBL_MAX;
                D[i][j] = is_last_window ? (ga_params->Co + ga_params->Cd * i) : 0;
            } else {
                M[i][j] = MIN3(M[i-1][j-1], I[i-1][j-1], D[i-1][j-1]) + (Q[i-1] == T[j-1] ? ga_params->Cm : ga_params->Cx);
                I[i][j] = MIN2(I[i-1][j], M[i-1][j] + ga_params->Co) + ga_params->Ci;
                D[i][j] = MIN2(D[i][j-1], M[i][j-1] + ga_params->Co) + ga_params->Cd;
            }
        }
    }
}

// Backtrace function to create CIGAR string
void backtrace_bucle(const char *Q, const char *T, GapAffine_Parameters *ga_params,
                     int end_i, int end_j, int *i, int *j, char *cigar, int *cigar_len) {
    int current_table;
    double current_value;

    if (MIN3(M[*i][*j], I[*i][*j], D[*i][*j]) == M[*i][*j]) {
        current_value = M[*i][*j];
        current_table = 'M';
    } else if (MIN2(I[*i][*j],D[*i][*j]) == I[*i][*j]) {
        current_value = I[*i][*j];
        current_table = 'I';
    } else {
        current_value = D[*i][*j];
        current_table = 'D';
    }

    while ((*i > end_i && *j > end_j) || (*i == 0 && *j > 0) || (*j == 0 && *i > 0)) {
        if (*i == 0) {
            cigar[(*cigar_len)++] = 'I';
            (*j)--;
        } else if (*j == 0) {
            cigar[(*cigar_len)++] = 'D';
            (*i)--;
        } else if (current_table == 'M') {
            double C = (Q[*i-1] == T[*j-1]) ? ga_params->Cm : ga_params->Cx;
            if (current_value == M[*i-1][*j-1] + C) {
                current_value = M[*i-1][*j-1];
            } else if (current_value == I[*i-1][*j-1] + C) {
                current_table = 'I';
                current_value = I[*i-1][*j-1];
            } else if (current_value == D[*i-1][*j-1] + C) {
                current_table = 'D';
                current_value = D[*i-1][*j-1];
            }
            cigar[(*cigar_len)++] = (Q[*i-1] == T[*j-1]) ? 'M' : 'X';
            (*i)--;
            (*j)--;
        } else if (current_table == 'I') {
            if (current_value == M[*i-1][*j] + ga_params->Co + ga_params->Ci) {
                current_value = M[*i-1][*j];
                current_table = 'M';
            } else if (current_value == I[*i-1][*j] + ga_params->Ci) {
                current_value = I[*i-1][*j];
                current_table = 'I';
            }
            cigar[(*cigar_len)++] = 'D';
            (*i)--;
        } else if (current_table == 'D') {
            if (current_value == M[*i][*j-1] + ga_params->Co + ga_params->Cd) {
                current_value = M[*i][*j-1];
                current_table = 'M';
            } else if (current_value == D[*i][*j-1] + ga_params->Cd) {
                current_value = D[*i][*j-1];
                current_table = 'D';
            }
            cigar[(*cigar_len)++] = 'I';
            (*j)--;
        }
    }
}

// Function to calculate the CIGAR score
int calculate_cigar_score(GapAffine_Parameters *ga_params, char *cigar, int cigar_len) {
    int score = 0;
    for (int i = 0; i < cigar_len; i++) {
        switch (cigar[i]) {
            case 'M':
                score += ga_params->Cm;
                break;
            case 'X':
                score += ga_params->Cx;
                break;
            case 'I':
                score += (i == 0 || cigar[i-1] != 'I') ? (ga_params->Ci + ga_params->Co) : ga_params->Ci;
                break;
            case 'D':
                score += (i == 0 || cigar[i-1] != 'D') ? (ga_params->Cd + ga_params->Co) : ga_params->Cd;
                break;
        }
    }
    return score;
}

// Function to build and print the CIGAR string from operations
char* print_cigar_windowed(char *cigar_ops, int cigar_len, int score) {
    // Dynamically allocate memory for formatted_cigar
    char *formatted_cigar = (char *)malloc(2 * cigar_len * sizeof(char));
    if (formatted_cigar == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
    }
    int pos = 0;

    char current_op = cigar_ops[0];
    int count = 0;
    for (int i = 0; i < cigar_len; i++) {
        if (cigar_ops[i] == current_op) {
            count++;
        } else {
            pos += sprintf(&formatted_cigar[pos], "%d%c", count, current_op);
            current_op = cigar_ops[i];
            count = 1;
        }
    }
    pos += sprintf(&formatted_cigar[pos], "%d%c", count, current_op);
    formatted_cigar[pos] = '\0';
    
    return formatted_cigar;
}

// Main function to handle windowed gap-affine bound
int windowed_gapAffine_bound(const char *Q, const char *T, GapAffine_Parameters *ga_params, char *cigar, int *cigar_len, int *len_query, int *len_target) {
    int current_i = *len_query;
    int current_j = *len_target;
    int is_last_window = 0;

    int start_i = (current_i - ga_params->ws > 0) ? current_i - ga_params->ws : 0;
    int start_j = (current_j - ga_params->ws > 0) ? current_j - ga_params->ws : 0;
    int end_i = current_i;
    int end_j = current_j;
    while (start_i >= 0 && start_j >= 0 && !is_last_window) {
        is_last_window = (start_i == 0 && start_j == 0);

        update_window_bound(Q, T, ga_params, start_i, start_j, end_i, end_j, is_last_window);
        int temp_end_i, temp_end_j;
        if (is_last_window) {
            temp_end_i = (start_i > 0) ? start_i : 0;
            temp_end_j = (start_j > 0) ? start_j : 0;
        } else {
            temp_end_i = (end_i - (ga_params->ws - ga_params->os) > 0) ? end_i - (ga_params->ws - ga_params->os) : 0;
            temp_end_j = (end_j - (ga_params->ws - ga_params->os) > 0) ? end_j - (ga_params->ws - ga_params->os) : 0;
        }
        backtrace_bucle(Q, T, ga_params, temp_end_i, temp_end_j, &current_i, &current_j, cigar, cigar_len);

        start_i = (current_i - ga_params->ws > 0) ? current_i - ga_params->ws : 0;
        start_j = (current_j - ga_params->ws > 0) ? current_j - ga_params->ws : 0;
        end_i = current_i;
        end_j = current_j;
    }

    for (int i = 0; i < *cigar_len / 2; i++) {
        char temp = cigar[i];
        cigar[i] = cigar[*cigar_len - i - 1];
        cigar[*cigar_len - i - 1] = temp;
    }
    int score = calculate_cigar_score(ga_params, cigar, *cigar_len);
    if (debug == 1) {
        char* formatted_cigar = print_cigar_windowed(cigar, *cigar_len, score);

        printf("Windowed Gap-Affine upper_bound:   %d    %s\n", score, formatted_cigar);
    }
    return score;

}

int windowed_gapAffine_align(int bound, const char *Q, const char *T, GapAffine_Parameters *ga_params,
                    char *cigar, int *cigar_len, int *len_query, int *len_target) {

    // Initialize the matrices M, I, and D
    for (int i = 0; i <= *len_query; i++) {
        for (int j = 0; j <= *len_target; j++) {
            if (i == 0 && j == 0) {
                M[i][j] = 0;  
                I[i][j] = DBL_MAX;
                D[i][j] = DBL_MAX;
            }
            else if (j == 0) {
                M[i][j] = ga_params->Co + ga_params->Cd * i;
                I[i][j] = DBL_MAX;
                D[i][j] = ga_params->Co + ga_params->Cd * i;
            }
            else if (i == 0) {
                M[i][j] = ga_params->Co + ga_params->Ci * j;
                I[i][j] = ga_params->Co + ga_params->Ci * j;
                D[i][j] = DBL_MAX;
            }
            else {
                M[i][j] = DBL_MAX;  
                I[i][j] = DBL_MAX;
                D[i][j] = DBL_MAX;
            }
        }
    }
    int bottom_j = 1;
    // Align Score Calculation
    for (int i = 1; i <= *len_query; i++) {
        for (int j = bottom_j; j <= *len_target; j++) {
            int Cost = (Q[i - 1] == T[j - 1]) ? ga_params->Cm : ga_params->Cx;

            M[i][j] = MIN3(M[i - 1][j - 1], I[i - 1][j - 1], D[i - 1][j - 1]) + Cost;
            I[i][j] = MIN2(I[i - 1][j] + ga_params->Ci, M[i - 1][j] + ga_params->Co + ga_params->Ci);
            D[i][j] = MIN2(D[i][j - 1] + ga_params->Cd, M[i][j - 1] + ga_params->Co + ga_params->Cd);

            // Threshold check
            if (M[i][j] > bound && I[i][j] > bound && D[i][j] > bound) {
                if ((i-1 == 0) || M[i-1][j] == DBL_MAX) { // Top
                    // i++;
                    M[i][j] = I[i][j] = D[i][j] = DBL_MAX;
                    break;
                } else if ((j-1 == 0) || M[i][j-1] == DBL_MAX) {
                    bottom_j++;
                    M[i][j] = I[i][j] = D[i][j] = DBL_MAX;
                }
            }
        }
    }

    backtrace_bucle(Q, T, ga_params, 0, 0, len_query, len_target, cigar, cigar_len);
    for (int i = 0; i < *cigar_len / 2; i++) {
        char temp = cigar[i];
        cigar[i] = cigar[*cigar_len - i - 1];
        cigar[*cigar_len - i - 1] = temp;
    }

    int score = calculate_cigar_score(ga_params, cigar, *cigar_len);
    if (debug == 1) {
        char* formatted_cigar = print_cigar_windowed(cigar, *cigar_len, score);

        printf("Windowed Gap-Affine align:         %d    %s\n", score, formatted_cigar);
    }
    return score;
}


// void GapAffine_windowed(const char *Q, const char *T, const int *ws, const int *os, int *Cm, int *Cx, int *Co, int *Ci, int *Cd, 
//                 int *score, int *bound, int *memory, double *elapsed, int alpha, int beta, int gamma) {
void GapAffine_windowed(const char *Q, const char *T, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res) {
    
    int start_time = clock();
    ga_res->memory = get_memory_usage();

    int len_query = strlen(Q);
    int len_target = strlen(T);
    // Initialize matrices
    
    M = create_matrix(len_query + 1, len_target + 1);
    I = create_matrix(len_query + 1, len_target + 1);
    D = create_matrix(len_query + 1, len_target + 1);

    char *cigar = (char *)malloc(2 * MAX_LEN * sizeof(char));
    if (cigar == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
    }
    int cigar_len = 0;

    ga_res->bound = windowed_gapAffine_bound(Q, T, ga_params, cigar, &cigar_len, &len_query, &len_target);
    free(cigar);
    cigar = (char *)malloc(2 * MAX_LEN * sizeof(char));
    if (cigar == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
    }
    cigar_len = 0;
    
    ga_res->score = windowed_gapAffine_align(ga_res->bound, Q, T, ga_params, cigar, &cigar_len, &len_query, &len_target);
    int end_time = clock();
    ga_res->elapsed = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000.0;
    ga_res->memory = get_memory_usage() - ga_res->memory;

    // Free matrices
    for (int i = 0; i <= len_query; i++) {
        free(M[i]);
        free(I[i]);
        free(D[i]);
    }
    free(M);
    free(I);
    free(D);
    free(cigar);
}

// int main() {
//     char *query  = "CAGATCTGATACCATCCATATGCCTTT";
//     char *target = "CAGATCGATACCATCCATTTGCCTT";
//     // char *query = "ATCAACTAGAGCAAAACTTTTGAAGTGGAGGACACGTAGGCCTTCTTGAGC";
//     // char *target = "ATCAAGTAGAGCAAAACTTTTAAGTTGGAGGACACGTAGGCCTCGTGAGC";
//     int ws = 1500,os = 50;
//     int Cm = 0,Cx = 3,Co = 4,Ci = 3,Cd = 3;
//     int score = 0, memory = 0; 
//     double elapsed = 0;
//     int alpha = 0, beta = 0, gamma = 0;
//     GapAffine_windowed(query, target, &ws, &os, &Cm, &Cx, &Co, &Ci, &Cd, &score, &memory, &elapsed, alpha, beta, gamma);
//     printf("Score = %d\n", score);
// }