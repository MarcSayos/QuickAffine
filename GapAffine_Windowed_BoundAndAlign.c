#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>  // Include this header for DBL_MAX
#include <math.h>
#include "main.h"
#include <sys/time.h>
#include <sys/resource.h>

#define MAX_LEN 20000  // Define a maximum sequence length, adjust as needed
#define max(a,b) \
({ __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b; })
    
// Matrix Initialization
double **M, **I, **D;

int debug = 0;

// Function to create a matrix with initial values
double **create_matrix(int rows, int cols, double initial_value) {
    double **matrix = (double **)malloc(rows * sizeof(double *));
    for (int i = 0; i < rows; i++) {
        matrix[i] = (double *)malloc(cols * sizeof(double));
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = initial_value;
        }
    }
    return matrix;
}


void print_memory_usage() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    
    // Memory usage (in kilobytes)
    printf("Memory usage: %ld KB\n", usage.ru_maxrss);
}

// Function to update the window bounds
void update_window_bound(const char *Q, const char *T, const int *Cm, const int *Cx, const int *Co, const int *Ce,
                        int start_i, int start_j, int end_i, int end_j, int is_last_window) {
    for (int i = start_i; i <= end_i; i++) {
        for (int j = start_j; j <= end_j; j++) {
            if (i == start_i && j == start_j) {
                M[i][j] = 0;
                I[i][j] = 0;
                D[i][j] = 0;
            } else if (i == start_i && j != start_j) {
                M[i][j] = is_last_window ? (*Co + *Ce * j) : 0;
                I[i][j] = is_last_window ? (*Co + *Ce * j) : 0;
                D[i][j] = -DBL_MAX;
            } else if (j == start_j && i != start_i) {
                M[i][j] = is_last_window ? (*Co + *Ce * i) : 0;
                I[i][j] = -DBL_MAX;
                D[i][j] = is_last_window ? (*Co + *Ce * i) : 0;
            } else {
                M[i][j] = max(max(M[i-1][j-1], I[i-1][j-1]), D[i-1][j-1]) + (Q[i-1] == T[j-1] ? *Cm : *Cx);
                I[i][j] = max(I[i-1][j], M[i-1][j] + *Co) + *Ce;
                D[i][j] = max(D[i][j-1], M[i][j-1] + *Co) + *Ce;
            }
        }
    }
}

// Backtrace function to create CIGAR string
void backtrace_bucle(const char *Q, const char *T, const int *Cm, const int *Cx, const int *Co, const int *Ce,
                     int end_i, int end_j, int *i, int *j, char *cigar, int *cigar_len) {
    int current_table;
    double current_value;

    if (M[*i][*j] >= I[*i][*j] && M[*i][*j] >= D[*i][*j]) {
        current_value = M[*i][*j];
        current_table = 'M';
    } else if (I[*i][*j] >= D[*i][*j]) {
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
            double C = (Q[*i-1] == T[*j-1]) ? *Cm : *Cx;
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
            if (current_value == M[*i-1][*j] + *Co + *Ce) {
                current_value = M[*i-1][*j];
                current_table = 'M';
            } else if (current_value == I[*i-1][*j] + *Ce) {
                current_value = I[*i-1][*j];
                current_table = 'I';
            }
            cigar[(*cigar_len)++] = 'D';
            (*i)--;
        } else if (current_table == 'D') {
            if (current_value == M[*i][*j-1] + *Co + *Ce) {
                current_value = M[*i][*j-1];
                current_table = 'M';
            } else if (current_value == D[*i][*j-1] + *Ce) {
                current_value = D[*i][*j-1];
                current_table = 'D';
            }
            cigar[(*cigar_len)++] = 'I';
            (*j)--;
        }
    }
}

// Function to calculate the CIGAR score
int calculate_cigar_score(const int *Cm, const int *Cx, const int *Co, const int *Ce,char *cigar, int cigar_len) {
    int score = 0;
    for (int i = 0; i < cigar_len; i++) {
        switch (cigar[i]) {
            case 'M':
                score += *Cm;
                break;
            case 'X':
                score += *Cx;
                break;
            case 'I':
                score += (i == 0 || cigar[i-1] != 'I') ? (*Ce + *Co) : *Ce;
                break;
            case 'D':
                score += (i == 0 || cigar[i-1] != 'D') ? (*Ce + *Co) : *Ce;
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
int windowed_gapAffine_bound(const char *Q, const char *T, const int * window_size, const int *overlap_size, 
                    const int *Cm, const int *Cx, const int *Co, const int *Ce,
                    char *cigar, int *cigar_len, int *len_query, int *len_target) {
    int current_i = *len_query;
    int current_j = *len_target;
    int is_last_window = 0;

    int start_i = (current_i - *window_size > 0) ? current_i - *window_size : 0;
    int start_j = (current_j - *window_size > 0) ? current_j - *window_size : 0;
    int end_i = current_i;
    int end_j = current_j;
    while (start_i >= 0 && start_j >= 0 && !is_last_window) {
        is_last_window = (start_i == 0 && start_j == 0);

        update_window_bound(Q, T, Cm, Cx, Co, Ce, start_i, start_j, end_i, end_j, is_last_window);
        int temp_end_i, temp_end_j;
        if (is_last_window) {
            temp_end_i = (start_i > 0) ? start_i : 0;
            temp_end_j = (start_j > 0) ? start_j : 0;
        } else {
            temp_end_i = (end_i - (*window_size - *overlap_size) > 0) ? end_i - (*window_size - *overlap_size) : 0;
            temp_end_j = (end_j - (*window_size - *overlap_size) > 0) ? end_j - (*window_size - *overlap_size) : 0;
        }
        backtrace_bucle(Q, T, Cm, Cx, Co, Ce, temp_end_i, temp_end_j, &current_i, &current_j, cigar, cigar_len);

        start_i = (current_i - *window_size > 0) ? current_i - *window_size : 0;
        start_j = (current_j - *window_size > 0) ? current_j - *window_size : 0;
        end_i = current_i;
        end_j = current_j;
    }

    for (int i = 0; i < *cigar_len / 2; i++) {
        char temp = cigar[i];
        cigar[i] = cigar[*cigar_len - i - 1];
        cigar[*cigar_len - i - 1] = temp;
    }
    int score = calculate_cigar_score(Cm, Cx, Co, Ce, cigar, *cigar_len);
    if (debug == 1) {
        char* formatted_cigar = print_cigar_windowed(cigar, *cigar_len, score);

        printf("Windowed Gap-Affine upper_bound:   %d    %s\n", score, formatted_cigar);
    }
    return score;

}

int windowed_gapAffine_align(int bound, const char *Q, const char *T, const int * window_size, const int *overlap_size, 
                    const int *Cm, const int *Cx, const int *Co, const int *Ce,
                    char *cigar, int *cigar_len, int *len_query, int *len_target) {

    // Initialize the matrices M, I, and D
    for (int i = 0; i <= *len_query; i++) {
        for (int j = 0; j <= *len_target; j++) {
            if (i == 0 && j == 0) {
                M[i][j] = 0;  
                I[i][j] = -DBL_MAX;
                D[i][j] = -DBL_MAX;
            }
            else if (j == 0) {
                M[i][j] = *Co + *Ce * i;
                I[i][j] = -DBL_MAX;
                D[i][j] = *Co + *Ce * i;
            }
            else if (i == 0) {
                M[i][j] = *Co + *Ce * j;
                I[i][j] = *Co + *Ce * j;
                D[i][j] = -DBL_MAX;
            }
            else {
                M[i][j] = -DBL_MAX;  
                I[i][j] = -DBL_MAX;
                D[i][j] = -DBL_MAX;
            }
        }
    }
    int bottom_j = 1;
    // Align Score Calculation
    for (int i = 1; i <= *len_query; i++) {
        for (int j = bottom_j; j <= *len_target; j++) {
            int Cost = (Q[i - 1] == T[j - 1]) ? *Cm : *Cx;

            M[i][j] = max(max(M[i - 1][j - 1], I[i - 1][j - 1]), D[i - 1][j - 1]) + Cost;
            I[i][j] = max(I[i - 1][j] + *Ce, M[i - 1][j] + *Co + *Ce);
            D[i][j] = max(D[i][j - 1] + *Ce, M[i][j - 1] + *Co + *Ce);

            // Threshold check
            if (M[i][j] < bound && I[i][j] < bound && D[i][j] < bound) {
                if ((i-1 == 0) || M[i-1][j] == -DBL_MAX) { // Top
                    // i++;
                    M[i][j] = I[i][j] = D[i][j] = -DBL_MAX;
                    break;
                } else if ((j-1 == 0) || M[i][j-1] == -DBL_MAX) {
                    bottom_j++;
                    M[i][j] = I[i][j] = D[i][j] = -DBL_MAX;
                }
            }
        }
    }

    backtrace_bucle(Q, T, Cm, Cx, Co, Ce, 0, 0, len_query, len_target, cigar, cigar_len);
    for (int i = 0; i < *cigar_len / 2; i++) {
        char temp = cigar[i];
        cigar[i] = cigar[*cigar_len - i - 1];
        cigar[*cigar_len - i - 1] = temp;
    }

    int score = calculate_cigar_score(Cm, Cx, Co, Ce, cigar, *cigar_len);
    if (debug == 1) {
        char* formatted_cigar = print_cigar_windowed(cigar, *cigar_len, score);

        printf("Windowed Gap-Affine align:         %d    %s\n", score, formatted_cigar);
    }
    return score;
}


void GapAffine_windowed(const char *Q, const char *T, const int *ws, const int *os, const int *Cm, const int *Cx, const int *Co, const int *Ce, 
                int *score, int *memory, double *elapsed) {
    
    int start_time = clock();
    *memory = get_memory_usage();

    int len_query = strlen(Q);
    int len_target = strlen(T);
    // Initialize matrices
    
    M = create_matrix(len_query + 1, len_target + 1, -DBL_MAX);
    I = create_matrix(len_query + 1, len_target + 1, -DBL_MAX);
    D = create_matrix(len_query + 1, len_target + 1, -DBL_MAX);

    char *cigar = (char *)malloc(2 * MAX_LEN * sizeof(char));
    if (cigar == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
    }
    int cigar_len = 0;

    int bound = windowed_gapAffine_bound(Q, T, ws, os, Cm, Cx, Co, Ce, cigar, &cigar_len, &len_query, &len_target);
    free(cigar);
    cigar = (char *)malloc(2 * MAX_LEN * sizeof(char));
    if (cigar == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
    }
    cigar_len = 0;

    *score = windowed_gapAffine_align(bound, Q, T, ws, os, Cm, Cx, Co, Ce, cigar, &cigar_len, &len_query, &len_target);
    int end_time = clock();
    *elapsed = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000.0;
    *memory = get_memory_usage() - *memory;

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
//     char *query = "ATCAACTAGAGCAAAACTTTTGAAGTGGAGGACACGTAGGCCTTCTTGAGC";
//     char *target = "ATCAAGTAGAGCAAAACTTTTAAGTTGGAGGACACGTAGGCCTCGTGAGC";
//     int ws = 1500;
//     int os = 50;
//     int Cm = 0;
//     int Cx = -3;
//     int Co = -4;
//     int Ce = -3;
//     int score = main_windowed(query, target, &ws, &os, &Cm, &Cx, &Co, &Ce);
// }