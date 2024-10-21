#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>  // Include this header for DBL_MAX
#include <time.h>
#include <string.h>
#include "main.h"


#define MAX_LEN 20000  // Define a maximum sequence length, adjust as needed
#define MIN2(a,b) a < b ? a : b
#define MIN3(a,b,c) a < b ? MIN2(a,c) : MIN2(b,c)

// Matrix Initialization
double **M, **I, **D;

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

// void GapAffine(const char *Q, const char *T, int Cm, int Cx, int Co, int Ci, int Cd, int *score, int *memory, double *elapsed, char *cigar, int *cigar_len) {
void GapAffine(const char *Q, const char *T, int Cm, int Cx, int Co, int Ci, int Cd, int *score, int *memory, double *elapsed) {

    int start_time = clock();
    // *memory = get_memory_usage();

    const int len_query = strlen(Q);
    const int len_target = strlen(T);


    M = create_matrix(len_query + 1, len_target + 1, DBL_MAX);
    I = create_matrix(len_query + 1, len_target + 1, DBL_MAX);
    D = create_matrix(len_query + 1, len_target + 1, DBL_MAX);
    // // Allocate memory for matrices M, D, and I
    // double **M = (double **)malloc((len_query + 1) * sizeof(double *));
    // double **D = (double **)malloc((len_query + 1) * sizeof(double *));
    // double **I = (double **)malloc((len_query + 1) * sizeof(double *));
    // for (int i = 0; i <= len_query; i++) {
    //     M[i] = (double *)malloc((len_target + 1) * sizeof(double));
    //     D[i] = (double *)malloc((len_target + 1) * sizeof(double));
    //     I[i] = (double *)malloc((len_target + 1) * sizeof(double));
    // }

    // Initialization
    for (int i = 1; i <= len_query; i++) {
        M[i][0] = Co + Cd * i;
        // I[i][0] = DBL_MAX;
        D[i][0] = Co + Cd * i;
    }
    for (int j = 1; j <= len_target; j++) {
        M[0][j] = Co + Ci * j;
        I[0][j] = Co + Ci * j;
        // D[0][j] = DBL_MAX;
    }

    // Align Score
    for (int i = 1; i <= len_query; i++) {
        for (int j = 1; j <= len_target; j++) {
            M[i][j] = MIN3(M[i-1][j-1], I[i-1][j-1], D[i-1][j-1]) + ((Q[i-1] == T[j-1]) ? Cm : Cx);
            I[i][j] = MIN2(I[i-1][j], M[i-1][j] + Co) + Ci;
            D[i][j] = MIN2(D[i][j-1], M[i][j-1] + Co) + Cd;
        }
    }
    *score = MIN3(M[len_query][len_target], I[len_query][len_target], D[len_query][len_target]);
    int end_time = clock();
    *elapsed = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000.0;
    // *memory = get_memory_usage() - *memory;
    // Free allocated memory
    for (int i = 0; i <= len_query; i++) {
        free(M[i]);
        free(D[i]);
        free(I[i]);
    }
    free(M);
    free(D);
    free(I);
    // backtrace_bucle(Q, T, Cm, Cx, Co, Ci, Cd, 0,0,len_query,len_target, cigar, cigar_len);
}


// Backtrace function to create CIGAR string
void backtrace_bucle(const char *Q, const char *T, const int *Cm, const int *Cx, const int *Co, const int *Ci, const int *Cd,
                     int end_i, int end_j, int *i, int *j, char *cigar, int *cigar_len) {
    int current_table;
    double current_value;

    if (M[*i][*j] <= I[*i][*j] && M[*i][*j] >= D[*i][*j]) {
        current_value = M[*i][*j];
        current_table = 'M';
    } else if (I[*i][*j] <= D[*i][*j]) {
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
            if (current_value == M[*i-1][*j] + *Co + *Ci) {
                current_value = M[*i-1][*j];
                current_table = 'M';
            } else if (current_value == I[*i-1][*j] + *Ci) {
                current_value = I[*i-1][*j];
                current_table = 'I';
            }
            cigar[(*cigar_len)++] = 'D';
            (*i)--;
        } else if (current_table == 'D') {
            if (current_value == M[*i][*j-1] + *Co + *Cd) {
                current_value = M[*i][*j-1];
                current_table = 'M';
            } else if (current_value == D[*i][*j-1] + *Cd) {
                current_value = D[*i][*j-1];
                current_table = 'D';
            }
            cigar[(*cigar_len)++] = 'I';
            (*j)--;
        }
    }
}

// int main() {
//     char *query =  "CAGATCGGATACCATCCATATGCCTTT";
//     char *target = "CAGATCGATACCATCCATTTGCCTT";
//     // int ws = 1500;
//     // int os = 50;
//     int Cm = 0;
//     int Cx = 3;
//     int Co = 4;
//     int Ci = 3;
//     int Cd = 3;
//     int score = 0, memory = 0;
//     double elapsed= 0;

//     char *cigar = (char *)malloc(2 * MAX_LEN * sizeof(char));
//     if (cigar == NULL) {
//         fprintf(stderr, "Memory allocation failed\n");
//     }
//     int cigar_len = 0;

//     GapAffine(query, target, Cm, Cx, Co, Ci, Cd, &score, &memory, &elapsed, cigar, &cigar_len);
//     printf("Gap-Affine real score:         %d   \n", score);
// }