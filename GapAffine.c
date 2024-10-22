#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>  // Include this header for DBL_MAX
#include <time.h>
#include <string.h>
#include "main.h"


#define MAX_LEN 20000  // Define a maximum sequence length, adjust as needed
#define min(a,b) \
({ __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b; })

// Matrix Initialization
double **M, **I, **D;

void GapAffine(const char *Q, const char *T, int Cm, int Cx, int Co, int Ci, int Cd, int *score, int *memory, double *elapsed) {

    int start_time = clock();
    *memory = get_memory_usage();

    const int len_query = strlen(Q);
    const int len_target = strlen(T);


    M = create_matrix(len_query + 1, len_target + 1);
    I = create_matrix(len_query + 1, len_target + 1);
    D = create_matrix(len_query + 1, len_target + 1);
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
        I[i][0] = DBL_MAX;
        D[i][0] = Co + Cd * i;
    }
    for (int j = 1; j <= len_target; j++) {
        M[0][j] = Co + Ci * j;
        I[0][j] = Co + Ci * j;
        D[0][j] = DBL_MAX;
    }

    // Align Score
    for (int i = 1; i <= len_query; i++) {
        for (int j = 1; j <= len_target; j++) {
            M[i][j] = min(min(M[i-1][j-1], I[i-1][j-1]), D[i-1][j-1]) + ((Q[i-1] == T[j-1]) ? Cm : Cx);
            I[i][j] = min(I[i-1][j], M[i-1][j] + Co) + Ci;
            D[i][j] = min(D[i][j-1], M[i][j-1] + Co) + Cd;
        }
    }
    *score = min(min(M[len_query][len_target], I[len_query][len_target]), D[len_query][len_target]);
    int end_time = clock();
    *elapsed = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000.0;
    *memory = get_memory_usage() - *memory;
    // Free allocated memory
    for (int i = 0; i <= len_query; i++) {
        free(M[i]);
        free(D[i]);
        free(I[i]);
    }
    free(M);
    free(D);
    free(I);
}

// int main() {
//     char *query = "CAGATCTGATACCATCCATATGCCTTT";
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

//     GapAffine(query, target, Cm, Cx, Co, Ci, Cd, &score, &memory, &elapsed);
//     printf("Gap-Affine real score:         %d   \n", score);
// }