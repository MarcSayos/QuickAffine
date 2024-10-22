#include "main.h"


// Matrix Initialization
__uint32_t **M, **I, **D;

void GapAffine(const char *Q, const char *T, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res) {

    int start_time = clock();
    ga_res->memory = get_memory_usage();

    const int len_query = strlen(Q);
    const int len_target = strlen(T);

    M = create_matrix(len_query + 1, len_target + 1);
    I = create_matrix(len_query + 1, len_target + 1);
    D = create_matrix(len_query + 1, len_target + 1);

    // Initialization
    for (int i = 1; i <= len_query; i++) {
        M[i][0] = ga_params->Co + ga_params->Cd * i;
        I[i][0] = __UINT16_MAX__;
        D[i][0] = ga_params->Co + ga_params->Cd * i;
    }
    for (int j = 1; j <= len_target; j++) {
        M[0][j] = ga_params->Co + ga_params->Ci * j;
        I[0][j] = ga_params->Co + ga_params->Ci * j;
        D[0][j] = __UINT16_MAX__;
    }

    // Align Score
    for (int i = 1; i <= len_query; i++) {
        for (int j = 1; j <= len_target; j++) {
            M[i][j] = MIN3(M[i-1][j-1], I[i-1][j-1], D[i-1][j-1]) + ((Q[i-1] == T[j-1]) ? ga_params->Cm : ga_params->Cx);
            I[i][j] = MIN2(I[i-1][j], M[i-1][j] + ga_params->Co) + ga_params->Ci;
            D[i][j] = MIN2(D[i][j-1], M[i][j-1] + ga_params->Co) + ga_params->Cd;
        }
    }

    ga_res->score = MIN3(M[len_query][len_target], I[len_query][len_target], D[len_query][len_target]);
    int end_time = clock();
    ga_res->elapsed = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000.0;
    ga_res->memory = get_memory_usage() - ga_res->memory;


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