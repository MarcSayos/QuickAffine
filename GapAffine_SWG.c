#include "GapAffine_SWG.h"

void GapAffine(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res) {

    int start_time = clock();
    ga_res->memory = get_memory_usage();

    ga_algn->M = create_matrix(ga_algn->len_query + 1, ga_algn->len_target + 1);
    ga_algn->I = create_matrix(ga_algn->len_query + 1, ga_algn->len_target + 1);
    ga_algn->D = create_matrix(ga_algn->len_query + 1, ga_algn->len_target + 1);

    ga_algn->M[0][0] = 0;
    ga_algn->I[0][0] = 0;
    ga_algn->D[0][0] = 0;

    // Initialization
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

    // Align Score
    for (int i = 1; i <= ga_algn->len_query; i++) {
        for (int j = 1; j <= ga_algn->len_target; j++) {
            ga_algn->M[i][j] = MIN3(ga_algn->M[i-1][j-1], ga_algn->I[i-1][j-1], ga_algn->D[i-1][j-1]) + ((ga_algn->query[i-1] == ga_algn->target[j-1]) ? ga_params->Cm : ga_params->Cx);
            ga_algn->I[i][j] = MIN2(ga_algn->I[i-1][j], ga_algn->M[i-1][j] + ga_params->Co) + ga_params->Ci;
            ga_algn->D[i][j] = MIN2(ga_algn->D[i][j-1], ga_algn->M[i][j-1] + ga_params->Co) + ga_params->Cd;
            ga_res->cells += 1;
        }
    }

    ga_res->score = MIN3(ga_algn->M[ga_algn->len_query][ga_algn->len_target], 
                         ga_algn->I[ga_algn->len_query][ga_algn->len_target], 
                         ga_algn->D[ga_algn->len_query][ga_algn->len_target]);
    int end_time = clock();
    ga_res->elapsed = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000.0;
    ga_res->memory = get_memory_usage() - ga_res->memory;

    if (DEBUG == 1) {
        printf(" Normal Gap-Affine:         %5d\n", ga_res->score);
    }
    // Free allocated memory
    for (int i = 0; i <= ga_algn->len_query; i++) {
        free(ga_algn->M[i]);
        free(ga_algn->D[i]);
        free(ga_algn->I[i]);
    }
    free(ga_algn->M);
    free(ga_algn->D);
    free(ga_algn->I);
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