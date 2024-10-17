#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// Import the external C code files
#include "GapAffine.c"
#include "GapAffine_Windowed_BoundAndAlign.c"

#define NEG_INF INT_MIN
#define MAX_LEN 20000  // Define a maximum sequence length, adjust as needed
#define max(a,b) \
({ __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b; })

int main(int argc, char *argv[]) {
    if (argc != 10) {
        fprintf(stderr, "Usage: %s <input_file> <output_file> <bases> <ws> <os> <Cm> <Cx> <Co> <Ce>\n", argv[0]);
        return 1;
    }

    char *input_file = argv[1];
    char *output_file = argv[2];
    const int bases = atoi(argv[3]);
    const int ws = atoi(argv[4]);
    const int os = atoi(argv[5]);
    const int Cm = atoi(argv[6]);
    const int Cx = atoi(argv[7]);
    const int Co = atoi(argv[8]);
    const int Ce = atoi(argv[9]);

    FILE *infile = fopen(input_file, "r");
    if (infile == NULL) {
        perror("Error opening input file");
        return 1;
    }

    FILE *outfile = fopen(output_file, "w");
    if (outfile == NULL) {
        perror("Error opening output file");
        fclose(infile);
        return 1;
    }

    char query[bases+bases/10], target[bases+bases/10];
    int real_score, bound;
    int total_real_score = 0, total_bound = 0;
    clock_t start_time, end_time;
    double elapsed, elapsed_2;
    double total_elapsed = 0, total_elapsed_2 = 0;

    while (fscanf(infile, ">%s\n<%s\n", query, target) == 2) {
        // Measure time and compute real score
        start_time = clock();
        real_score = AffineGap(query, target, Cm, Cx, Co, Ce);
        end_time = clock();
        elapsed = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000.0;
        // printf("Gap-Affine real score:             %d\n", real_score);
        // Measure time and compute bound using windowed method
        start_time = clock();
        bound = main_windowed(query, target, &ws, &os, &Cm, &Cx, &Co, &Ce);
        end_time = clock();
        elapsed_2 = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000.0;

        total_real_score += real_score;
        total_bound += bound;
        total_elapsed += elapsed;
        total_elapsed_2 += elapsed_2;
        // Write the results to the output file
        fprintf(outfile, "Gap-Affine real score:       %d    %.4f ms\n", real_score, elapsed);
        fprintf(outfile, "Windowed Gap-Affine score:   %d    %.4f ms\n", bound, elapsed_2);
        fprintf(outfile, "-----------------------------------------------\n");
    }

    fclose(infile);
    fclose(outfile);

    return 0;
}
