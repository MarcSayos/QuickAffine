#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>

// Import the external C code files
#include "GapAffine.c"
#include "GapAffine_Windowed_BoundAndAlign.c"

#define NEG_INF INT_MIN
#define MAX_LEN 20000  // Define a maximum sequence length, adjust as needed


int get_memory_usage() {
    FILE* file = fopen("/proc/self/status", "r");
    char line[128];
    int memory_usage_kb = -1;  // Default to -1 in case of an error

    while (fgets(line, sizeof(line), file) != NULL) {
        if (strncmp(line, "VmRSS:", 6) == 0) {
            // Extract the memory usage value in kilobytes (KB)
            sscanf(line, "VmRSS: %d kB", &memory_usage_kb);
            break;
        }
    }
    fclose(file);
    return memory_usage_kb;
}


void costs_transform(int *a, int *x, int *o, int *i, int *d) {

    int alpha = 1, beta, gamma;

    //o' = alpha * o >= 0
    if (alpha * (*o) < 0) alpha = -alpha;

    // a' = alpha * a + beta + gamma = 0 -> gamma = -alpha * a - beta
    beta = (alpha * (*d - *a) > -alpha * (*i)) ? alpha * (*d - *a) : -alpha * (*i); 
    gamma = -alpha * (*a) - beta; 
    
    *a = alpha * (*a) + beta + gamma;  // Match cost should be 0
    *x = alpha * (*x) + beta + gamma; 
    *o = alpha * (*o); 
    *i = alpha * (*i) + beta; 
    *d = alpha * (*d) + gamma; 
}


int main(int argc, char *argv[]) {
    if (argc != 11) {
        fprintf(stderr, "Usage: %s <input_file> <output_file> <bases> <ws> <os> <Cm> <Cx> <Co> <Ci> <Cd>\n", argv[0]);
        return 1;
    }

    char *input_file = argv[1];
    char *output_file = argv[2];
    const int bases = atoi(argv[3]);
    const int ws = atoi(argv[4]);
    const int os = atoi(argv[5]);
    int Cm = atoi(argv[6]);
    int Cx = atoi(argv[7]);
    int Co = atoi(argv[8]);
    int Ci = atoi(argv[9]);
    int Cd = atoi(argv[10]);

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
    int score = 0,   memory = 0;
    int score_2 = 0, memory_2 = 0;
    double elapsed = 0, elapsed_2 = 0;
    clock_t start_time, end_time;

    double total_elapsed = 0, total_elapsed_2 = 0;

    costs_transform(&Cm, &Cx, &Co, &Ci, &Cd);

    while (fscanf(infile, ">%s\n<%s\n", query, target) == 2) {

        GapAffine(query, target, Cm, Cx, Co, Ci, Cd, &score, &memory, &elapsed);

        GapAffine_windowed(query, target, &ws, &os, &Cm, &Cx, &Co, &Ci, &Cd, &score_2, &memory_2, &elapsed_2);

        // Write the results to the output file
        fprintf(outfile, "Gap-Affine real score:       %d    %.4f ms     %d KB\n", score, elapsed, memory);
        fprintf(outfile, "Windowed Gap-Affine score:   %d    %.4f ms     %d KB\n", score_2, elapsed_2, memory_2);
        fprintf(outfile, "-----------------------------------------------\n");
        total_elapsed += elapsed;
        total_elapsed_2 += elapsed_2;
    }
    fprintf(outfile, "Total elapsed Gap-Affine:              %.4f ms\n", total_elapsed);
    fprintf(outfile, "Total elapsed Gap-Affine_Windowed:     %.4f ms\n", total_elapsed_2);
    // fprintf(outfile, "-----------------------------------------------\n");

    fclose(infile);
    fclose(outfile);

    return 0;
}
