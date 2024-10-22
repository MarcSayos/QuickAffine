#include "main.h"
// Import the external C code files
#include "GapAffine.c"
#include "GapAffine_Windowed_BoundAndAlign.c"


// Function to create a matrix with initial values
__uint32_t **create_matrix(int rows, int cols) {
    __uint32_t **matrix = (__uint32_t **)malloc(rows * sizeof(__uint32_t *));
    for (__uint32_t i = 0; i < rows; i++) {
        matrix[i] = (__uint32_t *)malloc(cols * sizeof(__uint32_t));
    }
    return matrix;
}

void print_memory_usage() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    
    // Memory usage (in kilobytes)
    printf("Memory usage: %ld KB\n", usage.ru_maxrss);
}

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


void costs_transform(GapAffine_Parameters *ga_params) {
    int *alpha = &ga_params->alpha, *beta = &ga_params->beta, *gamma = &ga_params->gamma;

    if (ga_params->Cm >= 0 && ga_params->Cx >= 0 && ga_params->Co >= 0 && ga_params->Ci >= 0 && ga_params->Cd >= 0) {
        *alpha = 0, *beta = 0, *gamma = 0;
        return;
    }

    *alpha = 1;
    //o' = alpha * o >= 0
    if (*alpha * (ga_params->Co) < 0) *alpha = -*alpha;

    // a' = alpha * a + beta + gamma = 0 -> gamma = -alpha * a - beta
    *beta = (*alpha * (ga_params->Cd - ga_params->Cm) > -*alpha * (ga_params->Ci)) ? *alpha * (ga_params->Cd - ga_params->Cm) : -*alpha * (ga_params->Ci); 
    *gamma = -*alpha * (ga_params->Cm) - *beta; 
    
    ga_params->Cm = *alpha * (ga_params->Cm) + *beta + *gamma;  // Match cost should be 0
    ga_params->Cx = *alpha * (ga_params->Cx) + *beta + *gamma; 
    ga_params->Co = *alpha * (ga_params->Co); 
    ga_params->Ci = *alpha * (ga_params->Ci) + *beta; 
    ga_params->Cd = *alpha * (ga_params->Cd) + *gamma; 

}



int main(int argc, char *argv[]) {
    if (argc != 11) {
        fprintf(stderr, "Usage: %s <input_file> <output_file> <bases> <ws> <os> <Cm> <Cx> <Co> <Ci> <Cd>\n", argv[0]);
        return 1;
    }

    char *input_file = argv[1];
    char *output_file = argv[2];

    GapAffine_Parameters ga_params = {
        atoi(argv[3]),
        atoi(argv[4]),
        atoi(argv[5]),
        atoi(argv[6]),
        atoi(argv[7]),
        atoi(argv[8]),
        atoi(argv[9]),
        atoi(argv[10])
    };

    GapAffine_Results ga_res_1 = {0,0,0,0};
    GapAffine_Results ga_res_2 = {0,0,0,0};

    GapAffine_Alignment ga_algn;

    char query[ga_params.bases+ga_params.bases/10], target[ga_params.bases+ga_params.bases/10];
    double total_elapsed = 0, total_elapsed_2 = 0;
    costs_transform(&ga_params);

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
    
    while (fscanf(infile, ">%s\n<%s\n", query, target) == 2) {

        GapAffine(query, target, &ga_params, &ga_res_1);

        GapAffine_windowed(query, target, &ga_params, &ga_res_2);

        // Write the results to the output file
        fprintf(outfile, "Query: %s\nTarget: %s\nCm=%d, Cx=%d, Co=%d, Ci=%d, Cd=%d\n", query, target, ga_params.Cm, ga_params.Cx, ga_params.Co, ga_params.Ci, ga_params.Cd);
        fprintf(outfile, "Gap-Affine real score:       %d    %.4f ms     %d KB\n", ga_res_1.score, ga_res_1.elapsed, ga_res_1.memory);
        fprintf(outfile, "Windowed Gap-Affine score:   %d    %.4f ms     %d KB\n", ga_res_2.score, ga_res_2.elapsed, ga_res_2.memory);
        fprintf(outfile, "-----------------------------------------------\n");

        total_elapsed += ga_res_1.elapsed;
        total_elapsed_2 += ga_res_2.elapsed;
    }

    // free(ga_algn.query);
    // free(ga_algn.target);
    // free(ga_algn.cigar);
    // free(ga_algn.M);
    // free(ga_algn.I);
    // free(ga_algn.D);

    fprintf(outfile, "Total elapsed Gap-Affine:              %.4f ms\n", total_elapsed);
    fprintf(outfile, "Total elapsed Gap-Affine_Windowed:     %.4f ms\n", total_elapsed_2);

    fclose(infile);
    fclose(outfile);

    return 0;
}
