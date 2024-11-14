#include "main.h"
#include "GapAffine.h"
#include "GapAffine_Windowed_BoundAndAlign.h"

// Function to create a matrix with initial values
uint16_t **create_matrix(int rows, int cols) {
    uint16_t **matrix = (uint16_t **)malloc(rows * sizeof(uint16_t *));
    if (!matrix) {
        perror("Failed to allocate matrix");
        return NULL;  // Return NULL on error
    }
    for (int i = 0; i < rows; i++) {
        matrix[i] = malloc(cols * sizeof(uint16_t));
        if (!matrix[i]) {
            perror("Failed to allocate row");
            return NULL;  // Return NULL on error
        }
        memset(matrix[i], 0xFFFF, cols * sizeof(uint16_t)); // Fill with __UINT8_MAX__ by setting each byte to 0xFF
    }
    return matrix;  // Return the allocated and initialized matrix
}

// Function to create a matrix with initial values
// uint16_t **create_matrix(int rows, int cols) {
//     uint16_t **matrix = (uint16_t **)malloc(rows * sizeof(uint16_t *));
//     for (uint16_t i = 0; i < rows; i++) {
//         matrix[i] = (uint16_t *)malloc(cols * sizeof(uint16_t));
//     }
//     return matrix;
// }

void reset_matrices(GapAffine_Alignment *ga_algn) {

    for (int i = 0; i < ga_algn->len_query; i++) {
        memset(ga_algn->M[i], 0xFFFF, ga_algn->len_target * sizeof(uint16_t)); // Fill with __UINT8_MAX__ by setting each byte to 0xFF
        memset(ga_algn->I[i], 0xFFFF, ga_algn->len_target * sizeof(uint16_t)); // Fill with __UINT8_MAX__ by setting each byte to 0xFF
        memset(ga_algn->D[i], 0xFFFF, ga_algn->len_target * sizeof(uint16_t)); // Fill with __UINT8_MAX__ by setting each byte to 0xFF
    }
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
    // *beta = (*alpha * (ga_params->Cd - ga_params->Cm) > -*alpha * (ga_params->Ci)) ? *alpha * (ga_params->Cd - ga_params->Cm) : -*alpha * (ga_params->Ci); 
    // *beta = (*alpha * (ga_params->Ci) + 1 != 0 && *alpha * (ga_params->Cd) + (-*alpha * (ga_params->Cm) - 1) != 0) ? 1 : 2;
    *beta = ((ga_params->Cd) - (ga_params->Cm) == 0 && 1 < *alpha * (ga_params->Ci)) ? (-1) : (1);
    *gamma = -*alpha * (ga_params->Cm) - *beta; 
    
    ga_params->Cm = *alpha * (ga_params->Cm) + *beta + *gamma;  // Match cost should be 0
    ga_params->Cx = *alpha * (ga_params->Cx) + *beta + *gamma; 
    ga_params->Co = *alpha * (ga_params->Co); 
    ga_params->Ci = *alpha * (ga_params->Ci) + *beta; 
    ga_params->Cd = *alpha * (ga_params->Cd) + *gamma; 

}

void python_plot_print(GapAffine_Alignment *ga_algn) {
    FILE *file = fopen("matrix.csv", "w");
    if (file == NULL) {
        printf("Error opening file!\n");
    }

    for (int i = 0; i < ga_algn->len_query+1; i++) {
        for (int j = 0; j < ga_algn->len_target+1; j++) {
            if (ga_algn->M[i][j] == 65535) {fprintf(file, "%d,", 99);printf("%2d,", 99);}
            else {fprintf(file, "%2d,", ga_algn->M[i][j]); printf("%2d,", ga_algn->M[i][j]);}
        }
        fprintf(file, "\n");
        printf("\n");
    }

    fclose(file);
}


int main(int argc, char *argv[]) {
    if (DEBUG == 1) {
        if (argc != 11) {
            fprintf(stderr, "Usage:       %s <input_file> <output_file> <bases> <ws> <os> <Cm> <Cx> <Co> <Ci> <Cd>\n", argv[0]);
            fprintf(stderr, "Example run: %s test_datasets/test_dataset_1000b.seq res.out 1000 64 16 0 6 5 3 3\n", argv[0]);
            return 1;
        }
    } else {
        if (argc != 12) {
            fprintf(stderr, "Usage:       %s <input_file> <output_file> <bases> <ws> <os> <Cm> <Cx> <Co> <Ci> <Cd> <name>\n", argv[0]);
            return 1;
        }
    }

    char *input_file = argv[1];
    char *output_file = argv[2];
    char name[40];

    if (DEBUG != 1) {
        strncpy(name, argv[11], sizeof(name) - 1);  // Copy up to 39 characters
        name[39] = '\0';  // Null-terminate in case argv[11] is longer than 39 chars
    }

    GapAffine_Parameters ga_params = {
        atoi(argv[3]),
        atoi(argv[4]),
        atoi(argv[5]),
        atoi(argv[6]),
        atoi(argv[7]),
        atoi(argv[8]),
        atoi(argv[9]),
        atoi(argv[10]),
        atoi(argv[6]),
        atoi(argv[7]),
        atoi(argv[8]),
        atoi(argv[9]),
        atoi(argv[10])
    };


    GapAffine_Alignment ga_algn;
    ga_algn.query = (char *)malloc((ga_params.bases + 1) * sizeof(char));
    ga_algn.target = (char *)malloc((ga_params.bases + 1) * sizeof(char));
    
    double total_elapsed_1 = 0, total_elapsed_2 = 0;
    int total_memory_1 = 0, total_memory_2 = 0; 
    int total_cells_1 = 0, total_cells_2 = 0;
    int total_cells_banded = 0, total_cells_windowed = 0;
    costs_transform(&ga_params);

    FILE *infile = fopen(input_file, "r");
    if (infile == NULL) {
        perror("Error opening input file");
        return 1;
    }

    FILE *outfile;
    if (DEBUG == 1) outfile = fopen(output_file, "w");
    else outfile = fopen(output_file, "a");
    if (outfile == NULL) {
        perror("Error opening output file");
        fclose(infile);
        return 1;
    }
    
    while (fscanf(infile, ">%s\n<%s\n", ga_algn.query, ga_algn.target) == 2) {
        ga_algn.len_query = (int) strlen(ga_algn.query);
        ga_algn.len_target = (int) strlen(ga_algn.target);
        ga_algn.cigar_len = 0;
        GapAffine_Results ga_res_1 = {0,0,0,0,0,0,0,0};
        GapAffine_Results ga_res_2 = {0,0,0,0,0,0,0,0};
        
        GapAffine(&ga_algn, &ga_params, &ga_res_1);

        GapAffine_windowed(&ga_algn, &ga_params, &ga_res_2);

        // Write the results to the output file
        if (DEBUG == 1) {
            printf("Original Penalties: {%d,%d,%d,%d,%d}, New Penalties: {%d,%d,%d,%d,%d}\n", ga_params.or_Cm, ga_params.or_Cx, ga_params.or_Co, ga_params.or_Ci, ga_params.or_Cd, ga_params.Cm, ga_params.Cx, ga_params.Co, ga_params.Ci, ga_params.Cd);
            fprintf(outfile, "Algorithm            |     Score | Elapsed(ms) | Memory(KB) | Computed Cells | Computed Cells Banded\n");
            // fprintf(outfile, "Total cells: %63d\n", (int) (strlen(ga_algn.query)*strlen(ga_algn.target)));
            fprintf(outfile, "Gap-Affine real score:    %6d | %11.4f | %10d | %14d |\n",     ga_res_1.score, ga_res_1.elapsed, ga_res_1.memory, ga_res_1.computed_cells_score);
            fprintf(outfile, "Windowed Gap-Affine bound:%6d\n", ga_res_2.bound);
            fprintf(outfile, "Windowed Gap-Affine score:%6d | %11.4f | %10d | %14d | %d \n", ga_res_2.original_score, ga_res_2.elapsed, ga_res_2.memory, ga_res_2.computed_cells_windowed, ga_res_2.computed_cells_banded);
            // fprintf(outfile, "Original penalties  score:   %d\n", ga_res_2.original_score);
            fprintf(outfile, "-----------------------------------------------------------\n");
        } 
        if (DEBUG == 2) {
            fprintf(outfile, "RealScore = %d, Bound = %d, BandedScore = %d\n", ga_res_1.score, ga_res_2.bound, ga_res_2.original_score);
        } else if (DEBUG == 3) {
            python_plot_print(&ga_algn); 
        }

        total_elapsed_1 += ga_res_1.elapsed;
        total_elapsed_2 += ga_res_2.elapsed;
        total_memory_1 += ga_res_1.memory;
        total_memory_2 += ga_res_2.memory;
        total_cells_1 += ga_res_1.computed_cells_score;
        total_cells_2 += (ga_res_2.computed_cells_windowed+ga_res_2.computed_cells_banded);
        total_cells_windowed += ga_res_2.computed_cells_windowed;
        total_cells_banded += ga_res_2.computed_cells_banded;

        // Free matrices
        for (int i = 0; i <= ga_algn.len_query; i++) {
            free(ga_algn.M[i]);
            free(ga_algn.I[i]);
            free(ga_algn.D[i]);
        }
        free(ga_algn.M);
        free(ga_algn.I);
        free(ga_algn.D);
    }

    if (DEBUG == 1) {
        fprintf(outfile, "Total elapsed Gap-Affine:              %.4f ms\n", total_elapsed_1);
        fprintf(outfile, "Total elapsed Gap-Affine_Windowed:     %.4f ms\n", total_elapsed_2);
    } else {
        fprintf(outfile, "%s, Penalties {%d,%d,%d,%d,%d}, ", name, ga_params.or_Cm, ga_params.or_Cx, ga_params.or_Co, ga_params.or_Ci, ga_params.or_Cd);
        fprintf(outfile, "Window Size: %d, Overlap Size: %d\n", ga_params.ws, ga_params.os);
        fprintf(outfile, "Algorithm    | Elapsed(ms) | Memory(KB) | Computed Cells | Windowed |   Banded\n");
        fprintf(outfile, "SWG:            %10.4f | %10d | %14d |          |\n", total_elapsed_1, total_memory_1, total_cells_1);
        fprintf(outfile, "QuickedAffine:  %10.4f | %10d | %14d | %8d | %8d\n", total_elapsed_2, total_memory_2, total_cells_2, total_cells_windowed, total_cells_banded);
        fprintf(outfile, "------------------------------------------------------------------------------\n");
    }

    free(ga_algn.cigar);
    free(ga_algn.query);
    free(ga_algn.target);

    fclose(infile);
    fclose(outfile);

    return 0;
}
