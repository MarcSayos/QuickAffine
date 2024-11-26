#include "main.h"
#include "GapAffine_SWG.h"
#include "GapAffine_Windowed.h"
#include "GapAffine_Banded.h"

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

void reset_matrices(GapAffine_Alignment *ga_algn, int previous_length, int rows, int cols) {
    
    // Free matrices
    for (int i = 0; i <= previous_length; i++) {
        free(ga_algn->M[i]);
        free(ga_algn->I[i]);
        free(ga_algn->D[i]);
    }
    free(ga_algn->M);
    free(ga_algn->I);
    free(ga_algn->D);

    ga_algn->M = create_matrix(cols, rows);
    ga_algn->I = create_matrix(cols, rows);
    ga_algn->D = create_matrix(cols, rows);
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

    *beta = ((ga_params->Cd) - (ga_params->Cm) == 0 && 1 < *alpha * (ga_params->Ci)) ? (-1) : (1);
    *gamma = -*alpha * (ga_params->Cm) - *beta; 
    
    ga_params->Cm = *alpha * (ga_params->Cm) + *beta + *gamma;  // Match cost should be 0
    ga_params->Cx = *alpha * (ga_params->Cx) + *beta + *gamma; 
    ga_params->Co = *alpha * (ga_params->Co); 
    ga_params->Ci = *alpha * (ga_params->Ci) + *beta; 
    ga_params->Cd = *alpha * (ga_params->Cd) + *gamma; 

}


void read_inputs(int argc, char *argv[], GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Totals *ga_totals, char *name, char *input_file_path, char *output_file_path) {
    
    // Check if the number of arguments is valid
    if (argc < 2) {
        fprintf(stderr, "Usage:       %s \n    --input_file | -i <input_file> \n    --output_file | -o <output_file> \n    --algorithm | -a <SWG|Windowed|Windowed+Banded|SWG+Windowed+Banded> \n    --penalties | -p <Bowtie2|BWA-MEM|Personalized> [<Cm> <Cx> <Co> <Ci> <Cd>] \n    --dataset_name | -n <name> \n    --window_size | -ws <ws> \n    --overlap_size | -os <os>\n", argv[0]);
        fprintf(stderr, "Example run: %s ./quickaffine -i test_datasets/real/Nanopore_100s.seq -o res.out -a SWG+Windowed+Banded\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    // Default values
    ga_params->ws = 64;
    ga_params->os = 16;
    ga_params->swg = 0;
    ga_params->windowed = 0;
    ga_params->banded = 0;
    ga_params->Cm = ga_params->or_Cm = 0;
    ga_params->Cx = ga_params->or_Cx = 6;
    ga_params->Co = ga_params->or_Co = 5;
    ga_params->Ci = ga_params->or_Ci = 3;
    ga_params->Cd = ga_params->or_Cd = 3;
    name[0] = '\0'; // Initialize name to an empty string

    // Initialize input and output file paths to NULL
    input_file_path[0] = '\0';
    output_file_path[0] = '\0';
    char penalties_mode[20] = "default"; // Default penalty mode

    // Iterate over the input arguments
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--input_file") == 0 || strcmp(argv[i], "-i") == 0) {
            if (i + 1 < argc) {
                strncpy(input_file_path, argv[i + 1], PATH_MAX - 1);
                input_file_path[PATH_MAX - 1] = '\0'; // Ensure null-termination
                i++; // Skip the next argument since we've used it
            } else {
                fprintf(stderr, "Error: Missing argument for input_file\n");
                exit(EXIT_FAILURE);
            }
        } else if (strcmp(argv[i], "--output_file") == 0 || strcmp(argv[i], "-o") == 0) {
            if (i + 1 < argc) {
                strncpy(output_file_path, argv[i + 1], PATH_MAX - 1);
                output_file_path[PATH_MAX - 1] = '\0'; // Ensure null-termination
                i++; // Skip the next argument since we've used it
            } else {
                fprintf(stderr, "Error: Missing argument for output_file\n");
                exit(EXIT_FAILURE);
            }
        }  else if ((strcmp(argv[i], "--algorithm") == 0 || strcmp(argv[i], "-a") == 0) && i + 1 < argc) {
            i++;
            if (strcmp(argv[i], "SWG") == 0) {
                ga_params->swg = 1;
            } else if (strcmp(argv[i], "Windowed") == 0) {
                ga_params->windowed = 1;
            } else if (strcmp(argv[i], "SWG+Windowed") == 0) {
                ga_params->swg = 1;
                ga_params->windowed = 1;
            } else if (strcmp(argv[i], "Windowed+Banded") == 0) {
                ga_params->windowed = 1;
                ga_params->banded = 1;
            } else if (strcmp(argv[i], "SWG+Windowed+Banded") == 0) {
                ga_params->swg = 1;
                ga_params->windowed = 1;
                ga_params->banded = 1;
            } else {
                fprintf(stderr, "Error: Invalid algorithm specified.\n");
                exit(EXIT_FAILURE);
            }
        } else if ((strcmp(argv[i], "--penalties") == 0 || strcmp(argv[i], "-p") == 0) && i + 1 < argc) {
            i++;
            strncpy(penalties_mode, argv[i], sizeof(penalties_mode) - 1);
            penalties_mode[sizeof(penalties_mode) - 1] = '\0';

            if (strcmp(argv[i], "Bowtie2") == 0) {
                ga_params->Cm = ga_params->or_Cm = 0;
                ga_params->Cx = ga_params->or_Cx = 6;
                ga_params->Co = ga_params->or_Co = 5;
                ga_params->Ci = ga_params->or_Ci = 3;
                ga_params->Cd = ga_params->or_Cd = 3;
            } else if (strcmp(argv[i], "BWA-MEM") == 0) {
                ga_params->Cm = ga_params->or_Cm = -1;
                ga_params->Cx = ga_params->or_Cx = 4;
                ga_params->Co = ga_params->or_Co = 6;
                ga_params->Ci = ga_params->or_Ci = 2;
                ga_params->Cd = ga_params->or_Cd = 2;
            } else if (strcmp(argv[i], "Personalized") == 0 && i + 5 < argc) {
                ga_params->Cm = ga_params->or_Cm = atoi(argv[++i]);
                ga_params->Cx = ga_params->or_Cx = atoi(argv[++i]);
                ga_params->Co = ga_params->or_Co = atoi(argv[++i]);
                ga_params->Ci = ga_params->or_Ci = atoi(argv[++i]);
                ga_params->Cd = ga_params->or_Cd = atoi(argv[++i]);
            } else {
                fprintf(stderr, "Error: Invalid penalties specified or insufficient values for Personalized mode.\n");
                exit(EXIT_FAILURE);
            }
        } else if ((strcmp(argv[i], "--dataset_name") == 0 || strcmp(argv[i], "-n") == 0) && i + 1 < argc) {
            strncpy(name, argv[++i], 40);
        } else if ((strcmp(argv[i], "--window_size") == 0 || strcmp(argv[i], "-ws") == 0) && i + 1 < argc) {
            ga_params->ws = atoi(argv[++i]);
        } else if ((strcmp(argv[i], "--overlap_size") == 0 || strcmp(argv[i], "-os") == 0) && i + 1 < argc) {
            ga_params->os = atoi(argv[++i]);
        } else {
            fprintf(stderr, "Error: Unrecognized argument '%s'.\n", argv[i]);
            exit(EXIT_FAILURE);
        }
    }

    // Mandatory checks
    if (input_file_path[0] == '\0') {
        fprintf(stderr, "Error: --input_file is mandatory.\n");
        exit(EXIT_FAILURE);
    }

    if (!ga_params->swg && !ga_params->windowed && !ga_params->banded) {
        fprintf(stderr, "Error: --algorithm is mandatory.\n");
        exit(EXIT_FAILURE);
    }

    // Debug mode if no output file
    if (output_file_path[0] == '\0') {
        printf("Debug mode enabled.\n");
    } else {
        printf("Output file: %s\n", output_file_path);
    }

    // If no dataset name is provided, set to "none"
    if (name[0] == '\0') {
        strcpy(name, "none");
    }

    // Allocate memory for query, target, and cigar
    ga_algn->query = malloc(MAX_LEN);  // Adjust size as needed
    ga_algn->target = malloc(MAX_LEN);
    ga_algn->cigar = malloc(MAX_LEN);
    
    if (!ga_algn->query || !ga_algn->target || !ga_algn->cigar) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
}

void print_qt_results(FILE *outfile, GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res_swg, GapAffine_Results *ga_res_windowed, GapAffine_Results *ga_res_banded) {

        // Write the results to the output file
        if (DEBUG == 1) {
            printf("Original Penalties: {%d,%d,%d,%d,%d}, New Penalties: {%d,%d,%d,%d,%d}\n", ga_params->or_Cm, ga_params->or_Cx, ga_params->or_Co, ga_params->or_Ci, ga_params->or_Cd, ga_params->Cm, ga_params->Cx, ga_params->Co, ga_params->Ci, ga_params->Cd);
            fprintf(outfile, "Algorithm            |     Score | Elapsed(ms) | Memory(KB) | Computed Cells | Computed Cells Banded\n");
            // fprintf(outfile, "Total cells: %63d\n", (int) (strlen(ga_algn.query)*strlen(ga_algn.target)));
            fprintf(outfile, "Gap-Affine real score:    %6d | %11.4f | %10d | %14d |\n",     ga_res_swg->score, ga_res_swg->elapsed, ga_res_swg->memory, ga_res_swg->cells);
            fprintf(outfile, "Windowed Gap-Affine bound:%6d\n", ga_res_windowed->score);
            fprintf(outfile, "Windowed Gap-Affine score:%6d | %11.4f | %10d | %14d | %d \n", ga_res_banded->original_score, ga_res_banded->elapsed, ga_res_banded->memory, ga_res_windowed->cells, ga_res_banded->cells);
            // fprintf(outfile, "Original penalties  score:   %d\n", ga_res_banded->original_score);
            fprintf(outfile, "-----------------------------------------------------------\n");
        } 
        if (DEBUG == 2) {
            fprintf(outfile, "RealScore = %d, Bound = %d, BandedScore = %d\n", ga_res_swg->score, ga_res_windowed->score, ga_res_banded->score);
        } else if (DEBUG == 3) {
            python_plot_print(ga_algn); 
        }
}

void print_total_results(FILE *outfile, GapAffine_Parameters *ga_params, GapAffine_Totals *ga_totals, char *name) {

    if (DEBUG == 1) {
        fprintf(outfile, "Total elapsed Gap-Affine:              %.4f ms\n", ga_totals->elapsed_SWG);
        fprintf(outfile, "Total elapsed Gap-Affine_Windowed:     %.4f ms\n", ga_totals->elapsed_windowed);
        fprintf(outfile, "Total elapsed Gap-Affine_Banded:       %.4f ms\n", ga_totals->elapsed_banded);
    } else {
        fprintf(outfile, "%s, Penalties {%d,%d,%d,%d,%d}, ", name, ga_params->or_Cm, ga_params->or_Cx, ga_params->or_Co, ga_params->or_Ci, ga_params->or_Cd);
        fprintf(outfile, "Window Size: %d, Overlap Size: %d\n", ga_params->ws, ga_params->os);
        fprintf(outfile, "Algorithm    | Elapsed(ms) | Memory(KB) | Computed Cells \n");
        fprintf(outfile, "SWG:            %10.4f | %10d | %14d\n", ga_totals->elapsed_SWG, ga_totals->memory_SWG, ga_totals->cells_SWG);
        fprintf(outfile, "Windowed:       %10.4f | %10d | %14d\n", ga_totals->elapsed_windowed, ga_totals->memory_windowed, ga_totals->cells_windowed);
        fprintf(outfile, "Banded:         %10.4f | %10d | %14d\n", ga_totals->elapsed_banded, ga_totals->memory_banded, ga_totals->cells_banded);
        fprintf(outfile, "------------------------------------------------------------------------------\n");
    }
}

int main(int argc, char *argv[]) {

    GapAffine_Parameters ga_params;
    GapAffine_Alignment ga_algn;
    GapAffine_Totals ga_totals;

    
    char input_file[PATH_MAX];
    char output_file[PATH_MAX];
    char name[41];

    read_inputs(argc, argv, &ga_algn, &ga_params, &ga_totals, name, input_file, output_file);

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
        GapAffine_Results ga_res_windowed = {0,0,0,0,0};
        GapAffine_Results ga_res_banded = {0,0,0,0,0};
        GapAffine_Results ga_res_swg = {0,0,0,0,0};
        if (ga_params.swg) {
            GapAffine(&ga_algn, &ga_params, &ga_res_swg);
            ga_totals.elapsed_SWG += ga_res_swg.elapsed;
            ga_totals.memory_SWG += ga_res_swg.memory;
            ga_totals.cells_SWG += ga_res_swg.cells;
        }
        if (ga_params.windowed) {

            windowed_GapAffine(&ga_algn, &ga_params, &ga_res_windowed);

            ga_totals.elapsed_windowed += ga_res_windowed.elapsed;
            ga_totals.memory_windowed += ga_res_windowed.memory;
            ga_totals.cells_windowed += ga_res_windowed.cells;    
            
            if (ga_params.banded) {
                ga_params.upper_bound = ga_res_windowed.score;

                banded_GapAffine(&ga_algn, &ga_params, &ga_res_banded);

                ga_totals.elapsed_banded += ga_res_banded.elapsed;
                ga_totals.memory_banded += ga_res_banded.memory;
                ga_totals.cells_banded += ga_res_banded.cells;
            }
        }
        // print_qt_results(outfile, &ga_algn, &ga_params, &ga_res_swg, &ga_res_windowed, &ga_res_banded);
    }
    // print_total_results(outfile, &ga_params, &ga_totals, name);

    free(ga_algn.query);
    free(ga_algn.target);

    fclose(infile);
    fclose(outfile);

    return 0;
}
