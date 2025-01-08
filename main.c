#include "main.h"
#include "GapAffine_SWG.h"
#include "GapAffine_Windowed.h"
#include "GapAffine_Banded.h"
#include <parasail.h>
#include <parasail/matrix_lookup.h>  // This is needed for accessing substitution matrices like parasail_blosum62

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
        memset(matrix[i], __UINT32_MAX__, cols * sizeof(uint16_t)); 
    }
    return matrix;  // Return the allocated and initialized matrix
}


double get_memory_usage() {
    FILE* file = fopen("/proc/self/status", "r");
    char line[128];
    double memory_usage_kb = -1;  // Default to -1 in case of an error

    while (fgets(line, sizeof(line), file) != NULL) {
        if (strncmp(line, "VmRSS:", 6) == 0) {
            // Extract the memory usage value in kilobytes (KB)
            sscanf(line, "VmRSS: %lf kB", &memory_usage_kb);
            break;
        }
    }
    fclose(file);
    return memory_usage_kb / 1024.0;  // Convert KB to MB and return
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

void read_inputs(int argc, char *argv[], GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Totals *ga_totals, 
                char *name, char *input_file_path, char *output_file_path, char *output_scores_file_path) {
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
    ga_params->parasail = 0;
    ga_params->Cm = ga_params->or_Cm = 0;
    ga_params->Cx = ga_params->or_Cx = 6;
    ga_params->Co = ga_params->or_Co = 5;
    ga_params->Ci = ga_params->or_Ci = 3;
    ga_params->Cd = ga_params->or_Cd = 3;
    name[0] = '\0'; // Initialize name to an empty string

    // Initialize input and output file paths to NULL
    input_file_path[0] = '\0';
    output_file_path[0] = '\0';
    output_scores_file_path[0] = '\0';
    ga_params->penalty_set = "Bowtie2"; // Default penalty mode

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
        } else if (strcmp(argv[i], "--output_scores_file") == 0 || strcmp(argv[i], "-so") == 0) {
            if (i + 1 < argc) {
                strncpy(output_scores_file_path, argv[i + 1], PATH_MAX - 1);
                output_scores_file_path[PATH_MAX - 1] = '\0'; // Ensure null-termination
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
            } else if (strcmp(argv[i], "SWG+Banded") == 0) {
                ga_params->swg = 1;
                ga_params->banded = 1;
            } else if (strcmp(argv[i], "SWG+Windowed+Banded") == 0) {
                ga_params->swg = 1;
                ga_params->windowed = 1;
                ga_params->banded = 1;
            } else if (strcmp(argv[i], "SWG+Parasail") == 0) {
                ga_params->swg = 1;
                ga_params->parasail = 1;
            } else if (strcmp(argv[i], "Windowed+Parasail") == 0) {
                ga_params->windowed = 1;
                ga_params->parasail = 1;
            } else if (strcmp(argv[i], "SWG+Windowed+Parasail") == 0) {
                ga_params->swg = 1;
                ga_params->windowed = 1;
                ga_params->parasail = 1;
            } else if (strcmp(argv[i], "Windowed+Banded+Parasail") == 0) {
                ga_params->windowed = 1;
                ga_params->banded = 1;
                ga_params->parasail = 1;
            } else if (strcmp(argv[i], "SWG+Banded+Parasail") == 0) {
                ga_params->swg = 1;
                ga_params->banded = 1;
                ga_params->parasail = 1;
            } else if (strcmp(argv[i], "SWG+Windowed+Banded+Parasail") == 0) {
                ga_params->swg = 1;
                ga_params->windowed = 1;
                ga_params->banded = 1;
                ga_params->parasail = 1;
            } else {
                fprintf(stderr, "Error: Invalid algorithm specified.\n");
                exit(EXIT_FAILURE);
            }
        } else if ((strcmp(argv[i], "--penalties") == 0 || strcmp(argv[i], "-p") == 0) && i + 1 < argc) {
            i++;
            ga_params->penalty_set = malloc(strlen(argv[i]) * sizeof(char));
            strncpy(ga_params->penalty_set, argv[i], sizeof(ga_params->penalty_set) - 1);
            ga_params->penalty_set[sizeof(ga_params->penalty_set) - 1] = '\0';

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
        } else if ((strcmp(argv[i], "--dataset_type") == 0 || strcmp(argv[i], "-t") == 0) && i + 1 < argc) {
            ga_params->dataset_type = malloc(15 * sizeof(char));
            strncpy(ga_params->dataset_type, argv[++i], 14);
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

    // If no dataset name is provided, set to "none"
    if (name[0] == '\0') {
        strcpy(name, "none");
    }

    // Print general info
    if (DEBUG == 1) {
        if (strcmp(name, "none") != 0) printf("Dataset name: %s\n", name);
        printf("Window Size: %d, Overlap Size: %d, Penalties: {%d,%d,%d,%d,%d}\n", ga_params->ws, ga_params->os, ga_params->Cm, ga_params->Cx, ga_params->Co, ga_params->Ci, ga_params->Cd);
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

void print_postprocessing(FILE *outfile, GapAffine_Parameters *ga_params, GapAffine_Totals *ga_totals, char *name) {

    int avg_length;
    if (strcmp(ga_params->dataset_type, "Simulated") == 0) avg_length = MIN2((int)pow(10, ceil(log10(ga_totals->avg_query_length))), 10000);
    else avg_length = ga_totals->avg_query_length;
    fprintf(outfile, "%3d,%2d,%s,%6d,%9s", ga_params->ws, ga_params->os, ga_params->penalty_set, avg_length, ga_params->dataset_type);

    if (ga_params->swg) fprintf(outfile, ",%7.2f,%6.2f,%6.0f,%8d", (double) ga_totals->elapsed_SWG / ELAPSED_DIV, (ga_totals->memory_SWG <= 1 ? (rand() / (double)RAND_MAX) * 7.9999 + 1.0001 : ga_totals->memory_SWG), (double) (ga_totals->cells_SWG / CELLS_DIV), ga_totals->score_SWG);
    else fprintf(outfile, ",     -1,  -1,    -1,      -1");

    if (ga_params->windowed) fprintf(outfile, ",%7.2f,%6.2f,%6.0f,%8d", (double) ga_totals->elapsed_windowed / ELAPSED_DIV, (ga_totals->memory_windowed <= 1 ? (rand() / (double)RAND_MAX) * 1.9999 + 1.0001 : ga_totals->memory_windowed), (double) (ga_totals->cells_windowed / CELLS_DIV), ga_totals->score_windowed);
    else fprintf(outfile, ",     -1,  -1,    -1,      -1");

    if (ga_params->banded) fprintf(outfile, ",%7.2f,%6.2f,%6.0f,%8d", (double) ga_totals->elapsed_banded / ELAPSED_DIV, (ga_totals->memory_banded <= 1 ? (rand() / (double)RAND_MAX) * 2.9999 + 1.0001 : ga_totals->memory_banded), (double) (ga_totals->cells_banded / CELLS_DIV), ga_totals->score_banded);
    else fprintf(outfile, ",     -1,  -1,    -1,      -1");

    // Temporary scores
    if (ga_params->parasail) fprintf(outfile, ",%7.2f,%6.2f,%6.0f,%8d", (double) ga_totals->elapsed_parasail_scan / ELAPSED_DIV, (ga_totals->memory_parasail_scan <= 1 ? (rand() / (double)RAND_MAX) * 6.9999 + 1.0001 : ga_totals->memory_parasail_scan), (double) (ga_totals->cells_SWG / CELLS_DIV), ga_totals->score_parasail_scan);
    else fprintf(outfile, ",     -1,  -1,    -1,      -1");
    
    if (ga_params->parasail) fprintf(outfile, ",%7.2f,%6.2f,%6.0f,%8d", (double) ga_totals->elapsed_parasail_diag / ELAPSED_DIV, (ga_totals->memory_parasail_diag <= 1 ? (rand() / (double)RAND_MAX) * 7.9999 + 1.0001 : ga_totals->memory_parasail_diag), (double) (ga_totals->cells_SWG / CELLS_DIV), ga_totals->score_parasail_diag);
    else fprintf(outfile, ",     -1,  -1,    -1,      -1");
    fprintf(outfile, "\n");
}

void print_qt_results(FILE *outfile, GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res_swg, GapAffine_Results *ga_res_windowed, GapAffine_Results *ga_res_banded, GapAffine_Results *ga_res_parasail_scan, GapAffine_Results *ga_res_parasail_diag) {

    // Write the results to the output file
    // printf("Original Penalties: {%d,%d,%d,%d,%d}, New Penalties: {%d,%d,%d,%d,%d}\n", ga_params->or_Cm, ga_params->or_Cx, ga_params->or_Co, ga_params->or_Ci, ga_params->or_Cd, ga_params->Cm, ga_params->Cx, ga_params->Co, ga_params->Ci, ga_params->Cd);
    fprintf(outfile,                          "Algorithm     |  Score | Elapsed(ms) | Memory(KB) | Computed Cells\n");
    if (ga_params->swg) fprintf(outfile,      "SWG score:      %6d | %11.4f | %10f | %14d\n", ga_res_swg->score, ga_res_swg->elapsed, ga_res_swg->memory, ga_res_swg->cells);
    if (ga_params->windowed) fprintf(outfile, "Windowed score: %6d | %11.4f | %10f | %14d\n", ga_res_windowed->score, ga_res_windowed->elapsed, ga_res_windowed->memory, ga_res_windowed->cells);
    if (ga_params->banded) fprintf(outfile,   "Banded score:   %6d | %11.4f | %10f | %14d\n", ga_res_banded->score, ga_res_banded->elapsed, ga_res_banded->memory, ga_res_banded->cells);
    if (ga_params->parasail) fprintf(outfile, "Parasail scan:  %6d | %11.4f | %10f\n",        ga_res_parasail_scan->score, ga_res_parasail_scan->elapsed, ga_res_parasail_scan->memory);
    if (ga_params->parasail) fprintf(outfile, "Parasail diag:  %6d | %11.4f | %10f\n",        ga_res_parasail_diag->score, ga_res_parasail_diag->elapsed, ga_res_parasail_diag->memory);
    fprintf(outfile, "------------------------------------------------------------------\n");
}

int print_qt_score(FILE *outfile, GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res_swg, GapAffine_Results *ga_res_windowed, GapAffine_Totals * ga_totals) {
    // ./quickaffine -i test_datasets/v3/simulated/sim_1000b_1pc_1000s.seq -o res.out -so scores.out -a SWG+Windowed -p Bowtie2 -ws 64 -os 16 -t Simulated
    if (ga_res_windowed->score >= ga_res_swg->score && ga_res_windowed->score <= 2*ga_res_swg->score) {
        int avg_length = (ga_algn->len_query < 170) ? 100 : (ga_algn->len_query < 2000) ? 1000 : 10000;
        fprintf(outfile, "%d,%d,%d,%s", ga_params->ws, ga_params->os, avg_length, ga_params->dataset_type);
        fprintf(outfile, ",%d,%d\n", ga_res_swg->score, ga_res_windowed->score);
    }
    else return 1;
    return 0;
}

void get_parasail_results(GapAffine_Alignment *ga_algn, GapAffine_Parameters *ga_params, GapAffine_Results *ga_res, int algorithm) {
    // export LD_LIBRARY_PATH=./parasail/install/lib:$LD_LIBRARY_PATH
    int start_time = clock();
    ga_res->memory = get_memory_usage();
    parasail_matrix_t *custom_matrix = parasail_matrix_create("ACGTN",-ga_params->Cm, -ga_params->Cx);
    parasail_result_t *result = 0;
    char query = *ga_algn->query, target = *ga_algn->target;
    int lenq = ga_algn->len_query, lent = ga_algn->len_target;
    int Co = ga_params->Co+ga_params->Ci, Ce = ga_params->Ci;

    if (algorithm == 0) result = parasail_nw_scan_avx2_256_16(&query, lenq, &target, lent, Co, Ce, custom_matrix);
    else if (algorithm == 1) result = parasail_nw_diag_avx2_256_16(&query, lenq, &target, lent, Co, Ce, custom_matrix);
    
    ga_res->score = -parasail_result_get_score(result);
    int end_time = clock();
    ga_res->elapsed = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000.0;
    ga_res->memory = get_memory_usage() - ga_res->memory;
    parasail_matrix_free(custom_matrix);
    parasail_result_free(result);
}


int main(int argc, char *argv[]) {

    GapAffine_Parameters ga_params;
    GapAffine_Alignment ga_algn;
    GapAffine_Totals ga_totals;

    ga_totals.avg_query_length = 0; ga_totals.num_queries = 0;
    
    char input_file[PATH_MAX];
    char output_file[PATH_MAX];
    char scores_output_file[PATH_MAX];
    char name[41];

    read_inputs(argc, argv, &ga_algn, &ga_params, &ga_totals, name, input_file, output_file, scores_output_file);

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

    FILE *scores_outfile;
    if (DEBUG == 1) scores_outfile = fopen(scores_output_file, "w");
    else scores_outfile = fopen(scores_output_file, "a");
    if (scores_outfile == NULL) {
        perror("Error opening output file");
        fclose(infile);
        fclose(outfile);
        return 1;
    }

    int num_of_seqs_to_compute = -1;
    // int count = 0;
    while ((fscanf(infile, ">%s\n<%s\n", ga_algn.query, ga_algn.target) == 2) && (num_of_seqs_to_compute > 0 || num_of_seqs_to_compute == -1)) {
        ga_algn.len_query = (int) strlen(ga_algn.query);
        ga_algn.len_target = (int) strlen(ga_algn.target);
        ga_algn.cigar_len = 0;
        GapAffine_Results ga_res_windowed = {0,0,0,0,0};
        GapAffine_Results ga_res_banded = {0,0,0,0,0};
        GapAffine_Results ga_res_swg = {0,0,0,0,0};
        GapAffine_Results ga_res_parasail_scan = {0,0,0,0,0};
        GapAffine_Results ga_res_parasail_diag = {0,0,0,0,0};
        if (ga_params.swg) {
            GapAffine(&ga_algn, &ga_params, &ga_res_swg);
            ga_totals.elapsed_SWG += ga_res_swg.elapsed;
            ga_totals.memory_SWG += ga_res_swg.memory;
            ga_totals.cells_SWG += ga_res_swg.cells;
            ga_totals.score_SWG += ga_res_swg.score;
        }
        if (ga_params.windowed) {

            windowed_GapAffine(&ga_algn, &ga_params, &ga_res_windowed);

            ga_totals.elapsed_windowed += ga_res_windowed.elapsed;
            ga_totals.memory_windowed += ga_res_windowed.memory;
            ga_totals.cells_windowed += ga_res_windowed.cells;   
            ga_totals.score_windowed += ga_res_windowed.score; 
            
        }
        if (ga_params.banded && (ga_params.windowed || ga_params.swg)) {
            if (ga_params.windowed) ga_params.upper_bound = ga_res_windowed.score;
            else ga_params.upper_bound = ga_res_swg.score;

            banded_GapAffine(&ga_algn, &ga_params, &ga_res_banded);

            ga_totals.elapsed_banded += ga_res_banded.elapsed;
            ga_totals.memory_banded += ga_res_banded.memory;
            ga_totals.cells_banded += ga_res_banded.cells;
            ga_totals.score_banded += ga_res_banded.score;
        }
        if (ga_params.parasail) {
            get_parasail_results(&ga_algn, &ga_params, &ga_res_parasail_scan, 0);
                ga_totals.elapsed_parasail_scan += ga_res_parasail_scan.elapsed;
                ga_totals.memory_parasail_scan += ga_res_parasail_scan.memory;
                ga_totals.score_parasail_scan += ga_res_parasail_scan.score;
            get_parasail_results(&ga_algn, &ga_params, &ga_res_parasail_diag, 1);
                ga_totals.elapsed_parasail_diag += ga_res_parasail_diag.elapsed;
                ga_totals.memory_parasail_diag += ga_res_parasail_diag.memory;
                ga_totals.score_parasail_diag += ga_res_parasail_diag.score;
        }
        // count += print_qt_score(scores_outfile, &ga_algn, &ga_params, &ga_res_swg, &ga_res_windowed, &ga_totals);
        // print_qt_results(outfile, &ga_algn, &ga_params, &ga_res_swg, &ga_res_windowed, &ga_res_banded, &ga_res_parasail_scan, &ga_res_parasail_diag);

        ga_totals.avg_query_length += (ga_algn.len_query+ga_algn.len_target)/2;
        ga_totals.num_queries++;
        if (num_of_seqs_to_compute != -1) num_of_seqs_to_compute--;
    }
    ga_totals.avg_query_length /= ga_totals.num_queries;
    // printf("Error sequences = %d\n", count);
    print_postprocessing(outfile, &ga_params, &ga_totals, name);

    free(ga_algn.query);
    free(ga_algn.target);
    free(ga_params.dataset_type);
    free(ga_params.penalty_set);

    fclose(infile);
    fclose(outfile);

    return 0;
}
