#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>


void check_objectives(int a, int x, int o, int i, int d, int ai, int xi, int oi, int ii, int di, int alpha, int beta, int gamma) {
    if (ai != (alpha * (a) + beta + gamma)) {printf("a' is not transformed properly\n");return;}
    if (xi != (alpha * (x) + beta + gamma)) {printf("x' is not transformed properly\n");return;}
    if (oi != (alpha * (o)               )) {printf("o' is not transformed properly\n");return;}
    if (ii != (alpha * (i) + beta        )) {printf("i' is not transformed properly\n");return;}
    if (di != (alpha * (d)        + gamma)) {printf("d' is not transformed properly\n");return;}


    if (alpha * (x-a) < 0) {printf("Err1: input values are not correct\n");return;}
    if (alpha*(d-a) < -alpha * i) {printf("Err2: input values are not correct\n");return;}
    if (alpha > 0 && (d+i) < a) {printf("Err3: input values are not correct\n");return;}
    if (alpha < 0 && (d+i) > a) {printf("Err4: input values are not correct\n");return;}

    if (ai != 0) {printf("a' is not 0\n");return;}
    if (xi < 0)  {printf("x' is not > 0\n");return;}
    if (oi < 0)  {printf("o' is not > 0\n");return;}
    if (ii < 0)  {printf("i' is not > 0\n");return;}
    if (di < 0)  {printf("d' is not > 0\n");return;}
}

void costs_transform_5c(int a, int x, int o, int i, int d, int *ai, int *xi, int *oi, int *ii, int *di, int *alpha, int *beta, int *gamma) {
    *alpha = 1;

    //o' = alpha * o >= 0
    if (*alpha * o < 0) *alpha = -*alpha;

    // a' = alpha * a + beta + gamma = 0 -> gamma = -alpha * a - beta
    // *beta = (*alpha * (d - a) > -*alpha * (i)) ? *alpha * (d - a) : -*alpha * (i); 
    // *beta = (*alpha * i + 1 != 0 && *alpha * d + (-*alpha * a - 1) != 0) ? 1 : -1;
    *beta = (d - a == 0 && 1 < *alpha * i) ? (-1) : (1);
    *gamma = -*alpha * (a) - *beta; 
    
    *ai = *alpha * (a) + *beta + *gamma;  // Match cost should be 0
    *xi = *alpha * (x) + *beta + *gamma; 
    *oi = *alpha * (o); 
    *ii = *alpha * (i) + *beta; 
    *di = *alpha * (d) + *gamma; 

    check_objectives(a,x,o,i,d,*ai,*xi,*oi,*ii,*di, *alpha, *beta, *gamma);
}

void costs_transform_4c(int a, int x, int o, int e, int *ai, int *xi, int *oi, int *ei, int *alpha, int *beta) {

    *alpha = o < 0 ? -1 : 1;
    if (!((*alpha > 0 && x >= a && e >= a) || (*alpha < 0 && x <= a && e <= a))) printf("not possible\n");
    *beta = -*alpha * a;

    *ai = *alpha * (a) + *beta;
    *xi = *alpha * (x) + *beta; 
    *oi = *alpha * (o); 
    *ei = *alpha * (e) + *beta; 
}

// Compilar: gcc costs_transform.c -o costs_transform -lm
// Executar: ./costs_transform test_costs.seq
int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <input_file>\n", argv[0]);
        return 1;
    }
    char *input_file = argv[1];

    FILE *infile = fopen(input_file, "r");
    if (infile == NULL) {
        perror("Error opening input file");
        return 1;
    }

    int a,x,o,i,d;
    while (fscanf(infile, "{%d,%d,%d,%d,%d}\n", &a,&x,&o,&i,&d) == 5) {
        int ai = 0, xi = 0, oi = 0, ii = 0, di = 0;
        int alpha, beta, gamma;
        costs_transform_5c(a, x, o, i, d, &ai, &xi, &oi, &ii, &di, &alpha, &beta, &gamma);
        printf("{%2d,%2d,%2d,%2d,%2d} -> {%2d,%2d,%2d,%2d,%2d}, w/ alpha = %2d, beta = %2d, gamma = %2d.\n", a, x, o, i, d, ai, xi, oi, ii, di, alpha, beta, gamma);
        // costs_transform_4c(a, x, o, i, &ai, &xi, &oi, &ii, &alpha, &beta);
        // printf("{%2d,%2d,%2d,%2d} -> {%2d,%2d,%2d,%2d}, w/ alpha = %2d, beta = %2d\n", a, x, o, i, ai, xi, oi, ii, alpha, beta);
        printf("-----------------------------------------------\n");
    }
    
    return 0;
}