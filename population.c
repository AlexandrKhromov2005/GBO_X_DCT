#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "block_metrics.h"
#include "dct.h"
#include <float.h>

double get_s0_sum(double** block) {
    return fabs(block[5][1]) + fabs(block[6][1]) + fabs(block[7][1]) +
        fabs(block[3][3]) + fabs(block[4][3]) + fabs(block[5][3]) +
        fabs(block[1][5]) + fabs(block[2][5]) + fabs(block[3][5]) +
        fabs(block[0][7]) + fabs(block[1][7]);
}

double get_s1_sum(double** block) {
    return fabs(block[6][0]) + fabs(block[7][0]) + fabs(block[4][2]) +
        fabs(block[5][2]) + fabs(block[6][2]) + fabs(block[2][4]) +
        fabs(block[3][4]) + fabs(block[4][4]) + fabs(block[0][6]) +
        fabs(block[1][6]) + fabs(block[2][6]);
}

double rand_double() {
    return (double)rand() / RAND_MAX;
}

double sign(double n) {
    return (n > 0.0) ? 1.0 : -1.0;
}

double** create_population(int th, int d) {
    srand(time(NULL));

    double** population = (double**)malloc(d * sizeof(double*));

    for (int i = 0; i < d; i++) {
        population[i] = (double*)malloc(22 * sizeof(double));

        for (int j = 0; j < 22; j++) {
            population[i][j] = (double)-th + 2 * (double)th * rand_double();
        }
    }

    return population;
}

double of(double** original_dct_block, double** changed_dct_block, char b) {
    unsigned char** changed_block = (unsigned char**)malloc(8 * sizeof(unsigned char*));
    unsigned char** original_block = (unsigned char**)malloc(8 * sizeof(unsigned char*));
    
    rev_dct_func(changed_block, changed_dct_block);
    rev_dct_func(original_block, original_dct_block);


    return (b == 0) ? (get_s1_sum(changed_block) / get_s0_sum(changed_block)) - 0.01 * block_psnr(original_block, changed_block) : (get_s0_sum(changed_block) / get_s1_sum(changed_block)) - 0.01 * block_psnr(original_block, changed_block);
}

double** apply_x(double** block, double* x) {
    double** new_block = (double**)malloc(8 * sizeof(double*));

    for (int i = 0; i < 8; i++) {
        memcpy(new_block[i], block[i], 8 * sizeof(double));
    }

    int index = 0;
    for (int i = 21; i <= 42; i++) {
        int row = i / 8;
        int col = i % 8;
        new_block[row][col] = sign(new_block[row][col]) * fabs(new_block[row][col] + x[index]);
        index++;
    }

    return new_block;
}

#include <limits.h>

int find_x_best(double** population, double** original_dct_block, char b, int d) {
    double best_value = -DBL_MAX;
    int best_index = -1;

    for (int i = 0; i < d; i++) {
        double** changed_dct_block = apply_x(original_dct_block, population[i]);
        double value = of(original_dct_block, changed_dct_block, b);

        if (value > best_value) {
            best_value = value;
            best_index = i;
        }

        for (int j = 0; j < 8; j++) {
            free(changed_dct_block[j]);
        }
        free(changed_dct_block);
    }

    return best_index;
}

int find_x_worst(double** population, double** original_dct_block, char b, int d) {
    double worst_value = DBL_MAX;
    int worst_index = -1;

    for (int i = 0; i < d; i++) {
        double** changed_dct_block = apply_x(original_dct_block, population[i]);
        double value = of(original_dct_block, changed_dct_block, b);

        if (value < worst_value) {
            worst_value = value;
            worst_index = i;
        }

        for (int j = 0; j < 8; j++) {
            free(changed_dct_block[j]);
        }
        free(changed_dct_block);
    }

    return worst_index;
}

