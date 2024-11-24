#include <stdlib.h>
#include <time.h>
#include "stdio.h"
#include <math.h>
#include "block_metrics.h"
#include "dct_functions.h"
#include <float.h>
double get_s1_new_sum(double** block, double* x) {
    double result = 0.0;
    result += fabs(block[6][0] + x[0]);
    result += fabs(block[4][2] + x[2]);
    result += fabs(block[2][4] + x[4]);
    result += fabs(block[0][6] + x[6]);
    result += fabs(block[1][6] + x[8]);
    result += fabs(block[3][4] + x[10]);
    result += fabs(block[5][2] + x[12]);
    result += fabs(block[7][0] + x[14]);
    result += fabs(block[6][2] + x[16]);
    result += fabs(block[4][4] + x[18]);
    result += fabs(block[2][6] + x[20]);
    return result;
}

double get_s0_new_sum(double** block, double* x) {
    double result = 0.0;
    result += fabs(block[5][1] + x[1]);
    result += fabs(block[3][3] + x[3]);
    result += fabs(block[1][5] + x[5]);
    result += fabs(block[0][7] + x[7]);
    result += fabs(block[2][5] + x[9]);
    result += fabs(block[4][3] + x[11]);
    result += fabs(block[6][1] + x[13]);
    result += fabs(block[7][1] + x[15]);
    result += fabs(block[5][3] + x[17]);
    result += fabs(block[3][5] + x[19]);
    result += fabs(block[1][7] + x[21]);
    return result;
}

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

    double** population = (double**)malloc(d * sizeof(double*));

    for (int i = 0; i < d; i++) {
        population[i] = (double*)malloc(22 * sizeof(double));

        for (int j = 0; j < 22; j++) {
            population[i][j] = (double)-th + 2 * (double)th * rand_double();
        }
    }

    return population;
}

double of(double** original_dct_block, double** changed_dct_block, double* vector, const char b) {
    if (original_dct_block == NULL || changed_dct_block == NULL) {
        return -1.0; // Обработка ошибки
    }

    unsigned char** changed_block = (unsigned char**)malloc(8 * sizeof(unsigned char*));
    unsigned char** original_block = (unsigned char**)malloc(8 * sizeof(unsigned char*));

    unsigned char** new_dct_block = (double**)malloc(8 * sizeof(double*));


    if (changed_block == NULL || original_block == NULL) {
        return -1.0; // Обработка ошибки выделения памяти
    }

    for (int i = 0; i < 8; i++) {
        changed_block[i] = (unsigned char*)malloc(8 * sizeof(unsigned char));
        original_block[i] = (unsigned char*)malloc(8 * sizeof(unsigned char));
        new_dct_block[i] = (double*)malloc(8 * sizeof(double));


        if (changed_block[i] == NULL || original_block[i] == NULL) {
            return -1.0; // Обработка ошибки выделения памяти
        }
    }

    rev_dct_func(changed_block, changed_dct_block);
    rev_dct_func(original_block, original_dct_block);

    dct_func(changed_block, new_dct_block);

    double s1 = get_s1_sum(new_dct_block);
    double s0 = get_s0_sum(new_dct_block);
    double psnr = block_psnr(original_block, changed_block);

    printf("PSNR = %f", psnr);

    if ((b == 0) && (s0 == 0)) {
        s0 = 0.01;
    }
    else if ((b == 1) && (s1 == 0)) {
        s1 = 0.01;
    }

    double result = (b == 0) ? (s1 / s0) - (0.01 * psnr) : (s0 / s1) - (0.01 * psnr);

    printf(" F = %f\n", result);

    for (int i = 0; i < 8; i++) {
        free(changed_block[i]);
        free(original_block[i]);
        free(new_dct_block[i]);
    }
    free(changed_block);
    free(original_block);
    free(new_dct_block);

    return result;
}


double** apply_x(double** block, double* x) {
    if (block == NULL || x == NULL) {
        return NULL;
    }

    double** new_block = (double**)malloc(8 * sizeof(double*));
    if (new_block == NULL) {
        return NULL; // Обработка ошибки выделения памяти
    }

    for (int i = 0; i < 8; i++) {
        new_block[i] = (double*)malloc(8 * sizeof(double));
        if (new_block[i] == NULL) {
            // Освобождаем уже выделенную память
            for (int j = 0; j < i; j++) {
                free(new_block[j]);
            }
            free(new_block);
            return NULL; // Обработка ошибки выделения памяти
        }
    }



    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            new_block[i][j] = block[i][j];
        }
    }

    int index = 0;
    for (int i = 6; i >= 0; i--) {
        int row = i;
        int col = 6 - i;

        new_block[row][col] = sign(new_block[row][col]) * fabs(fabs(new_block[row][col]) + x[index]);
        index++;
    }
    for (int i = 0; i < 8; i++) {
        int row = i;
        int col = 7 - i;

        new_block[row][col] = sign(new_block[row][col]) * fabs(fabs(new_block[row][col]) + x[index]);
        index++;
    }
    for (int i = 7; i >= 1; i--) {
        int row = i;
        int col = 8 - i;

        new_block[row][col] = sign(new_block[row][col]) * fabs(fabs(new_block[row][col]) + x[index]);
        index++;
    }

    /*int index = 0;
    for (int i = 21; i <= 42; i++) {
        int row = i / 8;
        int col = i % 8;
        new_block[row][col] = sign(new_block[row][col]) * fabs(new_block[row][col] + x[index]);
        index++;
    }*/

    return new_block;
}

#include <limits.h>

int find_x_best(double** population, double** original_dct_block, int popul_size, const char b) {
    double best_value = DBL_MAX;
    int best_index = -1;

    for (int i = 0; i < popul_size; i++) {
        double** changed_dct_block = apply_x(original_dct_block, population[i]);
        double value = of(original_dct_block, changed_dct_block, population[i], b);

        if (value < best_value) {
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

int find_x_worst(double** population, double** original_dct_block, int popul_size, const char b) {
    double worst_value = -DBL_MAX;
    int worst_index = -1;

    for (int i = 0; i < popul_size; i++) {
        double** changed_dct_block = apply_x(original_dct_block, population[i]);
        double value = of(original_dct_block, changed_dct_block, population[i], b);

        if (value > worst_value) {
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

