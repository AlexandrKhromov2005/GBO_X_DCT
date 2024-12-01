#include <stdlib.h>
#include <time.h>
#include "stdio.h"
#include <math.h>
#include <float.h>
#include <limits.h>

#include "block_metrics.h"
#include "dct_functions.h"

/*typedef struct {
    int best;
    int worst;
}xind;*/

double get_s0_sum(double block[8][8]) {
    return fabs(block[5][1]) + fabs(block[6][1]) + fabs(block[7][1]) +
        fabs(block[3][3]) + fabs(block[4][3]) + fabs(block[5][3]) +
        fabs(block[1][5]) + fabs(block[2][5]) + fabs(block[3][5]) +
        fabs(block[0][7]) + fabs(block[1][7]);
}

/*#define GET_S0_SUM(block) (fabs((block)[5][1]) + fabs((block)[6][1]) + fabs((block)[7][1]) + \
                           fabs((block)[3][3]) + fabs((block)[4][3]) + fabs((block)[5][3]) + \
                           fabs((block)[1][5]) + fabs((block)[2][5]) + fabs((block)[3][5]) + \
                           fabs((block)[0][7]) + fabs((block)[1][7]))*/

double get_s1_sum(double block[8][8]) {
    return fabs(block[6][0]) + fabs(block[7][0]) + fabs(block[4][2]) +
        fabs(block[5][2]) + fabs(block[6][2]) + fabs(block[2][4]) +
        fabs(block[3][4]) + fabs(block[4][4]) + fabs(block[0][6]) +
        fabs(block[1][6]) + fabs(block[2][6]);
}

/*#define GET_S1_SUM(block) (fabs((block)[6][0]) + fabs((block)[7][0]) + fabs((block)[4][2]) + \
                           fabs((block)[5][2]) + fabs((block)[6][2]) + fabs((block)[2][4]) + \
                           fabs((block)[3][4]) + fabs((block)[4][4]) + fabs((block)[0][6]) + \
                           fabs((block)[1][6]) + fabs((block)[2][6]))*/

double rand_double() {
    return (double)rand() / RAND_MAX;
}

double sign(double n) {
    return (n > 0.0) ? 1.0 : -1.0;
}

void create_population(int th, int d, double population[][22]) {
    for (int i = 0; i < d; i++) {
        for (int j = 0; j < 22; j++) {
            population[i][j] = (double)-th + 2 * (double)th * rand_double();
        }
    }
}

double of(double original_dct_block[8][8], double changed_dct_block[8][8], double vector[22], const char b) {
    unsigned char changed_block[8][8];
    unsigned char original_block[8][8];
    double new_dct_block[8][8];

    rev_dct_func(changed_block, changed_dct_block);
    rev_dct_func(original_block, original_dct_block);

    dct_func(changed_block, new_dct_block);

    double s1 = get_s1_sum(new_dct_block);
    double s0 = get_s0_sum(new_dct_block);
    double psnr = block_psnr(original_block, changed_block);

    //printf("PSNR = %f", psnr);

    if ((b == 0) && (s0 == 0)) {
        s0 = 0.1;
    }
    else if ((b == 1) && (s1 == 0)) {
        s1 = 0.1;
    }

    double result = (b == 0) ? (s1 / s0) - (0.01 * psnr) : (s0 / s1) - (0.01 * psnr);

    //printf(" F = %f\n", result);

    return result;
}


void apply_x(double block[8][8], double x[22], double new_block[8][8]) {
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
}



int find_x_best(double population[][22], double original_dct_block[8][8], int popul_size, const char b) {
    double best_value = DBL_MAX;
    int best_index = -1;

    for (int i = 0; i < popul_size; i++) {
        double changed_dct_block[8][8];
        apply_x(original_dct_block, population[i], changed_dct_block);
        double value = of(original_dct_block, changed_dct_block, population[i], b);

        if (value < best_value) {
            best_value = value;
            best_index = i;
        }
    }
    return best_index;
}

int find_x_worst(double population[][22], double original_dct_block[8][8], int popul_size, const char b) {
    double worst_value = -DBL_MAX;
    int worst_index = -1;

    for (int i = 0; i < popul_size; i++) {
        double changed_dct_block[8][8];
        apply_x(original_dct_block, population[i], changed_dct_block);
        double value = of(original_dct_block, changed_dct_block, population[i], b);

        if (value > worst_value) {
            worst_value = value;
            worst_index = i;
        }
    }

    return worst_index;
}


/*xind find_x_bw(double population[][22], double original_dct_block[8][8], int popul_size, const char b){
    xind result;
    double worst_value = -DBL_MAX;
    double best_value = DBL_MAX;
    result.best = -1;
    result.worst = -1;
    for (int i = 0; i < popul_size; i++){
        double changed_dct_block[8][8];
        apply_x(original_dct_block, population[i], changed_dct_block);
        double value = of(original_dct_block, changed_dct_block, population[i], b);

        if (value > worst_value) {
            result.worst = i;
            worst_value = value;
        }

        if (value < best_value) {
            result.best = i;
            best_value = value;
        }
    }
    return result;
} */

