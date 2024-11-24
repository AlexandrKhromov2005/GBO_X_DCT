#include <iostream>
#include "stdio.h"
#include "image_processing.h"
#include "population.h"
#include "gbo.h"
#include "dct_functions.h"
#include <time.h>
#include <random>

double block_mse(unsigned char** original_block, unsigned char** changed_block) {
    double mse = 0.0;
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            mse += ((original_block[i][j] - changed_block[i][j]) * (original_block[i][j] - changed_block[i][j]));
        }
    }

    return mse / 64.0;
}

double block_psnr(unsigned char** original_block, unsigned char** changed_block) {
    return 10 * log10((255.0 * 255.0) / block_mse(original_block, changed_block));
}

int main() {
	srand(time(NULL));

	// Измерение времени выполнения
	clock_t start = clock();

	double** block = new double* [8];
	double** block_new = new double* [8];

	unsigned char** block_pic = new unsigned char* [8];
	unsigned char** block_new_pic = new unsigned char* [8];

	for (int i = 0; i < 8; ++i) {
		block[i] = new double[8];
		block_new[i] = new double[8];
		block_pic[i] = new unsigned char[8];
		block_new_pic[i] = new unsigned char[8];
	}

	// Инициализация значениями из матрицы
	block[0][0] = 180.0;  // Заменяем DC на 180.0
	block[0][1] = 3.87;
	block[0][2] = 3.65;
	block[0][3] = 0.13;
	block[0][4] = 0.38;
	block[0][5] = -1.87;
	block[0][6] = -3.27;
	block[0][7] = 5.24;

	block[1][0] = 5.75;
	block[1][1] = -1.09;
	block[1][2] = -0.20;
	block[1][3] = -2.72;
	block[1][4] = 0.26;
	block[1][5] = 0.46;
	block[1][6] = -3.80;
	block[1][7] = 2.44;

	block[2][0] = -4.83;
	block[2][1] = 0.18;
	block[2][2] = -1.33;
	block[2][3] = 0.17;
	block[2][4] = -0.28;
	block[2][5] = 0.38;
	block[2][6] = 3.21;
	block[2][7] = -3.11;

	block[3][0] = 2.04;
	block[3][1] = 2.16;
	block[3][2] = 2.71;
	block[3][3] = 2.27;
	block[3][4] = -2.62;
	block[3][5] = -1.19;
	block[3][6] = 2.27;
	block[3][7] = 0.33;

	block[4][0] = -1.13;
	block[4][1] = -1.91;
	block[4][2] = -1.37;
	block[4][3] = -2.32;
	block[4][4] = 1.63;
	block[4][5] = 1.29;
	block[4][6] = -0.38;
	block[4][7] = -0.49;

	block[5][0] = 1.80;
	block[5][1] = -0.23;
	block[5][2] = -1.00;
	block[5][3] = -0.81;
	block[5][4] = -1.87;
	block[5][5] = 1.75;
	block[5][6] = 0.39;
	block[5][7] = -0.47;

	block[6][0] = -2.19;
	block[6][1] = -0.69;
	block[6][2] = 1.21;
	block[6][3] = 1.28;
	block[6][4] = 1.61;
	block[6][5] = -2.98;
	block[6][6] = -2.92;
	block[6][7] = 0.74;

	block[7][0] = -1.01;
	block[7][1] = -2.18;
	block[7][2] = -1.40;
	block[7][3] = -1.33;
	block[7][4] = -1.13;
	block[7][5] = 2.68;
	block[7][6] = -0.18;
	block[7][7] = 1.57;

	double pr = 0.5;
	int th = 5; //разброс
	int n = 30; //кол-во векторов
	int m = 40; //кол-во итераций

	double** population = create_population(th, n);
	int best_ind, worst_ind = -1;

	for (int cur_iter = 0; cur_iter < m; cur_iter++) {

		printf("\nITERATION: %d\n", cur_iter);

		best_ind = find_x_best(population, block, n, 0);
		printf("\nBEST\n");
		for (int i = 0; i < 22; i++) {
			printf("%f ", population[best_ind][i]);
		}
		printf("\n");

		worst_ind = find_x_worst(population, block, n, 0);
		printf("\nWORST\n");
		for (int i = 0; i < 22; i++) {
			printf("%f ", population[worst_ind][i]);
		}
		printf("\n");

		gbo(population, best_ind, worst_ind, cur_iter, m, n, pr, th);
		}
		
	best_ind = find_x_best(population, block, n, 0);
	
	block_new = apply_x(block, population[best_ind]);

	double s1 = get_s1_sum(block_new);
	double s0 = get_s0_sum(block_new);

	printf("NEW S0 = %f , S1 = %f\n", s0, s1);

	rev_dct_func(block_pic, block);
	rev_dct_func(block_new_pic, block_new);

	double psnr = block_psnr(block_pic, block_new_pic);
	printf("\n FINAL PSNR = %f\n", psnr);
	
	// Завершение измерения времени
	clock_t end = clock();
	double time_spent = (double)(end - start) / CLOCKS_PER_SEC;

	printf("\nTime =  %f seconds\n", time_spent);

	// Освобождение памяти
	for (int i = 0; i < 8; ++i) {
		delete[] block[i];
		delete[] block_new[i];
		delete[] block_pic[i];
		delete[] block_new_pic[i];
	}
	delete[] block;
	delete[] block_new;
	delete[] block_pic;
	delete[] block_new_pic;

}
