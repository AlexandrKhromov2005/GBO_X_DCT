#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "population.h"

#define POPSIZE 100 

#define PI acos(-1)

//Вычисление rand от 0 до 1
double rand_num() {
	return (double)rand() / (double)RAND_MAX;
}

//Вычисление  нормально распределенного числа от 0 до 1
double randn() {
    static int hasSpare = 0;
    static double spare;
    if (hasSpare) {
        hasSpare = 0;
        return spare;
    }
    double u, v, s;
    do {
        u = (rand() / ((double)RAND_MAX)) * 2.0 - 1.0;
        v = (rand() / ((double)RAND_MAX)) * 2.0 - 1.0;
        s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);
    s = sqrt(-2.0 * log(s) / s);
    spare = v * s;
    hasSpare = 1;
    double result = u * s;
    result = (result + 3.0) / 6.0; 
    if (result < 0.0) result = 0.0;
    if (result > 1.0) result = 1.0;

    return result;
}

//Вычисление rho
double new_rho(double alpha) {
	return (2 * rand_num() * alpha) - alpha;
}

//Вычисление инедксов рандомных векторов
void gen_indexes(int indexes[5], int n, int cur_ind) {
	int cnt = 0;

	indexes[0] = -1;
	indexes[1] = -1;
	indexes[2] = -1;
	indexes[3] = -1;
	indexes[4] = -1;

	while (cnt < 5) {
		int temp = rand() % n;
		if ((indexes[0] != temp) && (indexes[1] != temp) && (indexes[2] != temp) && (indexes[3] != temp) && (indexes[4] != temp) && (cur_ind != temp)) {
			indexes[cnt] = temp;
			cnt++;
		}
	}
}

//Проверка x_next по значению OF
void check_x_next_for_of(double population[][22], double x_next[22], double f_x_next, xind f_list, int cur_vec) {
	if ((f_list.f_values[cur_vec] < f_x_next) && (f_list.f_values[f_list.worst] > f_x_next)) {
		return;
	}
	
	if (f_list.f_values[f_list.best] > f_x_next) {
		f_list.best = cur_vec;
	}	else if (f_list.f_values[f_list.worst] < f_x_next) {
				f_list.worst = cur_vec;
		}

	memcpy(population[cur_vec], x_next, 22);
	f_list.f_values[cur_vec] = f_x_next;
}

//Вычисление нового вектора через GSR
void gsr(double x_next[22], double x[22], double x1[22], double x2[22], double x_best[22], double x_worst[22], double xr1[22], double xr2[22], double xr3[22], double xr4[22], int cur_iter, int m, int n, double* alpha) {
	
	//Вычисление beta
	const double beta_min = 0.2;
	const double beta_max = 1.2;
	const double epsilon = 0.05;
	double beta = beta_min + (beta_max - beta_min) * pow(1.0 - pow(((double)cur_iter / (double)m), 3), 2);

	//Вычисление alpha
	*alpha = fabs(beta * sin(((3.0 * PI) / 2.0) + sin((3 * PI * beta) / 2.0)));


	double rho1 = new_rho(*alpha); //Инициализация rho1
	double rho2 = new_rho(*alpha); //Инициализация rho2
	double ra = rand_num();		   //Инициализация ra	
	double rb = rand_num();        //Инициализация rb	
	
	double x3[22];

	double rand_st;
	double randn_st;

	//Вычисление delta
	double delta[22];
	rand_st = rand_num();
	for (int i = 0; i < 22; i++) {
		delta[i] = 2.0 * rand_st * fabs(((xr1[i] + xr2[i] + xr3[i] + xr4[i]) / 4.0) - x[i]);
	}

	//Вычисление  delta X
	double x_delta[22];
	rand_st = ((double)(rand() % n) + 1.0);
	for (int i = 0; i < 22; i++) {
		x_delta[i] = rand_st * fabs(((x_best[i] - xr1[i]) + delta[i]) / 2.0);
	}

	//Вычисление Z
	double z [22];
	randn_st = randn();
	for (int i = 0; i < 22; i++) {
		z[i] = x[i] - randn_st * ((2.0 * x_delta[i] * x[i]) / (x_worst[i] - x_best[i] + epsilon));
	}

	//Вычисление P
	double p [22];
	rand_st = rand_num();
	randn_st = rand_num();
	for (int i = 0; i < 22; i++) {
		p[i] = rand_st * (((z[i] + x[i]) / 2.0) + randn_st * x_delta[i]);
	}

	//Вычисление Q
	double q [22];
	rand_st = rand_num();
	randn_st = rand_num();
	for (int i = 0; i < 22; i++) {
		q[i] = rand_st * (((z[i] + x[i]) / 2.0) - randn_st * x_delta[i]);
	}
	
	//Вычисление X1
	randn_st = randn();
	rand_st = rand_num();
	for (int i = 0; i < 22; i++) {
		x1[i] = x[i] - randn_st * rho1 * ((2 * x_delta[i] * x[i]) / (p[i] - q[i] + epsilon)) + rand_st * rho2 * (x_best[i] - x[i]);
	}
	
	//Вычисление X2
	randn_st = randn();
	rand_st = rand_num();
	for (int i = 0; i < 22; i++) {
		x2[i] = x_best[i] - randn_st * rho1 * ((2 * x_delta[i] * x[i]) / (p[i] - q[i] + epsilon)) + rand_st * rho2 * (xr1[i] - xr2[i]);
	}

	//Вычисление X3
	for (int i = 0; i < 22; i++) {
		x3[i] = x[i] - rho1 * (x2[i] - x1[i]);

		x_next[i] = ra * (rb * x1[i] + (1.0 - rb) * x2[i]) + (1.0 - ra) * x3[i];
	}

	//double rand_st = rand_num();
	////Вычисление вектора следующей популяции
	//for (int i = 0; i < 22; i++) {
	//	
	//	delta[i] = 2.0 * rand_st * fabs(((xr1[i] + xr2[i] + xr3[i] + xr4[i]) / 4.0) - x[i]);

	//	step[i] = ((x_best[i] - xr1[i]) + delta[i]) / 2.0;

	//	rand_st = rand_num();
	//	x_delta[i] = ((double)(rand() % n) + 1) * fabs(step[i]);

	//	z[i] = x[i] - randn() * ((2 * x_delta[i] * x[i]) / (x_worst[i] - x_best[i] + epsilon));

	//	p[i] = rand_num() * (((z[i] + x[i]) / 2.0) + rand_num() * x_delta[i]);

	//	q[i] = rand_num() * (((z[i] + x[i]) / 2.0) - rand_num() * x_delta[i]);

	//	x1[i] = x[i] - randn() * rho1 * ((2 * x_delta[i] * x[i]) / (p[i] - q[i] + epsilon)) + rand_num() * rho2 * (x_best[i] - x[i]);

	//	x2[i] = x_best[i] - randn() * rho1 * ((2 * x_delta[i] * x[i]) / (p[i] - q[i] + epsilon)) + rand_num() * rho2 * (xr1[i] - xr2[i]);

	//	x3[i] = x[i] - rho1 * (x2[i] - x1[i]);

	//	x_next[i] = ra * (rb * x1[i] + (1 - rb) * x2[i]) + (1 - ra) * x3[i];
	//}

}

//Вычисление нового вектора через LEO
void leo(double x_next[22], double x[22], double x_best[22], double x1[22], double x2[22], double xr1[22], double xr2[22], double x_p[22], double x_rand[22], double alpha, int th) {
	
	//Вычисление Y
	double y [22];
	if (rand_num() < 0.5) {
		memcpy(y, x, 22);
	}
	else{
		memcpy(y, x_best, 22);
	}

	//Инициализация mu1
	double mu1 = rand_num();

	//Вычисление u1, u2, u3
	double u1, u2, u3;
	if (mu1 < 0.5) {
		u1 = 2 * rand_num();
		u2 = rand_num();
		u3 = rand_num();
	}
	else {
		u1 = 1.0;
		u2 = 1.0;
		u3 = 1.0;
	}

	//Инициализация mu2
	double mu2 = rand_num();

	//Вычисление Xk
	double x_k [22];
	if (mu2 < 0.5) {
		memcpy(x_k, x_rand, 22);
	}
	else {
		memcpy(x_p, x_rand, 22);
	}

	double f1 = -1.0 + 2.0 * rand_num(); //Инициализация f1
	double f2 = -1.0 + 2.0 * rand_num(); //Инициализация f2

	double rho1 = new_rho(alpha); //Инициализация rho1

	//Вычисление вектора следующей популяции
	for (int i = 0; i < 22; i++) {
		x_next[i] = y[i] + f1 * (u1 * x_best[i] - u2 * x_k[i]) + f2 * rho1 * (u3 * (x2[i] - x1[i]) + u2 * (xr1[i] - xr2[i])) / 2.0;
	}
}

//Обновление популяции через GBO
void gbo(double population[][22], int m, int n, double pr, double th, double dct_block[8][8], unsigned char wm_bit) {
	//обработка популяции
	for (int cur_iter = 0; cur_iter < m; cur_iter++) {
		xind f_list = find_x_bw(population, dct_block, n, wm_bit);

		for (int cur_vec = 0; cur_vec < n; cur_vec++) {
			double alpha;		//инициализация alpha
			double x1[22];	   //вектор x1	
			double x2[22];	   //вектор x2
			double xr1[22];	   //1 рандомный вектор 	
			double xr2[22]; //2 рандомный вектор 	
			double xr3[22];    //3 рандомный вектор 	
			double xr4[22];    //4 рандомный вектор 	
			double x_p[22];    //рандомный вектор p	
			double x_rand[22]; //рандомно сгенерированный вектор
			for (int i = 0; i < 22; i++) {
				x_rand[i] = -th + 2.0 * th * rand_num();
			}
			double x_next[22];  //новый вектор

			int indexes[5];
			gen_indexes(indexes, n, cur_vec); //Массив индексов рандомных векторов

			//Копирование  рандомных векторов популяции
			memcpy(xr1, population[indexes[0]], 22);
			memcpy(xr2, population[indexes[1]], 22);
			memcpy(xr3, population[indexes[2]], 22);
			memcpy(xr4, population[indexes[3]], 22);
			memcpy(x_p, population[indexes[4]], 22);

			//Вычисление нового вектора через GSR
			gsr(x_next, population[cur_vec], x1, x2, population[f_list.best], population[f_list.worst], xr1, xr2, xr3, xr4, cur_iter, m, n, &alpha);

			if (rand_num() < pr) {
				//Вычисление нового вектора через LEO
				leo(x_next, x_next, population[best_ind], x1, x2, xr1, xr2, x_p, x_rand, alpha, th);
			}
		
			//Возвращаем значения вектора в передел [-Th, Th] 
			for (int i = 0; i < 22; i++) {
				if (x_next[i] > th) {
					x_next[i] = th;
				}
				else if (x_next[i] < -th) {
					x_next[i] = -th;
				}
			}

			//Расчитываем целевую функцию для нового вектора
			double changed_dct_block[8][8];
			apply_x(dct_block, x_next, changed_dct_block);
			double f_x_next = of(dct_block, changed_dct_block, x_next, wm_bit);

			check_x_next_for_of(population, x_next, f_x_next, f_list, cur_vec)
			
			}			
		}
	}
}
