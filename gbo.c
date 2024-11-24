#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define PI acos(-1)

//Вычисление rand от 0 до 1
double rand_num() {
	return (double)rand() / (double)RAND_MAX;
}

//Вычисление  нормально распределенного числа от 0 до 1
double randn() {
	static int use_last = 0;
	static double y2;

	if (use_last) {
		use_last = 0;
		return (y2 + 3) / 6.0;
	}

	double x1, x2, w, y1;

	do {
		x1 = 2.0 * rand_num() - 1.0;
		x2 = 2.0 * rand_num() - 1.0;
		w = x1 * x1 + x2 * x2;
	} while (w >= 1.0 || w == 0.0);

	w = sqrt((-2.0 * log(w)) / w);
	y1 = x1 * w;
	y2 = x2 * w;

	use_last = 1;

	return (y1 + 3) / 6.0;
}

//Вычисление rho
double new_rho(double alpha) {
	return (2 * rand_num() * alpha) - alpha;
}

//Вычисление инедксов рандомных векторов
int* gen_indexes(int n, int cur_ind, int best_ind, int worst_ind) {
	int* result = (int*)calloc(5, sizeof(int));
	int cnt = 0;
	while (cnt < 5) {
		int temp = rand() % n;
		if ((result[0] != temp) && (result[1] != temp) && (result[2] != temp) && (result[3] != temp) && (result[4] != temp) && (cur_ind != temp)) {
			result[cnt] = temp;
			cnt++;
		}
	}
	return result;
}

//Вычисление нового вектора через GSR
double* gsr(double* x, double* x1, double* x2, double* x_best, double* x_worst, double* xr1, double* xr2, double* xr3, double* xr4, int cur_iter, int m, int n, double* alpha) {
	
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
	
	//Вычисление delta
	double* delta = (double*)calloc(22, sizeof(double));
	for (int i = 0; i < 22; i++) {
		delta[i] = 2.0 * rand_num() * fabs(((xr1[i] + xr2[i] + xr3[i] + xr4[i]) / 4.0) - x[i]);
	}

	//Вычисление step
	double* step = (double*)calloc(22, sizeof(double));
	for (int i = 0; i < 22; i++) {
		step[i] = ((x_best[i] - xr1[i]) + delta[i]) / 2.0;
	}

	//Вычиследние delta X
	double* x_delta = (double*)calloc(22, sizeof(double));
	for (int i = 0; i < 22; i++) {
		x_delta[i] = ((double)(rand() % n) + 1) * fabs(step[i]);
	}

	//Вычисление Z
	double* z = (double*)calloc(22, sizeof(double));
	for (int i = 0; i < 22; i++) {
		z[i] = x[i] - randn() * ((2 * x_delta[i] * x[i]) / (x_worst[i] - x_best[i] + epsilon));
	}

	//Вычисление P
	double* p = (double*)calloc(22, sizeof(double));
	for (int i = 0; i < 22; i++) {
		p[i] = rand_num() * (((z[i] + x[i]) / 2.0) + rand() * x_delta[i]);
	}

	//Вычисление Q
	double* q = (double*)calloc(22, sizeof(double));
	for (int i = 0; i < 22; i++) {
		q[i] = rand() * (((z[i] + x[i]) / 2.0) - rand_num() * x_delta[i]);
	}

	//Вычисление X1
	for (int i = 0; i < 22; i++) {
		x1[i] = x[i] - randn() * rho1 * ((2 * x_delta[i] * x[i]) / (p[i] - q[i] + epsilon)) + rand_num() * rho2 * (x_best[i] - x[i]);
	}

	//Вычисление X2
	for (int i = 0; i < 22; i++) {
		x2[i] = x_best[i] - randn() * rho1 * ((2 * x_delta[i] * x[i]) / (p[i] - q[i] + epsilon)) + rand_num() * rho2 * (xr1[i] - xr2[i]);
	}
	
	//Вычисление X3
	double* x3 = (double*)calloc(22, sizeof(double));
	for (int i = 0; i < 22; i++) {
		x3[i] = x[i] - rho1 * (x2[i] - x1[i]);
	}

	//Вычисление вектора следующей популяции
	double* x_next = (double*)calloc(22, sizeof(double));
	for (int i = 0; i < 22; i++) {
		x_next[i] = ra * (rb * x1[i] + (1 - rb) * x[i]) + (1 -ra) * x3[i];
	}

	//Отчистка памяти
	free(delta);
	free(step);
	free(x_delta);
	free(z);
	free(p);
	free(q);
	free(x3);

	//Возврат вектора
	return x_next;
}

//Вычисление нового вектора через LEO
double* leo(double* x, double* x_best, double* x1, double* x2, double* xr1, double* xr2, double* x_p, double* x_rand, double alpha, int th) {
	
	//Вычисление Y
	double* y = (double*)malloc(22 * sizeof(double));
	if (rand_num() < 0.5) {
		memcpy(y, x, 22 * sizeof(double));
	}
	else{
		memcpy(y, x_best, 22 * sizeof(double));
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
	double* x_k = (double*)malloc(22 * sizeof(double));
	if (mu2 < 0.5) {
		memcpy(x_k, x_rand, 22 * sizeof(double));
	}
	else {
		memcpy(x_p, x_rand, 22 * sizeof(double));
	}

	double f1 = -1.0 + 2.0 * rand_num(); //Инициализация f1
	double f2 = -1.0 + 2.0 * rand_num(); //Инициализация f2

	double rho1 = new_rho(alpha); //Инициализация rho1

	//Вычисление вектора следующей популяции
	double* x_next = (double*)malloc(22 * sizeof(double));
	for (int i = 0; i < 22; i++) {
		x_next[i] = y[i] + f1 * (u1 * x_best[i] - u2 * x_k[i]) + f2 * rho1 * (u3 * (x2[i] - x1[i]) + u2 * (xr1[i] - xr2[i])) / 2.0;
		if (x_next[i] > th) {
			x_next[i] = th;
		}
		else if (x_next[i] < -th) {
			x_next[i] = -th;
		}
	}

	//Возврат вектора
	return x_next;
}

//Обновление популяции через GBO
void gbo(double** population, int best_ind, int worst_ind, int cur_iter, int m, int n, double pr, double th) {

	//выделяем память для новой популяции
	double** new_population = (double**)malloc(n * sizeof(double*));
	for (int i = 0; i < n; i++) {
		new_population[i] = (double*)malloc(22 * sizeof(double));
	}

	double epsilon = 0.1; //Инициализация epsilon

	//обработка популяции
	for (int cur_vec = 0; cur_vec < n; cur_vec++) {
		double alpha;										   //инициализация alpha
		double* x1 = (double*)malloc(22 * sizeof(double));	   //вектор x1	
		double* x2 = (double*)malloc(22 * sizeof(double));	   //вектор x2
		double* xr1 = (double*)malloc(22 * sizeof(double));	   //1 рандомный вектор 	
		double* xr2 = (double*)malloc(22 * sizeof(double));    //2 рандомный вектор 	
		double* xr3 = (double*)malloc(22 * sizeof(double));    //3 рандомный вектор 	
		double* xr4 = (double*)malloc(22 * sizeof(double));    //4 рандомный вектор 	
		double* x_p = (double*)malloc(22 * sizeof(double));    //рандомный вектор p	
		double* x_rand = (double*)malloc(22 * sizeof(double)); //рандомно сгенерированный вектор
		for (int i = 0; i < 22; i++) {
			x_rand[i] = (double)(-th) + 2 * (double)th * rand_num();
		}
		double* x_next = (double*)malloc(22 * sizeof(double));  //новый вектор

		
		int* indexes = gen_indexes(n, cur_vec, best_ind, worst_ind); //Массив индексов рандомных векторов

		//Копирование  рандомных векторов популяции
		memcpy(xr1, population[indexes[0]], 22 * sizeof(double)); 
		memcpy(xr2, population[indexes[1]], 22 * sizeof(double));
		memcpy(xr3, population[indexes[2]], 22 * sizeof(double));
		memcpy(xr4, population[indexes[3]], 22 * sizeof(double));
		memcpy(x_p, population[indexes[4]], 22 * sizeof(double));

		//Вычисление нового вектора через GSR
		x_next = gsr(population[cur_vec], x1, x2, population[best_ind], population[worst_ind], xr1, xr2, xr3, xr4, cur_iter, m, n, &alpha);

		if (rand_num() < pr) {
			//Вычисление нового вектора через LEO
			new_population[cur_vec] = leo(x_next, population[best_ind], x1, x2, xr1, xr2, x_p, x_rand, alpha, th);
		}
		else {
			//Если LEO не сработал, то просто копируем значение, полученно после GSR
			for (int i = 0; i < 22; i++) {
				if (x_next[i] > th) {
					x_next[i] = th;
				}
				else if (x_next[i] < -th) {
					x_next[i] = -th;
				}
			}
			memcpy(new_population[cur_vec], x_next, 22 * sizeof(double));
		}

		//Освобождение памяти
		free(x1);
		free(x2);
		free(xr1);
		free(xr2);
		free(xr3);
		free(xr4);
		free(x_p);
		free(x_rand);
		free(x_next);
		free(indexes);
	}

	//Обновление текущей популяции и освобождение выделенной памяти
	for (int cur_vec = 0; cur_vec < n; cur_vec++) {
		free(population[cur_vec]);
		population[cur_vec] = malloc(22 * sizeof(double));
		memcpy(population[cur_vec], new_population[cur_vec], 22 * sizeof(double));
		free(new_population[cur_vec]);
	}
	free(new_population);
}
