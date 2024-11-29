#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#define POPSIZE 100 

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
void gen_indexes(int indexes[5], int n, int cur_ind, int best_ind, int worst_ind) {
	int cnt = 0;
	while (cnt < 5) {
		int temp = rand() % n;
		if ((indexes[0] != temp) && (indexes[1] != temp) && (indexes[2] != temp) && (indexes[3] != temp) && (indexes[4] != temp) && (cur_ind != temp)) {
			indexes[cnt] = temp;
			cnt++;
		}
	}
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


	//Вычисление delta
	double delta[22];

	//Вычисление step
	double step [22];

	//Вычиследние delta X
	double x_delta [22];

	//Вычисление Z
	double z [22];

	//Вычисление P
	double p [22];

	//Вычисление Q
	double q [22];

	//Вычисление вектора следующей популяции
	for (int i = 0; i < 22; i++) {
		delta[i] = 2.0 * rand_num() * fabs(((xr1[i] + xr2[i] + xr3[i] + xr4[i]) / 4.0) - x[i]);

		step[i] = ((x_best[i] - xr1[i]) + delta[i]) / 2.0;

		x_delta[i] = ((double)(rand() % n) + 1) * fabs(step[i]);

		z[i] = x[i] - randn() * ((2 * x_delta[i] * x[i]) / (x_worst[i] - x_best[i] + epsilon));

		p[i] = rand_num() * (((z[i] + x[i]) / 2.0) + rand() * x_delta[i]);

		q[i] = rand() * (((z[i] + x[i]) / 2.0) - rand_num() * x_delta[i]);

		x1[i] = x[i] - randn() * rho1 * ((2 * x_delta[i] * x[i]) / (p[i] - q[i] + epsilon)) + rand_num() * rho2 * (x_best[i] - x[i]);

		x2[i] = x_best[i] - randn() * rho1 * ((2 * x_delta[i] * x[i]) / (p[i] - q[i] + epsilon)) + rand_num() * rho2 * (xr1[i] - xr2[i]);

		x3[i] = x[i] - rho1 * (x2[i] - x1[i]);

		x_next[i] = ra * (rb * x1[i] + (1 - rb) * x[i]) + (1 - ra) * x3[i];
	}
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
		if (x_next[i] > th) {
			x_next[i] = th;
		}
		else if (x_next[i] < -th) {
			x_next[i] = -th;
		}
	}
}

//Обновление популяции через GBO
void gbo(double population[][22], int best_ind, int worst_ind, int cur_iter, int m, int n, double pr, double th) {

	//выделяем память для новой популяции
	double new_population [POPSIZE][22];

	double epsilon = 0.1; //Инициализация epsilon

	//обработка популяции
	for (int cur_vec = 0; cur_vec < n; cur_vec++) {
		double alpha;										   //инициализация alpha
		double x1 [22];	   //вектор x1	
		double x2 [22];	   //вектор x2
		double xr1 [22];	   //1 рандомный вектор 	
		double xr2 [22]; //2 рандомный вектор 	
		double xr3 [22];    //3 рандомный вектор 	
		double xr4 [22];    //4 рандомный вектор 	
		double x_p [22];    //рандомный вектор p	
		double x_rand [22]; //рандомно сгенерированный вектор
		for (int i = 0; i < 22; i++) {
			x_rand[i] = (double)(-th) + 2 * (double)th * rand_num();
		}
		double x_next [22];  //новый вектор

		int indexes[5];
		gen_indexes(indexes, n, cur_vec, best_ind, worst_ind); //Массив индексов рандомных векторов

		//Копирование  рандомных векторов популяции
		memcpy(xr1, population[indexes[0]], 22); 
		memcpy(xr2, population[indexes[1]], 22);
		memcpy(xr3, population[indexes[2]], 22);
		memcpy(xr4, population[indexes[3]], 22);
		memcpy(x_p, population[indexes[4]], 22);

		//Вычисление нового вектора через GSR
		gsr(x_next, population[cur_vec], x1, x2, population[best_ind], population[worst_ind], xr1, xr2, xr3, xr4, cur_iter, m, n, &alpha);

		if (rand_num() < pr) {
			//Вычисление нового вектора через LEO
			leo(new_population[cur_vec], x_next, population[best_ind], x1, x2, xr1, xr2, x_p, x_rand, alpha, th);
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
			memcpy(new_population[cur_vec], x_next, 22);
		}
	}

	//Обновление текущей популяции и освобождение выделенной памяти
	for (int cur_vec = 0; cur_vec < n; cur_vec++) {
		memset(population[cur_vec], 0, 22);
		memcpy(population[cur_vec], new_population[cur_vec], 22);
		memset(new_population[cur_vec], 0, 22);
	}
}
