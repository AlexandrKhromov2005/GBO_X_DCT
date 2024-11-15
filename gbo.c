#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"



#define PI acos(-1)

double rand_gauss_calc(double mu, double sigma) {
    double u1, u2, z;
    do {
        u1 = (double)rand() / RAND_MAX;
        u2 = (double)rand() / RAND_MAX;
    } while (u1 <= 1e-7);  // Исключаем случай, когда u1 очень маленькое

    z = sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
    return mu + sigma * z;
}

void leo(double* x_next, double* x_best, double* x_rand1, double* x_rand2, double* x_rand_p, double* x_rand_gen, double* x1, double* x2, double rho1) {
    double* y = (double*)malloc(22 * sizeof(double));
    double* x_k = (double*)malloc(22 * sizeof(double));

    double rand_num = (double)rand() / RAND_MAX;
    if (rand_num < 0.5) {
        memcpy(y, x_next, 22 * sizeof(double));
    }
    else {
        memcpy(y, x_best, 22 * sizeof(double));
    }

    double f1 = -1.0 + 2.0 * ((double)rand() / RAND_MAX);
    double f2 = -1.0 + 2.0 * ((double)rand() / RAND_MAX);

    double mu1 = (double)rand() / RAND_MAX;
    double u1 = 0;
    double u2 = 0;
    double u3 = 0;

    rand_num = (double)rand() / RAND_MAX;
    if (mu1 < 0.5) {
        rand_num = (double)rand() / RAND_MAX;
        u1 = 2 * rand_num;

        rand_num = (double)rand() / RAND_MAX;
        u2 = rand_num;

        rand_num = (double)rand() / RAND_MAX;
        u3 = rand_num;
    }
    else {
        u1 = 1;
        u2 = 1;
        u3 = 1;
    }

    double mu2 = (double)rand() / RAND_MAX;
    if (mu2 < 0.5) {
        memcpy(x_k, x_rand_gen, 22 * sizeof(double));
    }
    else {
        memcpy(x_k, x_rand_p, 22 * sizeof(double));
    }

    for (int i = 0; i < 22; i++) {
        x_next[i] = y[i] + f1 * (u1 * x_best[i] - u2 * x_k[i]) + f2 * rho1 * (u3 * (x2[i] - x1[i]) + u2 * (x_rand1[i] - x_rand2[i])) / 2;
    }

    free(y);
    free(x_k);
}

void gbo(double* x_next, double* x, double* x_best, double* x_worst, double* x_rand1, double* x_rand2, double* x_rand_p, double* x_rand_gen, unsigned char cur_iter, unsigned char iter_number, double pr) {
   
    double beta_min = 0.2; 
    double beta_max = 1.2;
    double* x1 = (double*)malloc(22 * sizeof(double));
    double* x2 = (double*)malloc(22 * sizeof(double));
    double* x3 = (double*)malloc(22 * sizeof(double));
    double* x_k = (double*)malloc(22 * sizeof(double));
    double x_delta = 0;
    double v_delta = 0;
    double p = 0;
    double q = 0;

    double beta = beta_min + (beta_max - beta_min) * pow((1 - pow((double)cur_iter / iter_number, 3)), 2);
    double alpha = fabs(beta * sin((3 * PI) / 2 + sin((3 * PI * beta) / 2)));

    double rand_num = (double)rand() / RAND_MAX;
    double rho1 = 2 * rand_num * alpha - alpha;

    rand_num = (double)rand() / RAND_MAX;
    double rho2 = 2 * rand_num * alpha - alpha;

    double randn = rand_gauss_calc(0.0, 1.0);
    double epsilon = (double)rand() / (RAND_MAX / 0.1);

    double r_a = (double)rand() / RAND_MAX;
    double r_b = (double)rand() / RAND_MAX;

    //повторяем действие для всех элементов вектора
    for (int i = 0; i < 22; i++) {
        x_delta = fabs(x_best[i] - x_rand1[i]);

        v_delta = fabs(x_best[i] - x_worst[i]);
        p = x_delta + v_delta;
        q = x_delta - v_delta;

        /*p = x_delta - x_best[i];
        q = x_delta - x_worst[i];*/
        rand_num = (double)rand() / RAND_MAX;
        x1[i] = x[i] + randn * rho1 * ((2 * x_delta * x[i]) / (p - q + epsilon) + rand_num * rho2 * (x_best[i] - x[i]));

        rand_num = (double)rand() / RAND_MAX;
        x2[i] = x_best[i] + randn * rho1 * ((2 * x_delta * x[i]) / (p - q + epsilon) + rand_num * rho2 * (x_rand1[i] - x_rand2[i]));
        
        x3[i] = x[i] - rho1 * (x2[i] - x1[i]);

        x_next[i] = r_a * (r_b * x1[i] + (1 - r_b) * x2[i]) + (1 - r_a) * x3[i];
    }

    rand_num = (double)rand() / RAND_MAX;
    if (rand_num < pr) {
        leo(x_next, x_best, x_rand1, x_rand2, x_rand_p, x_rand_gen, x1, x2, rho1);
    }

    free(x1);
    free(x2);
    free(x3);
    free(x_k);
}

void x_next(double** population, int cur_vec, int best_ind, int worst_ind, unsigned char cur_iter, unsigned char iter_number, double pr, int th, int n) {
    double* gen_vec = (double*)malloc(22 * sizeof(double));

    for (int j = 0; j < 22; j++) {
        gen_vec[j] = ((double)rand() / RAND_MAX) * (2 * th) - th;
    }

    // Генерация случайного числа в диапазоне [0, n]
    int rand1_ind = rand() % (n);
    int rand2_ind = rand() % (n);
    int rand_p_ind = rand() % (n);

    gbo(population[cur_vec], population[cur_vec], population[best_ind], population[worst_ind], population[rand1_ind],
        population[rand2_ind], population[rand_p_ind], gen_vec, cur_iter, iter_number, pr);

    

    free(gen_vec);
}
