#include <math.h>

#define PI acos(-1)

double a_coef(int k) {
    return (k == 0) ? sqrt(1.0 / 8.0) : sqrt(2.0 / 8.0);
}

void dct_func(unsigned char** block, double** dct_block) {
    for (int u = 0; u < 8; u++) {
        for (int v = 0; v < 8; v++) {
            double dct_coef = 0.0;

            for (int x = 0; x < 8; x++) {
                for (int y = 0; y < 8; y++) {
                    dct_coef += block[x][y] * cos(((2 * x + 1) * u * PI) / 16) * cos(((2 * y + 1) * v * PI) / 16);
                }
            }

            dct_block[u][v] = dct_coef * (a_coef(v) * a_coef(u));
        }
    }
}

void rev_dct_func(unsigned char** block, double** dct_block) {
    for (int x = 0; x < 8; x++) {
        for (int y = 0; y < 8; y++) {
            double sum = 0.0;

            for (int u = 0; u < 8; u++) {
                for (int v = 0; v < 8; v++) {
                    sum += a_coef(u) * a_coef(v) * dct_block[u][v] * cos(((2 * x + 1) * u * PI) / 16) * cos(((2 * y + 1) * v * PI) / 16);
                }
            }

            block[x][y] = (unsigned char)(sum + 0.5);
        }
    }
}

