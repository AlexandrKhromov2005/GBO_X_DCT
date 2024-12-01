#include <math.h>

#define PI acos(-1)

double a_coef(int k) {
    return (k == 0) ? sqrt(1.0 / 8.0) : sqrt(2.0 / 8.0);
}

// Функция для вычисления коэффициента C(u) для DCT
double dct_coef(int u, int n) {
    if (u == 0) {
        return sqrt(1.0 / n);
    }
    else {
        return sqrt(2.0 / n);
    }
}

void dct_func(unsigned char block[8][8], double dct_block[8][8]) {

    double temp_block[8][8];
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            temp_block[i][j] = block[i][j] - 128;
        }
    }

    for (int u = 0; u < 8; ++u) {
        for (int v = 0; v < 8; ++v) {
            double sum = 0.0;
            for (int x = 0; x < 8; ++x) {
                for (int y = 0; y < 8; ++y) {
                    double thetaU = (2 * x + 1) * u * PI / (2 * 8);
                    double thetaV = (2 * y + 1) * v * PI / (2 * 8);
                    sum += temp_block[x][y] * cos(thetaU) * cos(thetaV);
                }
            }
            double cu = dct_coef(u, 8);
            double cv = dct_coef(v, 8);
            dct_block[u][v] = sum * cu * cv;
        }
    }
}

void rev_dct_func(unsigned char block[8][8], double dct_block[8][8]) {
    double temp_block[8][8];
    for (int x = 0; x < 8; ++x) {
        for (int y = 0; y < 8; ++y) {
            double sum = 0.0;
            for (int u = 0; u < 8; ++u) {
                for (int v = 0; v < 8; ++v) {
                    double thetaU = (2 * x + 1) * u * PI / (2 * 8);
                    double thetaV = (2 * y + 1) * v * PI / (2 * 8);
                    double cu = dct_coef(u, 8);
                    double cv = dct_coef(v, 8);
                    sum += cu * cv * dct_block[u][v] * cos(thetaU) * cos(thetaV);
                }
            }
            temp_block[x][y] = (int)(round(sum));
        }
    }

    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            block[i][j] = (unsigned char)temp_block[i][j] + 128;
        }
    }
}
