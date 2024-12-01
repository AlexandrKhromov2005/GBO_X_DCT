#include <math.h>
//#include <stdio.h>

#define PI acos(-1.0)
#define A_COEF(k) ((k) == 0 ? sqrt(1.0 / 8.0) : sqrt(2.0 / 8.0))

void dct_func(unsigned char block[8][8], double dct_block[8][8]) {

    char temp_block[8][8];
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
                    sum += (double)temp_block[x][y] * cos(thetaU) * cos(thetaV);
                }
            }
            double cu = A_COEF(u);
            double cv = A_COEF(v);
            dct_block[u][v] = sum * cu * cv;
        }
    }
}

void rev_dct_func(unsigned char block[8][8], double dct_block[8][8]) {
    char temp_block[8][8];
    for (int x = 0; x < 8; ++x) {
        for (int y = 0; y < 8; ++y) {
            double sum = 0.0;
            for (int u = 0; u < 8; ++u) {
                for (int v = 0; v < 8; ++v) {
                    double thetaU = (2 * x + 1) * u * PI / (2 * 8);
                    double thetaV = (2 * y + 1) * v * PI / (2 * 8);
                    double cu = A_COEF(u);
                    double cv = A_COEF(v);
                    sum += cu * cv * dct_block[u][v] * cos(thetaU) * cos(thetaV);
                }
            }
            temp_block[x][y] = (char)(round(sum));
        }
    }

    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            block[i][j] = temp_block[i][j] + 128;
        }
    }
}

/*void printBlock(unsigned char block[8][8]) {
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            printf("%3d ", block[i][j]);  // Выводим каждый элемент с отступом в 3 символа
        }
        printf("\n");  // Переход на новую строку после каждой строки массива
    }
}

void printDCTBlock(double block[8][8]) {
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            printf("%3f ", block[i][j]);  // Выводим каждый элемент с отступом в 3 символа
        }
        printf("\n");  // Переход на новую строку после каждой строки массива
    }
}

int main() {
    unsigned char block[8][8];
    unsigned char new_block[8][8];

    double dct_block[8][8];
    unsigned char cnt = 0;
    for(int i = 0; i < 8; i++){
        for(int j = 0; j < 8; j++){
            block[i][j] = cnt;
            cnt++;
        }
    }

    dct_func(block, dct_block);
    printBlock(block);
    printDCTBlock(dct_block);
    rev_dct_func(new_block, dct_block);
    printBlock(new_block);

}*/
