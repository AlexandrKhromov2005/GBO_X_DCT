#include <math.h>

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