#include <math.h>

double block_mse(unsigned char original_block[8][8], unsigned char changed_block[8][8]) {
    double mse = 0.0;
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            mse += ((original_block[i][j] - changed_block[i][j]) * (original_block[i][j] - changed_block[i][j]));
        }
    }

    return mse / 64.0;
}

double block_psnr(unsigned char original_block[8][8], unsigned char** changed_block[8][8]) {
    double mse = block_mse(original_block, changed_block);
    if (mse == 0) {
        mse = 0.00001;
    }

    return 10 * log10((255.0 * 255.0) / mse);
}
