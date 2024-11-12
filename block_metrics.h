#ifndef BLOCK_METRICS_H
#define BLOCK_METRICS_H

double block_mse(unsigned char** original_block, unsigned char** changed_block);
double block_psnr(unsigned char** original_block, unsigned char** changed_block);

#endif // BLOCK_METRICS_H