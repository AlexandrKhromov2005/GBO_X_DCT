#ifndef BLOCK_METRICS_H
#define BLOCK_METRICS_H

#ifdef __cplusplus
extern "C" {
#endif
	double block_mse(unsigned char original_block[8][8], unsigned char changed_block[8][8]);
	double block_psnr(unsigned char original_block[8][8], unsigned char changed_block[8][8]);
#ifdef __cplusplus
}
#endif

#endif // BLOCK_METRICS_H
