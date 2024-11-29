#ifndef DCT_FUNCTIONS_H
#define DCT_FUNCTIONS_H

#ifdef __cplusplus
extern "C" {
#endif
	double a_coef(int k);
	void dct_func(unsigned char block[8][8], double dct_block[8][8]);
	void rev_dct_func(unsigned char block[8][8], double dct_block[8][8]);

#ifdef __cplusplus
}
#endif

#endif // DCT_FUNCTIONS_H
