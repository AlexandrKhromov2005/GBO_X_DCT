#ifndef DCT_H
#define DCT_H

extern "C" double a_coef(int k);
extern "C" void dct_func(unsigned char** block, double** dct_block);
extern "C" void rev_dct_func(unsigned char** block, double** dct_block);

#endif // DCT_H