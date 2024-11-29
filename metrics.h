#ifndef METRICS_H
#define METRICS_H

#include <vector>
#define WM_SIZE 1024

void mse(std::vector<unsigned char> vecimg, std::vector<unsigned char> newvecimg);
void ssim(std::vector<unsigned char> vecimg, std::vector<unsigned char> newvecimg);
void ber(unsigned char cvz[WM_SIZE], unsigned char new_cvz[WM_SIZE]);
void ncc(unsigned char cvz[WM_SIZE], unsigned char new_cvz[WM_SIZE]);

#endif // IMAGE_PROCESSING_H
