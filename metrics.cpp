#include <iostream>
#include <vector>
#define WM_SIZE 1024

void mse(std::vector<unsigned char> vecimg, std::vector<unsigned char> newvecimg) {

    float pr = 0;

    for (int i = 0; i < size(vecimg); i++) {
        pr += (vecimg[i] - newvecimg[i]) * (vecimg[i] - newvecimg[i]);
    }

    float mse = pr / size(vecimg);

    float rmse = sqrt(mse);

    double psnr = 10 * log10((255 * 255) / mse);

    std::cout << "MSE = " << mse << "\n";

    std::cout << "RMSE = " << rmse << "\n";

    std::cout << "PSNR = " << psnr << "\n";
}

void ssim(std::vector<unsigned char> vecimg, std::vector<unsigned char> newvecimg) {

    float k1 = 0.01;

    float k2 = 0.03;

    float u_p = 0.0;

    float u_s = 0.0;

    float s1 = 0.0;

    float s2 = 0.0;

    float cov = 0.0;

    for (int i = 0; i < size(vecimg); i++) {

        u_p += float(vecimg[i]);

        u_s += float(newvecimg[i]);
    }

    u_p /= float(size(vecimg));

    u_s /= float(size(newvecimg));

    for (int i = 0; i < size(vecimg); i++) {
        s1 += (float(vecimg[i] - u_p)) * (float(vecimg[i] - u_p));

        s2 += (float(newvecimg[i] - u_s)) * (float(newvecimg[i] - u_s));
    }

    s1 /= float(size(vecimg));
    s2 /= float(size(newvecimg));

    for (int i = 0; i < size(vecimg); i++) {
                cov += (float((vecimg[i] - u_p)) * (float(newvecimg[i] - u_s))) / (float(size(vecimg)));
    }
    float ssim = (((2 * u_p * u_s) + k1) * ((2 * cov) + k2)) / (((u_p * u_p) + (u_s * u_s) + k1) * (s1 + s2 + k2));

    std::cout << "SSIM = " << ssim << "\n";
}

void ber(unsigned char cvz[WM_SIZE], unsigned char new_cvz[WM_SIZE]) {
    int err = 0;
    for (int i = 0; i < WM_SIZE; i++) {
        if (new_cvz[i] != cvz[i]) {
            err++;
        }
    }

    double ber = (double)err / WM_SIZE;

    std::cout << "BER = " << ber << "\n";
}

void ncc(unsigned char cvz[WM_SIZE], unsigned char new_cvz[WM_SIZE]) {

    float ncc = 0;
    float chisl = 0;
    float zn1 = 0;
    float zn2 = 0;

    for (int i = 0; i < WM_SIZE; i++) {
        chisl += cvz[i] * new_cvz[i];
        zn1 += cvz[i] * cvz[i];
        zn2 += new_cvz[i] * new_cvz[i];
    }
    ncc = chisl / (sqrt(zn1) * sqrt(zn2));

    std::cout << "NCC = " << ncc << "\n";
}
