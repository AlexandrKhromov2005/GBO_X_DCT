#include <vector>
#include <opencv2/opencv.hpp>
#include <string>
#include <iostream>
#include "dct.h"  

void format_image(const std::string& path) {
    cv::Mat image = cv::imread(path, cv::IMREAD_GRAYSCALE);
    if (image.empty()) {
        std::cerr << "Could not open or find the image: " << path << std::endl;
        return;
    }

    int new_width = (image.cols / 8) * 8;
    int new_height = (image.rows / 8) * 8;

    cv::Rect roi(0, 0, new_width, new_height);
    cv::Mat cropped_image = image(roi);
    cv::imwrite(path, cropped_image);
}

std::vector<double**> split_to_dct_blocks(const std::string& path) {
    std::vector<double**> dct_blocks;
    cv::Mat image = cv::imread(path, cv::IMREAD_GRAYSCALE);

    if (image.empty()) {
        std::cerr << "Could not open or find the image: " << path << std::endl;
        return dct_blocks;
    }

    int width = image.cols;
    int height = image.rows;

    if (width % 8 != 0 || height % 8 != 0) {
        std::cerr << "Image dimensions must be multiples of 8" << std::endl;
        return dct_blocks;
    }

    for (int i = 0; i < height; i += 8) {
        for (int j = 0; j < width; j += 8) {
            unsigned char** temp_block = new unsigned char* [8];
            for (int k = 0; k < 8; k++) {
                temp_block[k] = new unsigned char[8];
                for (int l = 0; l < 8; l++) {
                    temp_block[k][l] = image.at<unsigned char>(i + k, j + l);
                }
            }

            double** dct_block = new double* [8];
            for (int k = 0; k < 8; k++) {
                dct_block[k] = new double[8];
            }

            dct_func(temp_block, dct_block);
            dct_blocks.push_back(dct_block);

            for (int k = 0; k < 8; k++) {
                delete[] temp_block[k];
            }
            delete[] temp_block;
        }
    }

    return dct_blocks;
}

void save_image_from_dct_blocks(const std::vector<double**>& dct_blocks, const std::string& output_path) {
    if (dct_blocks.empty()) {
        std::cerr << "No DCT blocks provided" << std::endl;
        return;
    }

    int num_blocks_width = std::sqrt(dct_blocks.size());
    int num_blocks_height = num_blocks_width;

    int width = num_blocks_width * 8;
    int height = num_blocks_height * 8;

    cv::Mat image(height, width, CV_8UC1);
    int block_index = 0;

    for (int i = 0; i < height; i += 8) {
        for (int j = 0; j < width; j += 8) {
            unsigned char** temp_block = new unsigned char* [8];
            for (int k = 0; k < 8; k++) {
                temp_block[k] = new unsigned char[8];
            }

            rev_dct_func(temp_block, dct_blocks[block_index++]);

            for (int k = 0; k < 8; k++) {
                for (int l = 0; l < 8; l++) {
                    image.at<unsigned char>(i + k, j + l) = temp_block[k][l];
                }
            }

            for (int k = 0; k < 8; k++) {
                delete[] temp_block[k];
            }
            delete[] temp_block;
        }
    }

    cv::imwrite(output_path, image);
}