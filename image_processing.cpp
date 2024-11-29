#include <vector>
#include <array>
#include <opencv2/opencv.hpp>
#include <string>
#include <iostream>
#include "dct_functions.h"  

#define WM_SIZE 1024

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

std::vector<unsigned char> img_to_vec(const std::string& path) {
    std::vector<unsigned char> pixels;
    cv::Mat image = cv::imread(path, cv::IMREAD_GRAYSCALE);

    if (image.empty()) {
        std::cerr << "Could not open or find the image: " << path << std::endl;
        return pixels;
    }

    int width = image.cols;
    int height = image.rows;

    if (width % 8 != 0 || height % 8 != 0) {
        std::cerr << "Image dimensions must be multiples of 8" << std::endl;
        return pixels;
    }

    for (int i = 0; i < height; i += 8) {
        for (int j = 0; j < width; j += 8) {
            pixels.push_back(image.at<unsigned char>(j, j));         
        }
    }

    return pixels;
}

void get_wm_matrix(std::string& path, unsigned char wm_matrix[WM_SIZE]) {
    cv::Mat wm = cv::imread(path, cv::IMREAD_GRAYSCALE);
    for (int i = 0; i < sqrt(WM_SIZE); i++) {
        for (int j = 0; j < sqrt(WM_SIZE); j++) {
            wm_matrix[i * (int)sqrt(WM_SIZE) + j] = wm.at<unsigned char>(i, j) == 255 ? 1 : 0;
        }
    }
}

void save_wm_matrix(const std::string& path, const unsigned char wm_matrix[WM_SIZE]) {
    cv::Mat wm(sqrt(WM_SIZE), sqrt(WM_SIZE), CV_8UC1, cv::Scalar(0));

    for (int i = 0; i < sqrt(WM_SIZE); i++) {
        for (int j = 0; j < sqrt(WM_SIZE); j++) {
            wm.at<unsigned char>(i, j) = wm_matrix[i * (int)sqrt(WM_SIZE) + j] == 1 ? 255 : 0;
        }
    }

    cv::imwrite(path, wm);
}

std::vector<std::vector<std::vector<double>>> split_to_dct_blocks(const std::string& path) {
    std::vector<std::vector<std::vector<double>>> dct_blocks;
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
            unsigned char temp_block[8][8];
            for (int k = 0; k < 8; k++) {
                for (int l = 0; l < 8; l++) {
                    temp_block[k][l] = image.at<unsigned char>(i + k, j + l);
                }
            }

            double dct_block[8][8] = {};  // Инициализация нулями
            dct_func(temp_block, dct_block);
            std::vector<std::vector<double>> dct_block_array(8, std::vector<double>(8, 0.0));  // Инициализация вектора
            for (int k = 0; k < 8; k++) {
                for (int l = 0; l < 8; l++) {
                    dct_block_array[k][l] = dct_block[k][l];
                }
            }
            dct_blocks.push_back(dct_block_array);
        }
    }

    return dct_blocks;
}

void save_image_from_dct_blocks(const std::vector<std::vector<std::vector<double>>>& dct_blocks, const std::string& output_path) {
    if (dct_blocks.empty()) {
        std::cerr << "No DCT blocks provided" << std::endl;
        return;
    }

    int num_blocks_width = std::sqrt(dct_blocks.size());
    int num_blocks_height = num_blocks_width;

    int width = num_blocks_width * 8;
    int height = num_blocks_height * 8;

    cv::Mat image(height, width, CV_8UC1, cv::Scalar(0));
    int block_index = 0;

    for (int i = 0; i < height; i += 8) {
        for (int j = 0; j < width; j += 8) {
            // Преобразование std::array<std::array<double, 8>, 8> в double[8][8]
            double dct_block[8][8];
            for (int k = 0; k < 8; k++) {
                for (int l = 0; l < 8; l++) {
                    dct_block[k][l] = dct_blocks[block_index][k][l];
                }
            }

            unsigned char temp_block[8][8] = { 0 };
            rev_dct_func(temp_block, dct_block);

            for (int k = 0; k < 8; k++) {
                for (int l = 0; l < 8; l++) {
                    image.at<unsigned char>(i + k, j + l) = temp_block[k][l];
                }
            }

            block_index++;
        }
    }

    cv::imwrite(output_path, image);
}

void from_vec_to_list(double block[8][8], std::vector<std::vector<double>>& vec_block) {
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            block[i][j] = vec_block[i][j];
        }
    }
}

void save_new_block(std::vector<std::vector<std::vector<double>>>& new_dct_blocks, double block[8][8]) {
    std::vector<std::vector<double>> vec_block;
    for (int i = 0; i < 8; ++i) {
        std::vector<double> row;
        for (int j = 0; j < 8; ++j) {
            row.push_back(block[i][j]);
        }
        vec_block.push_back(row);
    }

    new_dct_blocks.push_back(vec_block);
}
