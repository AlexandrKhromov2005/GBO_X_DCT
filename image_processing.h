#ifndef IMAGE_PROCESSING_H
#define IMAGE_PROCESSING_H

#include <string>
#include <vector>
#define WM_SIZE 1024

void format_image(const std::string& path);
std::vector<std::vector<std::vector<double>>> split_to_dct_blocks(const std::string& path);
void save_image_from_dct_blocks(const std::vector<std::vector<std::vector<double>>>& dct_blocks, const std::string& output_path);
void save_new_block(std::vector<std::vector<std::vector<double>>>& new_dct_blocks, double block[8][8]);
void from_vec_to_list(double block[8][8], std::vector<std::vector<double>>& vec_block);
void get_wm_matrix(std::string& path, unsigned char wm_matrix[WM_SIZE]);
void save_wm_matrix(const std::string& path, const unsigned char wm_matrix[WM_SIZE]);
std::vector<unsigned char> img_to_vec(const std::string& path);

#endif // IMAGE_PROCESSING_H
