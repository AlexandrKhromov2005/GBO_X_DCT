#ifndef IMAGE_PROCESSING_H
#define IMAGE_PROCESSING_H

#include <string>
#include <vector>

void format_image(const std::string& path);
std::vector<double**> split_to_dct_blocks(const std::string& path);
void save_image_from_dct_blocks(const std::vector<double**>& dct_blocks, const std::string& output_path);

#endif // IMAGE_PROCESSING_H