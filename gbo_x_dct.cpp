#include <iostream>
#include "stdio.h"
#include "image_processing.h"
#include "population.h"
#include "gbo.h"
#include "metrics.h"
#include "dct_functions.h"
#include <time.h>
#include "block_metrics.h"
#include <random>
#include <vector>
#define POPSIZE 100
#define WM_SIZE 1024

int main() {
	srand(time(NULL));

	int mode;
	std::cout << "Enter mode: \n0)Embedding WM \n1)Extraction WM \n2)Metrics\n";
	std::cin >> mode;

	//Встраивание ЦВЗ
	if (mode == 0) {
		//Импорт изображения
		std::string image_path;
		std::cout << "Enter the path of the image:\n";
		std::cin >> image_path;
		std::vector<std::vector<std::vector<double>>> dct_blocks = split_to_dct_blocks(image_path);
		std::vector<std::vector<std::vector<double>>> new_dct_blocks;

		std::string new_image_path;
		std::string	wm_path;
		int n, m, th;
		double pr;

		//Импорт ЦВЗ
		std::cout << "Enter the path of the watermark:\n";
		std::cin >> wm_path;
		unsigned char wm[WM_SIZE];
		get_wm_matrix(wm_path, wm);

		std::cout << "Enter the path to save the image:\n";
		std::cin >> new_image_path;

		//Импорт  N, M, Th, pr
		std::cout << "Enter size of population:\n";
		std::cin >> n;
		std::cout << "Enter number of iterations:\n";
		std::cin >> m;
		std::cout << "Enter limit of values:\n";
		std::cin >> th;
		std::cout << "Enter probability:\n";
		std::cin >> pr;

		//Итерация по каждому ДКП блоку
		for (int cur_block = 0; cur_block < dct_blocks.size(); cur_block++) {
		
			//Создание инициализирующей популяции
			double population[POPSIZE][22];
			create_population(th, n, population);

			double block[8][8];
			double block_new[8][8];
			from_vec_to_list(block, dct_blocks[cur_block]);

			//Оптимизауия векторов с использованием GBO	
			gbo(population, m, n, pr, th, block, wm[cur_block % WM_SIZE]);

			//Встраивание лучшего вектора
			int best_ind = -1;
			best_ind = find_x_best(population, block, n, wm[cur_block % WM_SIZE]);
			apply_x(block, population[best_ind], block_new);

			//Записываем его в вектор нового изображения
			save_new_block(new_dct_blocks, block_new);
		
			//Вывод проделанной работы
			std::cout << "Completed " << ((double)(cur_block + 1 ) / dct_blocks.size()) * 100 << "%\n";
		}

		//Сохрание нового изображения
		save_image_from_dct_blocks(new_dct_blocks, new_image_path);
		std::cout << "Image saved.\n";
	}

	else if (mode == 1) {

		//Импорт изображения
		std::string image_path;
		std::cout << "Enter the path of the image:\n";
		std::cin >> image_path;
		std::vector<std::vector<std::vector<double>>> dct_blocks = split_to_dct_blocks(image_path);

		//Путь для сохранения извлеченного ЦВЗ
		std::string wm_path;
		std::cout << "Enter the path to save the watermark:\n";
		std::cin >> wm_path;
		std::vector<unsigned char> wm_bits;

		//Итерация по ДКП блокам изображения и извлечение битов ЦВЗ
		for (int cur_block = 0; cur_block < dct_blocks.size(); cur_block++) {

			double block[8][8];
			from_vec_to_list(block, dct_blocks[cur_block]);

			double s1 = get_s1_sum(block);
			double s0 = get_s0_sum(block);

			if (s1 < s0) {
				wm_bits.push_back(0);
			}
			else{ wm_bits.push_back(1); }
		}

		//Процесс голосования и формирование извлеченного ЦВЗ
		unsigned char wm[WM_SIZE];
		for(int i = 0; i < WM_SIZE; i++){
			unsigned char zero = 0;
			unsigned char one = 0;

			for (int j = 0; j < wm_bits.size(); j++) {
				if (j % WM_SIZE == i) {
					if (wm_bits[j] == 0) {
						zero++;
					}
					else { one++; }
				}
			}

			if (zero > one) {
				wm[i] = 0;
			}
			else { wm[i] = 1; }
		}

		//Вывод ЦВЗ
		std::cout << "WATERMARK\n";
		for (int i = 0; i < sqrt(WM_SIZE); i++) {
			for (int j = 0; j < sqrt(WM_SIZE); j++) {
				std::cout << (int)wm[i * (int)sqrt(WM_SIZE) + j] << " ";
			}
			std::cout << "\n";
		}

		//Сохранение извлеченного ЦВЗ
		save_wm_matrix(wm_path, wm);
		std::cout << "Watermark saved.\n";
	}

	else if (mode == 2) {

		//Импорт исходного изображения
		std::string image_path, new_image_path, wm_path, new_wm_path;
		std::cout << "Enter the path of the image:\n";
		std::cin >> image_path;
		std::vector<unsigned char> img_pixels = img_to_vec(image_path);

		//Импорт изображения, в которое встроен ЦВЗ
		std::cout << "Enter the path of the new image:\n";
		std::cin >> new_image_path;
		std::vector<unsigned char> new_img_pixels = img_to_vec(new_image_path);

		//Импорт исходного ЦВЗ
		std::cout << "Enter the path of the watermark:\n";
		std::cin >> wm_path;
		unsigned char wm[WM_SIZE];
		get_wm_matrix(wm_path, wm);

		//Импорт извлеченного ЦВЗ
		std::cout << "Enter the path of the extracted watermark:\n";
		std::cin >> new_wm_path;
		unsigned char new_wm[WM_SIZE];
		get_wm_matrix(new_wm_path, new_wm);

		//Вычисление метрик
		mse(img_pixels, new_img_pixels);
		ssim(img_pixels, new_img_pixels);
		ber(wm, new_wm);
		ncc(wm, new_wm);
	}
	std::cout << "Press any key to close this window...";
	std::getchar();
}
