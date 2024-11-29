#ifndef POPULATION_H
#define POPULATION_H

extern "C" double rand_double();
extern "C" double sign(double n);
extern "C" void create_population(int th, int d, double population[][22]);
extern "C" double get_s0_sum(double block[8][8]);
extern "C" double get_s1_sum(double block[8][8]);
extern "C" double of(double original_dct_block[8][8], double changed_dct_block[8][8], double vector[22], const char b);
extern "C" int find_x_best(double population[][22], double original_dct_block[8][8], int popul_size, const char b);
extern "C" int find_x_worst(double population[][22], double original_dct_block[8][8], int popul_size, const char b);
extern "C" void apply_x(double block[8][8], double x[22], double new_block[8][8]);

#endif // POPULATION_H
