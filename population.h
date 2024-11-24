#ifndef POPULATION_H
#define POPULATION_H

extern "C" double rand_double();
extern "C" double sign(double n);
extern "C" double** create_population(int th, int d);
extern "C" double get_s0_sum(double** block);
extern "C" double get_s1_sum(double** block);
extern "C" double of(double** original_block, double** changed_block, char b);
extern "C" int find_x_best(double** population, double** original_dct_block, int popul_size, const char b);
extern "C" int find_x_worst(double** population, double** original_dct_block, int popul_size, const char b);
extern "C" double** apply_x(double** block, double* x);

#endif // POPULATION_H
