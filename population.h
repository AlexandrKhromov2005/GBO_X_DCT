#ifndef POPULATION_H
#define POPULATION_H

#ifdef __cplusplus
extern "C" {
#endif
    typedef struct {
        int best;
        int worst;
        double f_values[22];
    }xind;

double rand_double();
double sign(double n);
void create_population(int th, int d, double population[][22]);
double get_s0_sum(double block[8][8]);
double get_s1_sum(double block[8][8]);
double of(double original_dct_block[8][8], double changed_dct_block[8][8], double vector[22], const char b);
int find_x_best(double population[][22], double original_dct_block[8][8], int popul_size, const char b);
int find_x_worst(double population[][22], double original_dct_block[8][8], int popul_size, const char b);
void apply_x(double block[8][8], double x[22], double new_block[8][8]);
xind find_x_bw(double population[][22], double original_dct_block[8][8], int popul_size, const char b);
#ifdef __cplusplus
}
#endif

#endif // POPULATION_H
