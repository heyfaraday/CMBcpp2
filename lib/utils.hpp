#pragma once

long double* n_vector_generator (unsigned int dim);

int* n_int_vector_generator(unsigned int dim);

void n_vector_destroyer(long double* ptr);

void n_int_vector_destroyer(int* ptr);

long double** n_matrix_generator(unsigned int dim_1, unsigned int dim_2);

int** n_int_matrix_generator(unsigned int dim_1, unsigned int dim_2);

void n_matrix_destroyer(long double** ptr, unsigned int dim_1);

void n_int_matrix_destroyer(int** ptr, unsigned int dim_1);

long double det(long double q1, long double q2, long double u1, long double u2);

long double mean_t_solver(long double q1, long double q2, long double u1, long double u2, long double x);

void d_solver(long double& root1, long double& root2,
              long double q1, long double q2, long double u1, long double u2);

bool is_equal(long double x, long double y);
