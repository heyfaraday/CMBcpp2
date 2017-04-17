#pragma once

long double* n_vector_generator (unsigned int dim);

void n_vector_destroyer(long double* ptr);

long double** n_matrix_generator(unsigned int dim_1, unsigned int dim_2);

void n_matrix_destroyer(long double** ptr, unsigned int dim_1);

long double det(long double q1, long double q2, long double u1, long double u2);

long double mean_t_solver(long double q1, long double q2, long double u1, long double u2, long double x);

void d_solver(long double* root1, long double* root2,
              long double q1, long double q2, long double u1, long double u2);
