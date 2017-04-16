#pragma once

long double* n_vector_generator (unsigned int dim);

void n_vector_destroyer(long double* ptr);

long double** n_matrix_generator(unsigned int dim_1, unsigned int dim_2);

void n_matrix_destroyer(long double** ptr, unsigned int dim_1);