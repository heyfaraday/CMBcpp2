#include "utils.hpp"

long double* n_vector_generator(unsigned int dim)
{
    long double* ptr = new long double [dim];

    for (unsigned int i = 0; i < dim; ++i) {
        ptr[i] = 0.0;
    }

    return ptr;
}

void n_vector_destroyer(long double* ptr)
{
    delete [] ptr;
}

long double** n_matrix_generator(unsigned int dim_1, unsigned int dim_2)
{
    long double** ptr = new long double* [dim_1];
    for (unsigned int i = 0; i < dim_1; ++i) {
        ptr[i] = new long double[dim_2];
    }

    for (unsigned int i = 0; i < dim_1; ++i) {
        for (unsigned int j = 0; j < dim_2; ++j) {
            ptr[i][j] = 0.0;
        }
    }

    return ptr;
}

void n_matrix_destroyer(long double** ptr, unsigned int dim_1)
{
    for (unsigned int i = 0; i < dim_1; ++i) {
        delete [] ptr[i];
    }
    delete ptr;
}
