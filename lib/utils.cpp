#include "utils.hpp"

double* n_vector_generator(unsigned int dim)
{
    double* ptr = new double [dim];

    for (unsigned int i = 0; i < dim; ++i) {
        ptr[i] = 0.0;
    }

    return ptr;
}

void n_vector_destroyer(double* ptr)
{
    delete [] ptr;
}

double** n_matrix_generator(unsigned int dim_1, unsigned int dim_2)
{
    double** ptr = new double* [dim_1];
    for (unsigned int i = 0; i < dim_1; ++i) {
        ptr[i] = new double[dim_2];
    }

    for (unsigned int i = 0; i < dim_1; ++i) {
        for (unsigned int j = 0; j < dim_2; ++j) {
            ptr[i][j] = 0.0;
        }
    }

    return ptr;
}

void n_matrix_destroyer(double** ptr, unsigned int dim_1)
{
    for (unsigned int i = 0; i < dim_1; ++i) {
        delete [] ptr[i];
    }
    delete ptr;
}
