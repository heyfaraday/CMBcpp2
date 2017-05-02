#include <cmath>
#include <limits>

#include "utils.hpp"

long double* n_vector_generator(unsigned int dim)
{
    long double* ptr = new long double [dim]();

//    for (unsigned int i = 0; i < dim; ++i) {
//        ptr[i] = 0.0L;
//    }

    return ptr;
}

int* n_int_vector_generator(unsigned int dim)
{
    int* ptr = new int [dim]();

//    for (unsigned int i = 0; i < dim; ++i) {
//        ptr[i] = 0;
//    }

    return ptr;
}

void n_vector_destroyer(long double* ptr)
{
    delete [] ptr;
}

void n_int_vector_destroyer(int* ptr)
{
    delete [] ptr;
}

long double** n_matrix_generator(unsigned int dim_1, unsigned int dim_2)
{
    long double** ptr = new long double* [dim_1];
    for (unsigned int i = 0; i < dim_1; ++i) {
        ptr[i] = new long double[dim_2]();
    }

//    for (unsigned int i = 0; i < dim_1; ++i) {
//        for (unsigned int j = 0; j < dim_2; ++j) {
//            ptr[i][j] = 0.0L;
//        }
//    }

    return ptr;
}

int** n_int_matrix_generator(unsigned int dim_1, unsigned int dim_2)
{
    int** ptr = new int* [dim_1];
    for (unsigned int i = 0; i < dim_1; ++i) {
        ptr[i] = new int[dim_2]();
    }

//    for (unsigned int i = 0; i < dim_1; ++i) {
//        for (unsigned int j = 0; j < dim_2; ++j) {
//            ptr[i][j] = 0.0L;
//        }
//    }

    return ptr;
}

void n_matrix_destroyer(long double** ptr, unsigned int dim_1)
{
    for (unsigned int i = 0; i < dim_1; ++i) {
        delete [] ptr[i];
    }
    delete ptr;
}

void n_int_matrix_destroyer(int** ptr, unsigned int dim_1)
{
    for (unsigned int i = 0; i < dim_1; ++i) {
        delete [] ptr[i];
    }
    delete ptr;
}

long double det(long double q1, long double q2, long double u1, long double u2) {
    return q1 * u2 - q2 * u1;
}

long double mean_t_solver(long double q1, long double q2, long double u1, long double u2, long double x) {
    return u2 * x * x * x + (u1 + 2.0L * q2) * x * x + (2.0L * q1 - u2) * x - u1;
}

void d_solver(long double& root1, long double& root2,
              long double q1, long double q2, long double u1, long double u2) {

    long double a = 3.0L * u2;
    long double b = 2.0L * (u1 + 2.0L * q2);
    long double c = (2.0L * q1 - u2);

    long double determinant = b * b - 4.0L * a * c;

    if (determinant > 0.0L) {
        root1 = (-b + sqrtl(determinant)) / (2.0L * a);
        root2 = (-b - sqrtl(determinant)) / (2.0L * a);
    } else {
        root1 = 0;
        root2 = 0;
    }
}

bool is_equal(long double x, long double y) {
    return std::fabsl(x - y) < std::numeric_limits<long double>::epsilon();
}
