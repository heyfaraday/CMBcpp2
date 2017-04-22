#include <iostream>
#include <fstream>
#include <chealpix.h>
#include <cmath>

#include <fft.hpp>
#include <utils.hpp>
#include <parameters.hpp>
#include "constants.hpp"
#include <io.hpp>
#include "pml.hpp"
#include "aml.hpp"
#include "functionals.hpp"

int main() {

    typedef std::numeric_limits<long double> dbl;
    std::cout.precision(dbl::max_digits10);

    long double** x_1_ml = n_matrix_generator(nmod, nmod);
    long double** x_2_ml = n_matrix_generator(nmod, nmod);
    long double** pml = n_matrix_generator(nmod, nmod);

    long double sum = 0.0L;
    long double norm = 0.0L;
    unsigned int m = 20;
    unsigned int l = 50;
    long double theta;
//    for (unsigned int j = 1; j < npix / 2; ++j) {
//        theta = 2.0L * j * PI / long_npix;
//        x_1_ml_gen(j, x_1_ml);
//        x_2_ml_gen(j, x_2_ml);
//        sum += (x_1_ml[m][l] + x_2_ml[m][l]) * (x_1_ml[m+1][l+1] + x_2_ml[m+1][l+1]) * sinl(theta);
//        norm += sinl(theta);
//    }


    x_1_ml_gen(30, x_2_ml);
    pml_gen(PI/3.0L, pml);

    std::cout << x_2_ml[4][5] << std::endl;
    std::cout << pml[4][5] << std::endl;
    std::cout << pml[4][5-1] << std::endl;

    n_matrix_destroyer(x_1_ml, nmod);
    n_matrix_destroyer(x_2_ml, nmod);
    n_matrix_destroyer(pml, nmod);

    return 0;
}