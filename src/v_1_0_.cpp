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


    for (unsigned int global_monter_carlo = 0; global_monter_carlo < nmonte; ++global_monter_carlo) {
        std::cout << "Global start: " << global_monter_carlo << std::endl;



        std::cout << "Global end: " << global_monter_carlo << std::endl;
    }

    long double** cos_ml = n_matrix_generator(nback, nback);
    i_aml("../data/alm_Q_norm_2048.dat", cos_ml);
    n_matrix_destroyer(cos_ml, nback);
    std::cout << cos_ml[5][10];


    return 0;
}
