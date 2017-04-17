#include <iostream>
#include <fstream>
#include <chealpix.h>
#include "math.h"

#include <fft.hpp>
#include <utils.hpp>
#include <parameters.hpp>
#include "constants.hpp"
#include <io.hpp>
#include "pml.hpp"
#include "aml.hpp"
#include "functionals.hpp"

int main(int argc, char const *argv[]) {

    typedef std::numeric_limits<long double> dbl;
    std::cout.precision(dbl::max_digits10);

    long double** map_1 = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** map_2 = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** cos_lm = n_matrix_generator(nmod, nmod);
    long double** sin_lm = n_matrix_generator(nmod, nmod);
    long double** pml = n_matrix_generator(nmod, nmod);

    aml_gasdev(cos_lm, sin_lm, 0.0L, 1.0L);

    for (unsigned int m = 0; m < nmod; ++m) {
        for (unsigned int l = 0; l < nmod; ++l) {
            std::cout << "m = " << m << ", l = " << l << " alm: " << cos_lm[m][l] << ", " << sin_lm[m][l] << std::endl;
        }
    }

    fft_map_forward(map_1, cos_lm, sin_lm);
    fft_map_yy_forward(map_2, cos_lm, sin_lm);

//    for (unsigned int i = 0; i <= npix; ++i) {
//        for (unsigned int j = 0; j <= npix / 2; ++j) {
//            map_1[i][j] = map_1[i][j] + map_2[i][j];
//        }
//    }

    fft_map_backward(map_1, cos_lm, sin_lm);

    io_map("out.dat", map_1);

    std::cout << area(map_1, 0.0L) << std::endl;
    std::cout << length(map_1, 0.0L) << std::endl;

    std::cout << PI << std::endl;

    n_matrix_destroyer(pml, nmod);
    n_matrix_destroyer(cos_lm, nmod);
    n_matrix_destroyer(sin_lm, nmod);
    n_matrix_destroyer(map_1, npix + 1);
    n_matrix_destroyer(map_2, npix + 1);

    return 0;
}
