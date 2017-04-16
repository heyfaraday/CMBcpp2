#include <iostream>
#include <fstream>
#include <chealpix.h>

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

    long double** map = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** cos_lm = n_matrix_generator(nmod + 1, nmod + 1);
    long double** sin_lm = n_matrix_generator(nmod + 1, nmod + 1);
    long double** pml = n_matrix_generator(nmod + 1, nmod + 1);

    aml_gasdev(cos_lm, sin_lm, 0.0L, 1.0L);

    fft_map_forward(map, cos_lm, sin_lm);

    fft_map_backward(map, cos_lm, sin_lm);

    io_map("out.dat", map);

    std::cout << area(map, 0.0L) << std::endl;

    std::cout << PI << std::endl;

    n_matrix_destroyer(pml, nmod + 1);
    n_matrix_destroyer(cos_lm, nmod + 1);
    n_matrix_destroyer(sin_lm, nmod + 1);
    n_matrix_destroyer(map, npix + 1);

    return 0;
}
