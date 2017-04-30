#include "utils.hpp"
#include "parameters.hpp"
#include "fft.hpp"
#include "io.hpp"
#include <iostream>

int main() {

    typedef std::numeric_limits<long double> dbl;
    std::cout.precision(dbl::max_digits10);

    long double** map = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** cos_ml = n_matrix_generator(nmod, nmod);
    long double** sin_ml = n_matrix_generator(nmod, nmod);

    sin_ml[1][1] = 1.0;
    cos_ml[0][1] = 1.0;
    cos_ml[0][0] = 1.0;
    sin_ml[1][1] = 1.0;

    fft_map_forward(map, cos_ml, sin_ml);



    std::cout << map[50][100] << std::endl;

    n_matrix_destroyer(map, npix + 1);
    n_matrix_destroyer(cos_ml, nmod);
    n_matrix_destroyer(sin_ml, nmod);

    return(0);
}