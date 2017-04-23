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

    long double** map = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** map_x = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** map_y = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** map_xx = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** map_yy = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** map_xy = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** cos_ml = n_matrix_generator(nmod, nmod);
    long double** sin_ml = n_matrix_generator(nmod, nmod);
    long double** pml = n_matrix_generator(nmod, nmod);

    aml_gasdev(cos_ml, sin_ml, 0.0L, 1.0L);

    fft_map_forward(map, cos_ml, sin_ml);
    fft_map_x_forward(map_x, cos_ml, sin_ml);
    fft_map_y_forward(map_y, cos_ml, sin_ml);
    fft_map_xx_forward(map_xx, cos_ml, sin_ml);
    fft_map_yy_forward(map_yy, cos_ml, sin_ml);
    fft_map_xy_forward(map_xy, cos_ml, sin_ml);

    points_classifier(map_x, map_y, cos_ml, sin_ml, "points_out.dat");

    int i = 200;
    int j = 200;
    std::cout << map_yy[i][j] << " " << fft_point_yy_forward(2.0L * PI * i / long_npix, 2.0L * PI * j / long_npix, cos_ml, sin_ml);

    o_map("out.dat", map);

    n_matrix_destroyer(pml, nmod);
    n_matrix_destroyer(cos_ml, nmod);
    n_matrix_destroyer(sin_ml, nmod);
    n_matrix_destroyer(map, npix + 1);
    n_matrix_destroyer(map_x, npix + 1);
    n_matrix_destroyer(map_y, npix + 1);
    n_matrix_destroyer(map_xx, npix + 1);
    n_matrix_destroyer(map_yy, npix + 1);
    n_matrix_destroyer(map_xy, npix + 1);

    return 0;
}

