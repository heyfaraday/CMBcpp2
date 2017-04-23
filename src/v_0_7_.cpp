#include <iostream>
#include <cmath>

#include "utils.hpp"
#include "parameters.hpp"
#include "fft.hpp"
#include "io.hpp"
#include "pml.hpp"
#include "constants.hpp"
#include "functionals.hpp"
#include "aml.hpp"

int main() {

    typedef std::numeric_limits<long double> dbl;
    std::cout.precision(dbl::max_digits10);

    long double** cos_ml = n_matrix_generator(nmod, nmod);
    long double** sin_ml = n_matrix_generator(nmod, nmod);
    long double** cos_ml_back = n_matrix_generator(nmod, nmod);
    long double** sin_ml_back = n_matrix_generator(nmod, nmod);
    long double** map = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** map_x = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** map_y = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** map_xx = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** map_yy = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** map_xy = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** pml = n_matrix_generator(nmod, nmod);
    long double** pml_plus = n_matrix_generator(nmod, nmod);
    long double** pml_minus = n_matrix_generator(nmod, nmod);

//    i_aml("alm_T0B_T_cos.dat", cos_ml);
//    i_aml("alm_T0B_T_sin.dat", sin_ml);

//    healpix_transform(cos_ml);
//    healpix_transform(sin_ml);


    add_dipole(cos_ml, sin_ml, 1.0L, 2.0L, 3.0L);
    cos_ml[0][2] = 12.0L;
    cos_ml[1][2] = 12.4L;
    cos_ml[2][2] = -5.0L;
    sin_ml[1][2] = 3.0L;
    sin_ml[2][2] = 0.5L;

    std::cout << "map_start" << std::endl;
    fft_map_forward(map, cos_ml, sin_ml);
    std::cout << "map_end" << std::endl;

    for (int m = 0; m < 3; ++m) {
        for (int l = 0; l < 3; ++l) {
            std::cout << "m: " << m << " l: " << l << " " << cos_ml[m][l] << " " << sin_ml[m][l] << std::endl;
        }
    }

    fft_map_backward(map, cos_ml, sin_ml);

    for (int m = 0; m < 3; ++m) {
        for (int l = 0; l < 3; ++l) {
            std::cout << "m: " << m << " l: " << l << " " << cos_ml[m][l] << " " << sin_ml[m][l] << std::endl;
        }
    }

    o_map("map_dohua.dat", map);

    n_matrix_destroyer(cos_ml, nmod);
    n_matrix_destroyer(sin_ml, nmod);
    n_matrix_destroyer(cos_ml_back, nmod);
    n_matrix_destroyer(sin_ml_back, nmod);
    n_matrix_destroyer(map, npix + 1);
    n_matrix_destroyer(map_x, npix + 1);
    n_matrix_destroyer(map_y, npix + 1);
    n_matrix_destroyer(map_xx, npix + 1);
    n_matrix_destroyer(map_yy, npix + 1);
    n_matrix_destroyer(map_xy, npix + 1);
    n_matrix_destroyer(pml, nmod);
    n_matrix_destroyer(pml_plus, nmod);
    n_matrix_destroyer(pml_minus, nmod);

    return 0;
}
