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

    long double** p = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** q = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** u = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** map_xx = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** map_yy = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** map_xy = n_matrix_generator(npix + 1, npix / 2 + 1);

    long double** cos_ml = n_matrix_generator(nback, nback);
    long double** sin_ml = n_matrix_generator(nback, nback);
    long double** q_cos_ml = n_matrix_generator(nback, nback);
    long double** q_sin_ml = n_matrix_generator(nback, nback);
    long double** u_cos_ml = n_matrix_generator(nback, nback);
    long double** u_sin_ml = n_matrix_generator(nback, nback);

    aml_gasdev(cos_ml, sin_ml, 0.0L, 1.0L);

    fft_map_xx_forward(map_xx, cos_ml, sin_ml);
    fft_map_yy_forward(map_yy, cos_ml, sin_ml);
    fft_map_xy_forward(map_xy, cos_ml, sin_ml);

    for (unsigned int i = 0; i <= npix; ++i) {
        for (unsigned int j = 0; j <= npix / 2; ++j) {
            q[i][j] = map_xx[i][j] - map_yy[i][j];
            u[i][j] = 2.0L * map_xy[i][j];
            p[i][j] = sqrtl(q[i][j] * q[i][j] + u[i][j] * u[i][j]);
        }
    }

    fft_map_backward(q, q_cos_ml, q_sin_ml, nback);
    fft_map_backward(u, u_cos_ml, u_sin_ml, nback);

    o_map("out.dat", p);

    singular_points_classifier(q, u, q_cos_ml, q_sin_ml, u_cos_ml, u_sin_ml, "singular_points_out.dat");

    n_matrix_destroyer(u_cos_ml, nback);
    n_matrix_destroyer(u_sin_ml, nback);
    n_matrix_destroyer(q_cos_ml, nback);
    n_matrix_destroyer(q_sin_ml, nback);
    n_matrix_destroyer(cos_ml, nback);
    n_matrix_destroyer(sin_ml, nback);

    n_matrix_destroyer(p, npix + 1);
    n_matrix_destroyer(q, npix + 1);
    n_matrix_destroyer(u, npix + 1);
    n_matrix_destroyer(map_xx, npix + 1);
    n_matrix_destroyer(map_yy, npix + 1);
    n_matrix_destroyer(map_xy, npix + 1);

    return 0;
}
