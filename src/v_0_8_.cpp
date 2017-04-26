#include <iostream>
#include <cmath>
#include <random>

#include "utils.hpp"
#include "parameters.hpp"
#include "aml.hpp"
#include "fft.hpp"
#include "io.hpp"
#include "functionals.hpp"
#include "spectra.hpp"

int main() {

    typedef std::numeric_limits<long double> dbl;
    std::cout.precision(dbl::max_digits10);

    std::default_random_engine generator;

    long double** map = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** map_x = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** map_y = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** map_xx = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** map_yy = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** map_xy = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** q = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** q_x = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** q_y = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** u = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** u_x = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** u_y = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** p = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** p_x = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** p_y = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** cos_ml = n_matrix_generator(2 * nmod, 2 * nmod);
    long double** sin_ml = n_matrix_generator(2 * nmod, 2 * nmod);
    long double** cos_ml_q = n_matrix_generator(2 * nmod, 2 * nmod);
    long double** sin_ml_q = n_matrix_generator(2 * nmod, 2 * nmod);
    long double** cos_ml_u = n_matrix_generator(2 * nmod, 2 * nmod);
    long double** sin_ml_u = n_matrix_generator(2 * nmod, 2 * nmod);
    long double** whitelist = n_matrix_generator(npix + 1, npix / 2 + 1);

    std::cout << "1" << std::endl;

    aml_gasdev(generator, cos_ml, sin_ml, 0.0L, 1.0L);

    fft_map_forward(map, cos_ml, sin_ml);
    fft_map_x_forward(map_x, cos_ml, sin_ml);
    fft_map_y_forward(map_y, cos_ml, sin_ml);
    fft_map_xx_forward(map_xx, cos_ml, sin_ml);
    fft_map_yy_forward(map_yy, cos_ml, sin_ml, 2 * nmod);
    fft_map_xy_forward(map_xy, cos_ml, sin_ml);

    std::cout << "2" << std::endl;

    for (unsigned int i = 0; i < npix; ++i) {
        for (unsigned int j = 1; j < npix / 2; ++j) {
            q[i][j] = map_xx[i][j] - map_yy[i][j];
            u[i][j] = 2.0L * map_xy[i][j];
        }
    }

    std::cout << "3" << std::endl;

    fft_map_backward(q, cos_ml_q, sin_ml_q, 2 * nmod);
    fft_map_backward(u, cos_ml_u, sin_ml_u, 2 * nmod);

    fft_map_x_forward(q_x, cos_ml_q, sin_ml_q, 2 * nmod);
    fft_map_y_forward(q_y, cos_ml_q, sin_ml_q, 2 * nmod);
    fft_map_x_forward(u_x, cos_ml_u, sin_ml_u, 2 * nmod);
    fft_map_y_forward(u_y, cos_ml_u, sin_ml_u, 2 * nmod);

    fft_map_forward(q, cos_ml_q, sin_ml_q, 2 * nmod);
    fft_map_forward(u, cos_ml_u, sin_ml_u, 2 * nmod);

    std::cout << "4" << std::endl;

    for (unsigned int i = 0; i < npix; ++i) {
        for (unsigned int j = 1; j < npix / 2; ++j) {
            p[i][j] = sqrtl(q[i][j] * q[i][j] + u[i][j] * u[i][j]);
            p_x[i][j] = (q[i][j] * q_x[i][j] + u[i][j] * u_x[i][j]);
            p_y[i][j] = (q[i][j] * q_y[i][j] + u[i][j] * u_y[i][j]);
        }
    }

    std::cout << "5" << std::endl;

    o_map("v_0_8_1_.dat", p);

    singular_points_classifier(q, u, cos_ml_q, sin_ml_q, cos_ml_u, sin_ml_u, "v_0_8_1_sing_points.dat", whitelist, 2 * nmod);


    std::cout << "6" << std::endl;

    long double sigma_p = sigma_0_map(p);

    points_classifier_p(p_x, p_y, cos_ml_q, sin_ml_q, cos_ml_u, sin_ml_u, "v_0_8_1_points.dat", whitelist, 2 * nmod, sigma_p);

    std::cout << "7" << std::endl;


    std::cout << sigma_0_map(map) << std::endl;
    std::cout << sigma_0_aml(cos_ml, sin_ml) << std::endl;
    std::cout << "point: " << fft_point_x_forward(50 * map_parameter, 50 * map_parameter, cos_ml, sin_ml, 2 * nmod) << std::endl;
    std::cout << "map: " << map_x[50][50];

    std::cout << std::endl;
    std::cout << u_x[50][50] << std::endl;
    std::cout << fft_point_x_forward(50 * map_parameter, 50 * map_parameter, cos_ml_u,  sin_ml_u, 2 * nmod);

    std::cout << std::endl;
    std::cout << map_x[327][30] << std::endl;
    std::cout << 4.018463030851293594 / map_parameter << " " << 30 * map_parameter;

    n_matrix_destroyer(map, npix + 1);
    n_matrix_destroyer(map_x, npix + 1);
    n_matrix_destroyer(map_y, npix + 1);
    n_matrix_destroyer(map_xx, npix + 1);
    n_matrix_destroyer(map_yy, npix + 1);
    n_matrix_destroyer(map_xy, npix + 1);
    n_matrix_destroyer(q, npix + 1);
    n_matrix_destroyer(q_x, npix + 1);
    n_matrix_destroyer(q_y, npix + 1);
    n_matrix_destroyer(u, npix + 1);
    n_matrix_destroyer(u_x, npix + 1);
    n_matrix_destroyer(u_y, npix + 1);
    n_matrix_destroyer(p, npix + 1);
    n_matrix_destroyer(p_x, npix + 1);
    n_matrix_destroyer(p_y, npix + 1);
    n_matrix_destroyer(cos_ml, nmod);
    n_matrix_destroyer(sin_ml, nmod);
    n_matrix_destroyer(cos_ml_q, nmod);
    n_matrix_destroyer(sin_ml_q, nmod);
    n_matrix_destroyer(cos_ml_u, nmod);
    n_matrix_destroyer(sin_ml_u, nmod);
    n_matrix_destroyer(whitelist, npix + 1);

    return 0;
}
