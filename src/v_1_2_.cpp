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
#include "distance.hpp"
#include "constants.hpp"

int main() {

    typedef std::numeric_limits<long double> dbl;
    std::cout.precision(dbl::max_digits10);

    long double** u = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** q = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** u_x = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** q_x = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** u_y = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** q_y = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** cos_lm_q = n_matrix_generator(nmod, nmod);
    long double** sin_lm_q = n_matrix_generator(nmod, nmod);
    long double** cos_lm_u = n_matrix_generator(nmod, nmod);
    long double** sin_lm_u = n_matrix_generator(nmod, nmod);

    long double** whitelist = n_matrix_generator(npix + 1, npix / 2 + 1);

    std::default_random_engine generator;

    std::cout << "1" << std::endl;

//    i_aml("q_cos_hfi.dat", cos_lm_q);
//    i_aml("q_sin_hfi.dat", sin_lm_q);
//    i_aml("u_cos_hfi.dat", cos_lm_u);
//    i_aml("u_sin_hfi.dat", sin_lm_u);

    aml_gasdev(generator, cos_lm_q, sin_lm_q, 0.0l, 1.0l);
    std::default_random_engine generator2;
    aml_gasdev(generator2, cos_lm_u, sin_lm_u, 0.0l, 1.0l);

//    std::cout << sin_lm_q[0][0] << std::endl;
//    std::cout << sin_lm_q[0][1] << std::endl;
//    std::cout << sin_lm_q[1][1] << std::endl;
//    std::cout << sin_lm_q[0][2] << std::endl;
//    std::cout << sin_lm_q[1][2] << std::endl;
//    std::cout << sin_lm_q[2][2] << std::endl;

//    i_aml("u_cos_hfi.dat", cos_lm_u);
//    i_aml("u_sin_hfi.dat", sin_lm_u);

    std::cout << "2" << std::endl;
//
    remove_monopole(cos_lm_q, sin_lm_q);
    remove_dipole(cos_lm_q, sin_lm_q);
    remove_monopole(cos_lm_u, sin_lm_u);
    remove_dipole(cos_lm_u, sin_lm_u);

//    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
//    aml_gasdev(generator, cos_lm_q, sin_lm_q, 0.0L, 1.0L);
//    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
//    aml_gasdev(generator, cos_lm_u, sin_lm_u, 0.0L, 1.0L);

    fft_map_forward(q, cos_lm_q, sin_lm_q);
    fft_map_forward(u, cos_lm_u, sin_lm_u);
    fft_map_x_forward(q_x, cos_lm_q, sin_lm_q);
    fft_map_y_forward(q_y, cos_lm_q, sin_lm_q);
    fft_map_x_forward(u_x, cos_lm_u, sin_lm_u);
    fft_map_y_forward(u_y, cos_lm_u, sin_lm_u);

    singular_points_classifier(q, u, q_x, u_x, q_y, u_y, "hfi_sing.dat", whitelist, nmod);

    long double** u_xx = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** q_xx = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** u_yy = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** q_yy = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** u_xy = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** q_xy = n_matrix_generator(npix + 1, npix / 2 + 1);

    fft_map_xx_forward(q_xx, cos_lm_q, sin_lm_q);
    fft_map_yy_forward(q_yy, cos_lm_q, sin_lm_q);
    fft_map_xy_forward(q_xy, cos_lm_q, sin_lm_q);
    fft_map_xx_forward(u_xx, cos_lm_q, sin_lm_q);
    fft_map_yy_forward(u_yy, cos_lm_q, sin_lm_q);
    fft_map_xy_forward(u_xy, cos_lm_q, sin_lm_q);

    std::cout << "3" << std::endl;

    long double** p = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** p_x = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** p_y = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** p_xx = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** p_yy = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** p_xy = n_matrix_generator(npix + 1, npix / 2 + 1);

    for (unsigned int i = 0; i < npix; ++i) {
        for (unsigned int j = 1; j < npix / 2; ++j) {
            p[i][j] = sqrtl(q[i][j] * q[i][j] + u[i][j] * u[i][j]);
            p_x[i][j] = (q[i][j] * q_x[i][j] + u[i][j] * u_x[i][j]);
            p_y[i][j] = (q[i][j] * q_y[i][j] + u[i][j] * u_y[i][j]);
            p_xx[i][j] = (q_x[i][j] * q_x[i][j] + u_x[i][j] * u_x[i][j] + q[i][j] * q_xx[i][j] + u[i][j] * u_xx[i][j]);
            p_yy[i][j] = (q_y[i][j] * q_y[i][j] + u_y[i][j] * u_y[i][j] + q[i][j] * q_yy[i][j] + u[i][j] * u_yy[i][j]);
            p_xy[i][j] = (q_y[i][j] * q_x[i][j] + u_y[i][j] * u_x[i][j] + q[i][j] * q_xy[i][j] + u[i][j] * u_xy[i][j]);
        }
    }

    o_map("map_test.dat",p);

    n_matrix_destroyer(u, npix + 1);
    n_matrix_destroyer(q, npix + 1);
    n_matrix_destroyer(u_x, npix + 1);
    n_matrix_destroyer(q_x, npix + 1);
    n_matrix_destroyer(u_y, npix + 1);
    n_matrix_destroyer(q_y, npix + 1);
    n_matrix_destroyer(u_xx, npix + 1);
    n_matrix_destroyer(q_xx, npix + 1);
    n_matrix_destroyer(u_yy, npix + 1);
    n_matrix_destroyer(q_yy, npix + 1);
    n_matrix_destroyer(u_xy, npix + 1);
    n_matrix_destroyer(q_xy, npix + 1);
    n_matrix_destroyer(cos_lm_q, nmod);
    n_matrix_destroyer(sin_lm_q, nmod);
    n_matrix_destroyer(cos_lm_u, nmod);
    n_matrix_destroyer(sin_lm_u, nmod);
//
    std::cout << "4" << std::endl;
//


    long double sigma = sigma_0_map(p);



    std::cout << "5" << std::endl;

    points_classifier_p(p, p_x, p_y, p_xx, p_yy, p_xy, "hfi_points.dat", whitelist, nmod, sigma);

    std::cout << q[100][100] << std::endl;


//
    n_matrix_destroyer(p, npix + 1);
    n_matrix_destroyer(p_x, npix + 1);
    n_matrix_destroyer(p_y, npix + 1);
    n_matrix_destroyer(p_xx, npix + 1);
    n_matrix_destroyer(p_yy, npix + 1);
    n_matrix_destroyer(p_xy, npix + 1);

    n_matrix_destroyer(whitelist, nmod);

    return 0;
}
