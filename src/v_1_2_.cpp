#include <iostream>
#include <random>
// #include "functionals_p.hpp"
#include "functionals.hpp"

#include "parameters.hpp"
#include "utils.hpp"
#include "aml.hpp"
#include "fft.hpp"
#include "io.hpp"
#include "spectra.hpp"
#include "monte.hpp"

int main() {

    typedef std::numeric_limits<long double> dbl;
    std::cout.precision(dbl::max_digits10);

    std::cout << "1" << std::endl;

    long double** u = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** q = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** u_x = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** q_x = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** u_y = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** q_y = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** u_xx = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** q_xx = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** u_yy = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** q_yy = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** u_xy = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** q_xy = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** cos_ml_u = n_matrix_generator(nmod, nmod);
    long double** cos_ml_q = n_matrix_generator(nmod, nmod);
    long double** sin_ml_u = n_matrix_generator(nmod, nmod);
    long double** sin_ml_q = n_matrix_generator(nmod, nmod);

    long double* cl = n_vector_generator(nmod);

    long double** p = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** p_x = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** p_y = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** p_xx = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** p_yy = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** p_xy = n_matrix_generator(npix + 1, npix / 2 + 1);

    long double** whitelist = n_matrix_generator(npix + 1, npix / 2 + 1);

//    std::default_random_engine generator;
//    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
//    aml_gasdev(generator, cos_ml_q, sin_ml_q, 0.0l, 1.0l);
//    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
//    aml_gasdev(generator, cos_ml_u, sin_ml_u, 0.0l, 1.0l);

    std::cout << "2" << std::endl;

    i_aml("q_cos_hfi.dat", cos_ml_q);
    i_aml("q_sin_hfi.dat", sin_ml_q);
//    i_aml("u_cos_hfi.dat", cos_ml_u);
//    i_aml("u_sin_hfi.dat", sin_ml_u);




    o_cl("test_cl_healpix", cl);

    std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
    aml_gasdev(generator, cos_ml_q, sin_ml_q, 0.0l, 1.0l);
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
    aml_gasdev(generator, cos_ml_u, sin_ml_u, 0.0l, 1.0l);

    std::cout << "3" << std::endl;

    fft_map_forward(q, cos_ml_q, sin_ml_q);
    fft_map_forward(u, cos_ml_u, sin_ml_u);

    std::cout << "4" << std::endl;

    fft_map_x_forward(q_x, cos_ml_q, sin_ml_q);
    fft_map_x_forward(u_x, cos_ml_u, sin_ml_u);
    fft_map_y_forward(q_y, cos_ml_q, sin_ml_q);
    fft_map_y_forward(u_y, cos_ml_u, sin_ml_u);

    std::cout << "5" << std::endl;

    fft_map_xx_forward(q_xx, cos_ml_q, sin_ml_q);
    fft_map_xx_forward(u_xx, cos_ml_u, sin_ml_u);
    fft_map_yy_forward(q_yy, cos_ml_q, sin_ml_q);
    fft_map_yy_forward(u_yy, cos_ml_u, sin_ml_u);
    fft_map_xy_forward(q_xy, cos_ml_q, sin_ml_q);
    fft_map_xy_forward(u_xy, cos_ml_u, sin_ml_u);

    std::cout << "6" << std::endl;


//    aml_to_cl(cos_ml_q, sin_ml_q, cl);

//    aml_to_healpix_cl(cos_ml_q, sin_ml_q, cl);
//
//    std::cout << "sigma_from_map " << sigma_0_map(q) << std::endl;
//    std::cout << "sigma from cl_healpix " << sigma_0_cl(cl) << std::endl;
//    std::cout << "sigma from normal aml " << sigma_0_aml(cos_ml_q, sin_ml_q) << std::endl;

    for (unsigned int i = 0; i < npix + 1; ++i) {
        for (unsigned int j = 0; j < npix / 2 + 1; ++j) {
            p[i][j] = sqrtl(q[i][j] * q[i][j] + u[i][j] * u[i][j]);
            p_x[i][j] = (q[i][j] * q_x[i][j] + u[i][j] * u_x[i][j]);
            p_y[i][j] = (q[i][j] * q_y[i][j] + u[i][j] * u_y[i][j]);
            p_xx[i][j] = (q_x[i][j] * q_x[i][j] + u_x[i][j] * u_x[i][j] + q[i][j] * q_xx[i][j] + u[i][j] * u_xx[i][j]);
            p_yy[i][j] = (q_y[i][j] * q_y[i][j] + u_y[i][j] * u_y[i][j] + q[i][j] * q_yy[i][j] + u[i][j] * u_yy[i][j]);
            p_xy[i][j] = (q_y[i][j] * q_x[i][j] + u_y[i][j] * u_x[i][j] + q[i][j] * q_xy[i][j] + u[i][j] * u_xy[i][j]);
        }
    }

    std::cout << "7" << std::endl;

    long double sigma = sigma_0_map(p);

    singular_points_classifier(q, u, q_x, u_x, q_y, u_y, "hfi_sing.dat", whitelist, nmod);

    std::cout << "8" << std::endl;

    points_classifier_p(p, p_x, p_y, p_xx, p_yy, p_xy, "hfi_points.dat", whitelist, nmod, sigma);

    o_map("map_test.dat", p);

    std::cout << "9" << std::endl;

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
    n_matrix_destroyer(cos_ml_u, nmod);
    n_matrix_destroyer(cos_ml_q, nmod);
    n_matrix_destroyer(sin_ml_u, nmod);
    n_matrix_destroyer(sin_ml_q, nmod);

    n_vector_destroyer(cl);

    n_matrix_destroyer(p, npix + 1);
    n_matrix_destroyer(p_x, npix + 1);
    n_matrix_destroyer(p_y, npix + 1);
    n_matrix_destroyer(p_xx, npix + 1);
    n_matrix_destroyer(p_yy, npix + 1);
    n_matrix_destroyer(p_xy, npix + 1);

    n_matrix_destroyer(whitelist, npix + 1);

    return 0;
}
