#include <iostream>
#include <random>

#include "parameters.hpp"
#include "utils.hpp"
#include "aml.hpp"
#include "fft.hpp"
#include "io.hpp"
#include "spectra.hpp"
#include "monte.hpp"


int main() {

    std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());

    typedef std::numeric_limits<long double> dbl;
    std::cout.precision(dbl::max_digits10);

    long double** p = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** q = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** u = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** map_out = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** cos_q = n_matrix_generator(nmod, nmod);
    long double** sin_q = n_matrix_generator(nmod, nmod);
    long double** cos_u = n_matrix_generator(nmod, nmod);
    long double** sin_u = n_matrix_generator(nmod, nmod);

    aml_gasdev(generator, cos_q, sin_q, 0.0l, 1.0l);
    aml_gasdev(generator, cos_u, sin_u, 0.0l, 1.0l);

//    i_aml("../data/CMB_q_cos_smooth.dat", cos_q);
//    i_aml("../data/CMB_q_sin_smooth.dat", sin_q);
//    i_aml("../data/CMB_u_cos_smooth.dat", cos_u);
//    i_aml("../data/CMB_u_sin_smooth.dat", sin_u);

    fft_map_forward(q, cos_q, sin_q);
    fft_map_forward(u, cos_u, sin_u);

    for (int i = 0; i <= npix; ++i) {
        for (int j = 0; j <= npix / 2; ++j) {
            p[i][j] = sqrtl(q[i][j] * q[i][j] + u[i][j] * u[i][j]);
        }
    }

    long double sigma = sigma_0_map(p);

//    circles_area_map(p, 0.03 * sigma, map_out);
//    o_map("map_test_area1.dat", map_out);
//    circles_area_map(p, 0.04 * sigma, map_out);
//    o_map("map_test_area2.dat", map_out);
//    circles_area_map(p, 0.05 * sigma, map_out);
//    o_map("map_test_area3.dat", map_out);
//    circles_area_map(p, 0.06 * sigma, map_out);
//    o_map("map_test_area4.dat", map_out);
//    circles_area_map(p, 0.07 * sigma, map_out);
//    o_map("map_test_area5.dat", map_out);
//    circles_area_map(p, 0.08 * sigma, map_out);
//    o_map("map_test_area6.dat", map_out);
//    circles_area_map(p, 0.09 * sigma, map_out);
//    o_map("map_test_area7.dat", map_out);
//    circles_area_map(p, 0.1 * sigma, map_out);
//    o_map("map_test_area8.dat", map_out);
//    circles_area_map(p, 0.11 * sigma, map_out);
//    o_map("map_test_area9.dat", map_out);
//    circles_area_map(p, 0.12 * sigma, map_out);
//    o_map("map_test_area10.dat", map_out);
//    circles_area_map(p, 0.13 * sigma, map_out);
//    o_map("map_test_area11.dat", map_out);
//    circles_area_map(p, 0.14 * sigma, map_out);
//    o_map("map_test_area12.dat", map_out);
//    circles_area_map(p, 0.15 * sigma, map_out);
//    o_map("map_test_area13.dat", map_out);
//    circles_area_map(p, 0.16 * sigma, map_out);
//    o_map("map_test_area14.dat", map_out);
//    circles_area_map(p, 0.17 * sigma, map_out);
//    o_map("map_test_area15.dat", map_out);
//    circles_area_map(p, 0.18 * sigma, map_out);
//    o_map("map_test_area16.dat", map_out);
//    circles_area_map(p, 0.19 * sigma, map_out);
//    o_map("map_test_area17.dat", map_out);
//    circles_area_map(p, 0.20 * sigma, map_out);
//    o_map("map_test_area18.dat", map_out);
//    circles_area_map(p, 0.22 * sigma, map_out);
//    o_map("map_test_area19.dat", map_out);
//    circles_area_map(p, 0.24 * sigma, map_out);
//    o_map("map_test_area20.dat", map_out);
//    circles_area_map(p, 0.26 * sigma, map_out);
//    o_map("map_test_area21.dat", map_out);
//    circles_area_map(p, 0.28 * sigma, map_out);
//    o_map("map_test_area22.dat", map_out);
//    circles_area_map(p, 0.30 * sigma, map_out);
//    o_map("map_test_area23.dat", map_out);

//    circles_area_map(p, 0.32 * sigma, map_out);
//    o_map("map_test_area24.dat", map_out);
//    circles_area_map(p, 0.34 * sigma, map_out);
//    o_map("map_test_area25.dat", map_out);
//    circles_area_map(p, 0.36 * sigma, map_out);
//    o_map("map_test_area26.dat", map_out);
//    circles_area_map(p, 0.38 * sigma, map_out);
//    o_map("map_test_area27.dat", map_out);
//    circles_area_map(p, 0.40 * sigma, map_out);
//    o_map("map_test_area28.dat", map_out);
//    circles_area_map(p, 0.42 * sigma, map_out);
//    o_map("map_test_area29.dat", map_out);
//    circles_area_map(p, 0.44 * sigma, map_out);
//    o_map("map_test_area30.dat", map_out);
//    circles_area_map(p, 0.46 * sigma, map_out);
//    o_map("map_test_area31.dat", map_out);
//    circles_area_map(p, 0.48 * sigma, map_out);
//    o_map("map_test_area32.dat", map_out);
//    circles_area_map(p, 0.50 * sigma, map_out);
//    o_map("map_test_area33.dat", map_out);
//    circles_area_map(p, 0.54 * sigma, map_out);
//    o_map("map_test_area34.dat", map_out);
//    circles_area_map(p, 0.58 * sigma, map_out);
//    o_map("map_test_area35.dat", map_out);
//    circles_area_map(p, 0.62 * sigma, map_out);
//    o_map("map_test_area36.dat", map_out);
//    circles_area_map(p, 0.66 * sigma, map_out);
//    o_map("map_test_area37.dat", map_out);
//    circles_area_map(p, 0.70 * sigma, map_out);
//    o_map("map_test_area38.dat", map_out);
//    circles_area_map(p, 0.74 * sigma, map_out);
//    o_map("map_test_area39.dat", map_out);
//    circles_area_map(p, 0.78 * sigma, map_out);
//    o_map("map_test_area40.dat", map_out);
//    circles_area_map(p, 0.82 * sigma, map_out);
//    o_map("map_test_area41.dat", map_out);
//    circles_area_map(p, 0.86 * sigma, map_out);
//    o_map("map_test_area42.dat", map_out);
//    circles_area_map(p, 0.90 * sigma, map_out);
//    o_map("map_test_area43.dat", map_out);
//    circles_area_map(p, 0.94 * sigma, map_out);
//    o_map("map_test_area44.dat", map_out);
//    circles_area_map(p, 0.98 * sigma, map_out);
//    o_map("map_test_area45.dat", map_out);
//    circles_area_map(p, 1.02 * sigma, map_out);
//    o_map("map_test_area46.dat", map_out);

//    circles_area_map(p, 1.06 * sigma, map_out);
//    o_map("map_test_area47.dat", map_out);
//    circles_area_map(p, 1.1 * sigma, map_out);
//    o_map("map_test_area48.dat", map_out);
//    circles_area_map(p, 1.16 * sigma, map_out);
//    o_map("map_test_area49.dat", map_out);
//    circles_area_map(p, 1.22 * sigma, map_out);
//    o_map("map_test_area50.dat", map_out);
//    circles_area_map(p, 1.28 * sigma, map_out);
//    o_map("map_test_area51.dat", map_out);
//    circles_area_map(p, 1.34 * sigma, map_out);
//    o_map("map_test_area52.dat", map_out);
//    circles_area_map(p, 1.4 * sigma, map_out);
//    o_map("map_test_area53.dat", map_out);
//    circles_area_map(p, 1.46 * sigma, map_out);
//    o_map("map_test_area54.dat", map_out);
//    circles_area_map(p, 1.52 * sigma, map_out);
//    o_map("map_test_area55.dat", map_out);
//    circles_area_map(p, 1.58 * sigma, map_out);
//    o_map("map_test_area56.dat", map_out);
//    circles_area_map(p, 1.64 * sigma, map_out);
//    o_map("map_test_area57.dat", map_out);
//    circles_area_map(p, 1.7 * sigma, map_out);
//    o_map("map_test_area58.dat", map_out);
//    circles_area_map(p, 1.76 * sigma, map_out);
//    o_map("map_test_area59.dat", map_out);
//    circles_area_map(p, 1.82 * sigma, map_out);
//    o_map("map_test_area60.dat", map_out);
//    circles_area_map(p, 1.88 * sigma, map_out);
//    o_map("map_test_area61.dat", map_out);
//    circles_area_map(p, 1.94 * sigma, map_out);
//    o_map("map_test_area62.dat", map_out);
//    circles_area_map(p, 2.0 * sigma, map_out);
//    o_map("map_test_area63.dat", map_out);
//    circles_area_map(p, 2.06 * sigma, map_out);
//    o_map("map_test_area64.dat", map_out);
//    circles_area_map(p, 2.12 * sigma, map_out);
//    o_map("map_test_area65.dat", map_out);
//    circles_area_map(p, 2.18 * sigma, map_out);
//    o_map("map_test_area66.dat", map_out);
//    circles_area_map(p, 2.24 * sigma, map_out);
//    o_map("map_test_area67.dat", map_out);
//    circles_area_map(p, 2.3 * sigma, map_out);
//    o_map("map_test_area68.dat", map_out);
//    circles_area_map(p, 2.4 * sigma, map_out);
//    o_map("map_test_area69.dat", map_out);

    circles_length_map(p, 0.03 * sigma, map_out);
    o_map("map_test_area1.dat", map_out);
    circles_length_map(p, 0.04 * sigma, map_out);
    o_map("map_test_area2.dat", map_out);
    circles_length_map(p, 0.05 * sigma, map_out);
    o_map("map_test_area3.dat", map_out);
    circles_length_map(p, 0.06 * sigma, map_out);
    o_map("map_test_area4.dat", map_out);
    circles_length_map(p, 0.07 * sigma, map_out);
    o_map("map_test_area5.dat", map_out);
    circles_length_map(p, 0.08 * sigma, map_out);
    o_map("map_test_area6.dat", map_out);
    circles_length_map(p, 0.09 * sigma, map_out);
    o_map("map_test_area7.dat", map_out);
    circles_length_map(p, 0.1 * sigma, map_out);
    o_map("map_test_area8.dat", map_out);
    circles_length_map(p, 0.11 * sigma, map_out);
    o_map("map_test_area9.dat", map_out);
    circles_length_map(p, 0.12 * sigma, map_out);
    o_map("map_test_area10.dat", map_out);
    circles_length_map(p, 0.13 * sigma, map_out);
    o_map("map_test_area11.dat", map_out);
    circles_length_map(p, 0.14 * sigma, map_out);
    o_map("map_test_area12.dat", map_out);
    circles_length_map(p, 0.15 * sigma, map_out);
    o_map("map_test_area13.dat", map_out);
    circles_length_map(p, 0.16 * sigma, map_out);
    o_map("map_test_area14.dat", map_out);
    circles_length_map(p, 0.17 * sigma, map_out);
    o_map("map_test_area15.dat", map_out);
    circles_length_map(p, 0.18 * sigma, map_out);
    o_map("map_test_area16.dat", map_out);
    circles_length_map(p, 0.19 * sigma, map_out);
    o_map("map_test_area17.dat", map_out);
    circles_length_map(p, 0.20 * sigma, map_out);
    o_map("map_test_area18.dat", map_out);
    circles_length_map(p, 0.22 * sigma, map_out);
    o_map("map_test_area19.dat", map_out);
    circles_length_map(p, 0.24 * sigma, map_out);
    o_map("map_test_area20.dat", map_out);
    circles_length_map(p, 0.26 * sigma, map_out);
    o_map("map_test_area21.dat", map_out);
    circles_length_map(p, 0.28 * sigma, map_out);
    o_map("map_test_area22.dat", map_out);
    circles_length_map(p, 0.30 * sigma, map_out);
    o_map("map_test_area23.dat", map_out);

    circles_length_map(p, 0.32 * sigma, map_out);
    o_map("map_test_area24.dat", map_out);
    circles_length_map(p, 0.34 * sigma, map_out);
    o_map("map_test_area25.dat", map_out);
    circles_length_map(p, 0.36 * sigma, map_out);
    o_map("map_test_area26.dat", map_out);
    circles_length_map(p, 0.38 * sigma, map_out);
    o_map("map_test_area27.dat", map_out);
    circles_length_map(p, 0.40 * sigma, map_out);
    o_map("map_test_area28.dat", map_out);
    circles_length_map(p, 0.42 * sigma, map_out);
    o_map("map_test_area29.dat", map_out);
    circles_length_map(p, 0.44 * sigma, map_out);
    o_map("map_test_area30.dat", map_out);
    circles_length_map(p, 0.46 * sigma, map_out);
    o_map("map_test_area31.dat", map_out);
    circles_length_map(p, 0.48 * sigma, map_out);
    o_map("map_test_area32.dat", map_out);
    circles_length_map(p, 0.50 * sigma, map_out);
    o_map("map_test_area33.dat", map_out);
    circles_length_map(p, 0.54 * sigma, map_out);
    o_map("map_test_area34.dat", map_out);
    circles_length_map(p, 0.58 * sigma, map_out);
    o_map("map_test_area35.dat", map_out);
    circles_length_map(p, 0.62 * sigma, map_out);
    o_map("map_test_area36.dat", map_out);
    circles_length_map(p, 0.66 * sigma, map_out);
    o_map("map_test_area37.dat", map_out);
    circles_length_map(p, 0.70 * sigma, map_out);
    o_map("map_test_area38.dat", map_out);
    circles_length_map(p, 0.74 * sigma, map_out);
    o_map("map_test_area39.dat", map_out);
    circles_length_map(p, 0.78 * sigma, map_out);
    o_map("map_test_area40.dat", map_out);
    circles_length_map(p, 0.82 * sigma, map_out);
    o_map("map_test_area41.dat", map_out);
    circles_length_map(p, 0.86 * sigma, map_out);
    o_map("map_test_area42.dat", map_out);
    circles_length_map(p, 0.90 * sigma, map_out);
    o_map("map_test_area43.dat", map_out);
    circles_length_map(p, 0.94 * sigma, map_out);
    o_map("map_test_area44.dat", map_out);
    circles_length_map(p, 0.98 * sigma, map_out);
    o_map("map_test_area45.dat", map_out);
    circles_length_map(p, 1.02 * sigma, map_out);
    o_map("map_test_area46.dat", map_out);

    circles_length_map(p, 1.06 * sigma, map_out);
    o_map("map_test_area47.dat", map_out);
    circles_length_map(p, 1.1 * sigma, map_out);
    o_map("map_test_area48.dat", map_out);
    circles_length_map(p, 1.16 * sigma, map_out);
    o_map("map_test_area49.dat", map_out);
    circles_length_map(p, 1.22 * sigma, map_out);
    o_map("map_test_area50.dat", map_out);
    circles_length_map(p, 1.28 * sigma, map_out);
    o_map("map_test_area51.dat", map_out);
    circles_length_map(p, 1.34 * sigma, map_out);
    o_map("map_test_area52.dat", map_out);
    circles_length_map(p, 1.4 * sigma, map_out);
    o_map("map_test_area53.dat", map_out);
    circles_length_map(p, 1.46 * sigma, map_out);
    o_map("map_test_area54.dat", map_out);
    circles_length_map(p, 1.52 * sigma, map_out);
    o_map("map_test_area55.dat", map_out);
    circles_length_map(p, 1.58 * sigma, map_out);
    o_map("map_test_area56.dat", map_out);
    circles_length_map(p, 1.64 * sigma, map_out);
    o_map("map_test_area57.dat", map_out);
    circles_length_map(p, 1.7 * sigma, map_out);
    o_map("map_test_area58.dat", map_out);
    circles_length_map(p, 1.76 * sigma, map_out);
    o_map("map_test_area59.dat", map_out);
    circles_length_map(p, 1.82 * sigma, map_out);
    o_map("map_test_area60.dat", map_out);
    circles_length_map(p, 1.88 * sigma, map_out);
    o_map("map_test_area61.dat", map_out);
    circles_length_map(p, 1.94 * sigma, map_out);
    o_map("map_test_area62.dat", map_out);
    circles_length_map(p, 2.0 * sigma, map_out);
    o_map("map_test_area63.dat", map_out);
    circles_length_map(p, 2.06 * sigma, map_out);
    o_map("map_test_area64.dat", map_out);
    circles_length_map(p, 2.12 * sigma, map_out);
    o_map("map_test_area65.dat", map_out);
    circles_length_map(p, 2.18 * sigma, map_out);
    o_map("map_test_area66.dat", map_out);
    circles_length_map(p, 2.24 * sigma, map_out);
    o_map("map_test_area67.dat", map_out);
    circles_length_map(p, 2.3 * sigma, map_out);
    o_map("map_test_area68.dat", map_out);
    circles_length_map(p, 2.4 * sigma, map_out);
    o_map("map_test_area69.dat", map_out);

    o_map("map_test.dat", p);

    n_matrix_destroyer(p, npix + 1);
    n_matrix_destroyer(q, npix + 1);
    n_matrix_destroyer(u, npix + 1);
    n_matrix_destroyer(map_out, npix + 1);
    n_matrix_destroyer(cos_q, nmod);
    n_matrix_destroyer(sin_q, nmod);
    n_matrix_destroyer(cos_u, nmod);
    n_matrix_destroyer(sin_u, nmod);

    return 0;
}

