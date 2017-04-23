#include <iostream>
#include <cmath>

#include "utils.hpp"
#include "parameters.hpp"
#include "fft.hpp"
#include "io.hpp"
#include "pml.hpp"
#include "constants.hpp"

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

    // ! forw and back test
    cos_ml[1][2] = 1.0L;

    fft_map_forward(map, cos_ml, sin_ml);

    fft_map_backward(map, cos_ml_back, sin_ml_back);

    std::cout << "back_1_2: " << cos_ml_back[1][2] << std::endl;

    o_map("map.dat", map);
    //


    // ! norm test
    long double sum = 0.0L;
    long double norm = 0.0L;

    unsigned int my_m = 3;
    unsigned int my_l = 3;

    for (unsigned int j = 1; j < npix / 2; ++j) {

        y_plus_gen(j, pml_plus);
        y_minus_gen(j, pml_minus);

        sum = sum + pml_plus[my_m-1][my_l] * pml_plus[my_m][my_l] * sinl(map_parameter * static_cast<long double>(j));
        norm = norm + sinl(map_parameter * static_cast<long double>(j));

    }

    std::cout << "sum test: " << sum / norm * 4.0L * PI << std::endl;
//
//    y_plus_gen(70, pml);
//
//    long double long_j = 70;
//
//    std::cout << "test_y " << pml[2][2] << std::endl;
//    std::cout << "test_y_2 " << 1.0L / 8.0L * sqrtl(5.0L / PI) * (1 - cosl(2.0L * long_j * PI / long_npix)) * (1 - cosl(2.0L * long_j * PI / long_npix)) << std::endl;
//

    // след тес
    y_plus_gen(20, pml);
    y_plus_gen(20, pml);

    long double long_j = 20;

    std::cout << "test_y " << pml[1][2] << std::endl;
    std::cout << "test_y_2 " << 1.0L / 4.0L * sqrtl(5.0 / PI) * sinl(2.0L * long_j * PI / long_npix) * (1 - cosl(2.0L * long_j * PI / long_npix)) << std::endl;



    //



    std::cout << "sum test: " << sum / norm * 4.0L * PI << std::endl;
    //


    // ! read and write with planck test

    //

    // ! laplace test

        fft_map_x_forward(map_x, cos_ml, sin_ml);
        fft_map_y_forward(map_y, cos_ml, sin_ml);
        fft_map_xx_forward(map_xx, cos_ml, sin_ml);
        fft_map_yy_forward(map_yy, cos_ml, sin_ml);
        fft_map_xy_forward(map_xy, cos_ml, sin_ml);

        std::cout << "map " << map[10][20] << std::endl;
        std::cout << "map_x " << map_x[10][20] << std::endl;
        std::cout << "map_y " << map_y[10][20] << std::endl;
        std::cout << "map_xx " << map_xx[10][20] << std::endl;
        std::cout << "map_yy " << map_yy[10][20] << std::endl;
        std::cout << "map_xy " << map_xy[10][20] << std::endl;
        std::cout << "map_xx + map_yy / 12 " << -(map_xx[10][20] + map_yy[10][20]) / 6.0L << std::endl;
    //

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
