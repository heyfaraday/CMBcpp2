#include <iostream>
#include <fstream>
#include <chealpix.h>
#include <cmath>
#include <random>

#include <fft.hpp>
#include <utils.hpp>
#include <parameters.hpp>
#include "constants.hpp"
#include <io.hpp>
#include <monte.hpp>
#include "pml.hpp"
#include "aml.hpp"
#include "functionals.hpp"

int main() {

    // unsigned int seed1 = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
    //  std::default_random_engine generator(std::random_device{}());
    // std::default_random_engine generator(1023515);
    //  or

    typedef std::numeric_limits<long double> dbl;
    std::cout.precision(dbl::max_digits10);

    long double** map = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** map_x = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** map_y = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** cos_ml = n_matrix_generator(nmod, nmod);
    long double** sin_ml = n_matrix_generator(nmod, nmod);

    void level_points(long double** map_x, long double** map_y, long double** cos_ml, long double** sin_ml,
                      std::string name);

    i_aml("alm_T0B_T_cos.dat", cos_ml);
    i_aml("alm_T0B_T_sin.dat", sin_ml);

    remove_dipole(cos_ml, sin_ml);
    remove_monopole(cos_ml, sin_ml);

    std::cout << "start" << std::endl;

    fft_map_forward(map, cos_ml, sin_ml);
//    fft_map_x_forward(map_x, cos_ml, sin_ml);
//    fft_map_y_forward(map_y, cos_ml, sin_ml);

    level_area(map, "gasdev_area_test_2.dat");
    level_length(map, "gasdev_length_test_2.dat");

    std::cout << "end" << std::endl;

//    level_points(map_x, map_y, cos_ml, sin_ml, "gasdev_points_test_1.dat");

//    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
//
//    aml_gasdev(generator, cos_ml, sin_ml, 0.0L, 1.0L);
//
//    fft_map_forward(map, cos_ml, sin_ml);
//    fft_map_x_forward(map_x, cos_ml, sin_ml);
//    fft_map_y_forward(map_y, cos_ml, sin_ml);
//
//    level_area(map, "gasdev_area_test_1.dat");
//    level_length(map, "gasdev_length_test_1.dat");
//    level_points(map_x, map_y, cos_ml, sin_ml, "gasdev_points_test_1.dat");
//
//    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
//
//    aml_gasdev(generator, cos_ml, sin_ml, 0.0L, 1.0L);
//
//    fft_map_forward(map, cos_ml, sin_ml);
//    fft_map_x_forward(map_x, cos_ml, sin_ml);
//    fft_map_y_forward(map_y, cos_ml, sin_ml);
//
//    level_area(map, "gasdev_area_test_1.dat");
//    level_length(map, "gasdev_length_test_1.dat");
//    level_points(map_x, map_y, cos_ml, sin_ml, "gasdev_points_test_1.dat");

    n_matrix_destroyer(map, npix + 1);
    n_matrix_destroyer(map_x, npix + 1);
    n_matrix_destroyer(map_y, npix + 1);
    n_matrix_destroyer(cos_ml, nmod);
    n_matrix_destroyer(sin_ml, nmod);

    return 0;
}
