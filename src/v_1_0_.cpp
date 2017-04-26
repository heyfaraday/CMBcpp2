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
#include "pml.hpp"
#include "aml.hpp"
#include "functionals.hpp"
#include "monte.hpp"
#include "spectra.hpp"

int main() {

    typedef std::numeric_limits<long double> dbl;
    std::cout.precision(dbl::max_digits10);

    long double** map = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** cos_ml = n_matrix_generator(nmod, nmod);
    long double** sin_ml = n_matrix_generator(nmod, nmod);
    long double* cl = n_vector_generator(nmod);

    i_aml("alm_T0B_T_cos.dat", cos_ml);
    i_aml("alm_T0B_T_sin.dat", sin_ml);

    remove_dipole(cos_ml, sin_ml);
    remove_monopole(cos_ml, sin_ml);

    aml_to_cl(cos_ml, sin_ml, cl);

    std::default_random_engine generator;


    for (unsigned int global_monter_carlo = 0; global_monter_carlo < nmonte; ++global_monter_carlo) {
        std::cout << "Global start: " << global_monter_carlo << std::endl;

        generator.seed(std::chrono::system_clock::now().time_since_epoch().count());

        aml_from_cl(generator, cos_ml, sin_ml, cl);

        fft_map_forward(map, cos_ml, sin_ml);
        level_area(map, "gasdev_area_test_3.dat");
        level_length(map, "gasdev_length_test_3.dat");

        std::cout << "Global end: " << global_monter_carlo << std::endl;
    }

    n_matrix_destroyer(map, npix + 1);
    n_matrix_destroyer(cos_ml, nmod);
    n_matrix_destroyer(sin_ml, nmod);

    return 0;
}
