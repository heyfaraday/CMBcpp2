#include <iostream>
#include <fstream>
#include <chealpix.h>
#include <cmath>
#include <fftw3.h>

#include <fft.hpp>
#include <utils.hpp>
#include <parameters.hpp>
#include "constants.hpp"
#include <io.hpp>
#include "pml.hpp"
#include "aml.hpp"
#include "functionals.hpp"
#include "monte.hpp"
#include "io.hpp"

int main() {

    typedef std::numeric_limits<long double> dbl;
    std::cout.precision(dbl::max_digits10);

    long double** map = n_matrix_generator(npix + 1, npix / 2 + 1);
    long double** pml = n_matrix_generator(nmod, nmod);
    long double** cos_ml = n_matrix_generator(nmod, nmod);
    long double** sin_ml = n_matrix_generator(nmod, nmod);
    aml_gasdev(cos_ml, sin_ml, 0.0L, 1.0L);
    fft_map_forward(map, cos_ml, sin_ml);

    std::cout << area(map, 0.0L) << std::endl;

    circles_area(map, 0.0L, "area.dat");

    std::cout << length(map, 0.0L) << std::endl;

    o_map("lol.dat", map);

    fftwl_complex* in = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*in) * 8));
    fftwl_complex* out = static_cast<fftwl_complex*>(fftwl_malloc(sizeof(*out) * 8));
    fftwl_plan p = fftwl_plan_dft_1d(int(8), in, out, FFTW_BACKWARD, FFTW_PATIENT);

    for (int i = 0; i < 8; ++i) {
        in[i][0] = 1.0L;
        out[i][0] = 0.0L;
        in[i][1] = 0.0L;
        out[i][1] = 0.0L;
    }

    in[0][0] = 1.0L;

    for (int i = 0; i < 8; ++i) {
        std::cout << "in. real " << in[i][0] << std::endl;
    }

    fftwl_execute(p);

    for (int i = 0; i < 8; ++i) {
        std::cout << "out. real " << out[i][0] << std::endl;
    }



    fftwl_destroy_plan(p);
    fftwl_free(in);
    fftwl_free(out);

    n_matrix_destroyer(map, npix + 1);
    n_matrix_destroyer(pml, nmod);
    n_matrix_destroyer(cos_ml, nmod);
    n_matrix_destroyer(sin_ml, nmod);
    return 0;
}
