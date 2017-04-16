#include <fftw3.h>
#include <string>

#include "fft.hpp"

#include "parameters.hpp"

void fft_map_forward(double** map, double** cos_lm, double** sin_lm)
{

    for (int i = 0; i <= npix; ++i) {
        for (int j = 0; j <= npix/2; ++j) {
            map[i][j] = 0.0;
        }
    }
}
