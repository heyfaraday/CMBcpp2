#include "math.h"

#include "pml.hpp"

#include "constants.hpp"
#include "parameters.hpp"

void pml_gen(unsigned int j, long double** array) {

    for (unsigned int m = 0; m < nmod; ++m) {
        for (unsigned int l = 0; l < nmod; ++l) {
            array[m][l] = 0.0;
        }
    }

    double theta = 2.0 * j * PI / npix;

    array[0][0] = 1.0 / sqrt(4.0 * PI);

    for (unsigned int m = 0; m < nmod - 1; ++m) {
        array[m + 1][m + 1] = - array[m][m] * sin(theta) * sqrt(2.0 * m + 3.0) / sqrt(2.0 * m + 2.0);
    }

    for (unsigned int m = 0; m < nmod - 1; ++m) {
        array[m][m + 1] = array[m][m] * cos(theta) * sqrt(2.0 * m + 3.0);
    }

    for (unsigned int m = 0; m < nmod - 2; ++m) {
        for (unsigned int l = m + 2; l < nmod; ++l) {
            array[m][l] = ((2.0 * l - 1.0) * sqrt((l - m) * (2.0 * l + 1.0)) / sqrt((l + m) * (2.0 * l - 1.0))
                           * array[m][l - 1] * cos(theta) - (l + m - 1.0) * sqrt((l - m) * (l - 1.0 - m)
                           * (2.0 * l + 1.0)) / sqrt((l + m) * (l - 1.0 + m) * (2.0 * l - 3.0))
                           * array[m][l - 2]) / (l - m);
        }
    }

    for (unsigned int m = 1; m < nmod ; ++m) {
        for (unsigned int l = m; l < nmod; ++l) {
            array[m][l] *= sqrt(2.0);
        }
    }
}


void pml_gen(double theta, long double** array) {

    for (unsigned int m = 0; m < nmod; ++m) {
        for (unsigned int l = 0; l < nmod; ++l) {
            array[m][l] = 0.0;
        }
    }

    array[0][0] = 1.0 / sqrt(4.0 * PI);

    for (unsigned int m = 0; m < nmod - 1; ++m) {
        array[m + 1][m + 1] = - array[m][m] * sin(theta) * sqrt(2.0 * m + 3.0) / sqrt(2.0 * m + 2.0);
    }

    for (unsigned int m = 0; m < nmod - 1; ++m) {
        array[m][m + 1] = array[m][m] * cos(theta) * sqrt(2.0 * m + 3.0);
    }

    for (unsigned int m = 0; m < nmod - 2; ++m) {
        for (unsigned int l = m + 2; l < nmod; ++l) {
            array[m][l] = ((2.0 * l - 1.0) * sqrt((l - m) * (2.0 * l + 1.0)) / sqrt((l + m) * (2.0 * l - 1.0))
                           * array[m][l - 1] * cos(theta) - (l + m - 1.0) * sqrt((l - m) * (l - 1.0 - m)
                           * (2.0 * l + 1.0)) / sqrt((l + m) * (l - 1.0 + m) * (2.0 * l - 3.0))
                           * array[m][l - 2]) / (l - m);
        }
    }

    for (unsigned int m = 1; m < nmod ; ++m) {
        for (unsigned int l = m; l < nmod; ++l) {
            array[m][l] *= sqrt(2.0);
        }
    }
}
