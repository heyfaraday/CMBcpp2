#include "math.h"

#include "pml.hpp"

#include "constants.hpp"
#include "parameters.hpp"

void pml_gen(unsigned int j, long double** array) {

    long double long_j = static_cast<long double>(j);

    for (unsigned int m = 0; m < nmod; ++m) {
        for (unsigned int l = 0; l < nmod; ++l) {
            array[m][l] = 0.0L;
        }
    }

    long double theta = 2.0L * long_j * PI / long_npix;

    array[0][0] = 1.0L / sqrtl(4.0L * PI);

    for (unsigned int m = 0; m < nmod - 1; ++m) {

        long double long_m = static_cast<long double>(m);

        array[m + 1][m + 1] = - array[m][m] * sinl(theta) * sqrtl(2.0L * long_m + 3.0L)
                              / sqrtl(2.0L * long_m + 2.0L);
    }

    for (unsigned int m = 0; m < nmod - 1; ++m) {

        long double long_m = static_cast<long double>(m);

        array[m][m + 1] = array[m][m] * cosl(theta) * sqrtl(2.0L * long_m + 3.0L);
    }

    for (unsigned int m = 0; m < nmod - 2; ++m) {
        for (unsigned int l = m + 2; l < nmod; ++l) {

            long double long_m = static_cast<long double>(m);
            long double long_l = static_cast<long double>(l);

            array[m][l] = ((2.0L * long_l - 1.0L) * sqrtl((long_l - long_m) * (2.0L * long_l + 1.0L))
                           / sqrtl((long_l + long_m) * (2.0L * l - 1.0L)) * array[m][l - 1] * cosl(theta)
                           - (long_l + long_m - 1.0L) * sqrtl((long_l - long_m) * (long_l - 1.0L - long_m)
                           * (2.0L * long_l + 1.0L)) / sqrtl((long_l + long_m) * (long_l - 1.0L + long_m)
                           * (2.0L * long_l - 3.0L)) * array[m][l - 2]) / (long_l - long_m);
        }
    }

    for (unsigned int m = 1; m < nmod ; ++m) {
        for (unsigned int l = m; l < nmod; ++l) {
            array[m][l] *= sqrtl(2.0L);
        }
    }
}


void pml_gen(long double theta, long double** array) {

    for (unsigned int m = 0; m < nmod; ++m) {
        for (unsigned int l = 0; l < nmod; ++l) {
            array[m][l] = 0.0L;
        }
    }

    array[0][0] = 1.0L / sqrtl(4.0L * PI);

    for (unsigned int m = 0; m < nmod - 1; ++m) {

        long double long_m = static_cast<long double>(m);

        array[m + 1][m + 1] = - array[m][m] * sinl(theta) * sqrtl(2.0L * long_m + 3.0L)
                              / sqrtl(2.0L * long_m + 2.0L);
    }

    for (unsigned int m = 0; m < nmod - 1; ++m) {

        long double long_m = static_cast<long double>(m);

        array[m][m + 1] = array[m][m] * cosl(theta) * sqrtl(2.0L * long_m + 3.0L);
    }

    for (unsigned int m = 0; m < nmod - 2; ++m) {
        for (unsigned int l = m + 2; l < nmod; ++l) {

            long double long_m = static_cast<long double>(m);
            long double long_l = static_cast<long double>(l);

            array[m][l] = ((2.0L * long_l - 1.0L) * sqrtl((long_l - long_m) * (2.0L * long_l + 1.0L))
                           / sqrtl((long_l + long_m) * (2.0L * l - 1.0L)) * array[m][l - 1] * cosl(theta)
                           - (long_l + long_m - 1.0L) * sqrtl((long_l - long_m) * (long_l - 1.0L - long_m)
                           * (2.0L * long_l + 1.0L)) / sqrtl((long_l + long_m) * (long_l - 1.0L + long_m)
                           * (2.0L * long_l - 3.0L)) * array[m][l - 2]) / (long_l - long_m);
        }
    }

    for (unsigned int m = 1; m < nmod ; ++m) {
        for (unsigned int l = m; l < nmod; ++l) {
            array[m][l] *= sqrtl(2.0L);
        }
    }
}
