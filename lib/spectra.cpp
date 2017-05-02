#include <cmath>
#include <iostream>
#include "spectra.hpp"

#include "parameters.hpp"
#include "constants.hpp"
#include "utils.hpp"

long double sigma_0_map(long double** map) {

    long double sum_without_norm = 0.0L;
    long double norm = 0.0L;

    for (unsigned int i = 0; i < npix; ++i) {
        for (unsigned int j = 1; j < npix / 2 - 1; ++j) {

            if (!is_equal(map[i][j], 0.0L) and !is_equal(map[i + 1][j], 0.0L)
                and !is_equal(map[i][j + 1], 0.0L) and !is_equal(map[i + 1][j + 1], 0.0L)) {

                sum_without_norm = sum_without_norm + sinl(map_parameter * static_cast<long double>(j))
                                                      * map[i][j] * map[i][j];
                norm = norm + sinl(map_parameter * static_cast<long double>(j));
            }
        }
    }
    return sqrtl(sum_without_norm / norm);
}

long double sigma_0_cl(long double* cl) {

    long double sum = 0.0L;

    for (unsigned int l = 0; l < nmod; ++l) {
        sum = sum + (2.0L * l + 1.0L) * cl[l];
    }

    return sqrtl(sum / 4.0L / PI);
}

long double sigma_1_cl(long double* cl) {

    long double sum = 0.0L;

    for (unsigned int l = 0; l < nmod; ++l) {
        sum = sum + l * (l + 1.0L) * (2.0L * l + 1.0L) * cl[l];
    }

    return sqrtl(sum / 4.0L / PI);
}

long double sigma_2_cl(long double* cl) {

    long double sum = 0.0L;

    for (unsigned int l = 0; l < nmod; ++l) {
        sum = sum + (l + 2.0L) * (l - 1.0L) * l * (l + 1.0L) * (2.0L * l + 1.0L) * cl[l];
    }

    return sqrtl(sum / 4.0L / PI);
}

long double sigma_0_aml(long double** cos_ml, long double** sin_ml) {

    long double sum = cos_ml[0][0] * cos_ml[0][0];

    std::cout << std::endl << "I'm here " << sum << std::endl;

    for (unsigned int l = 1; l < nmod; ++l) {
        for (unsigned int m = 1; m < nmod; ++m) {
            sum = sum + (cos_ml[m][l] * cos_ml[m][l] + sin_ml[m][l] * sin_ml[m][l]) * 2.0L;
        }
        sum = sum + (cos_ml[0][l] * cos_ml[0][l] + sin_ml[0][l] * sin_ml[0][l]);
    }
    return sqrtl(sum / 4.0L / PI);
}

long double sigma_1_aml(long double** cos_ml, long double** sin_ml) {

    long double sum = 0.0L;

    for (unsigned int l = 0; l < nmod; ++l) {
        for (unsigned int m = 0; m < nmod; ++m) {
            sum = sum + l * (l + 1.0L) * cos_ml[m][l] * cos_ml[m][l] + sin_ml[m][l] * sin_ml[m][l];
        }
    }

    return sqrtl(sum / 4.0L / PI);
}

long double sigma_2_aml(long double** cos_ml, long double** sin_ml) {

    long double sum = 0.0L;

    for (unsigned int l = 0; l < nmod; ++l) {
        for (unsigned int m = 0; m < nmod; ++m) {
            sum = sum + (l + 2.0L) * (l - 1.0L) * l * (l + 1.0L) *
                                cos_ml[m][l] * cos_ml[m][l] + sin_ml[m][l] * sin_ml[m][l];
        }
    }

    return sqrtl(sum / 4.0L / PI);
}

void aml_to_healpix_cl(long double** cos_ml, long double** sin_ml, long double* cl) {

    for (unsigned int l = 0; l < nmod; ++l) {
        cl[l] = 0.0L;
    }

    for (unsigned int l = 0; l < nmod; ++l) {
        for (unsigned int m = 1; m < nmod; ++m) {
            cl[l] = cl[l] + (cos_ml[m][l] * cos_ml[m][l] + sin_ml[m][l] * sin_ml[m][l]) * 2.0L;  // * 2.0 from healpix
        }
        cl[l] = cl[l] + (cos_ml[0][l] * cos_ml[0][l]);
        cl[l] = cl[l] / (2.0L * l + 1.0L);
    }
}

void aml_to_cl(long double** cos_ml, long double** sin_ml, long double* cl) {

    for (unsigned int l = 0; l < nmod; ++l) {
        cl[l] = 0.0L;
    }

    for (unsigned int l = 0; l < nmod; ++l) {
        for (unsigned int m = 0; m < nmod; ++m) {
            cl[l] = cl[l] + (cos_ml[m][l] * cos_ml[m][l] + sin_ml[m][l] * sin_ml[m][l]);
        }
        cl[l] = cl[l] / (2.0L * l + 1.0L);
    }
}
