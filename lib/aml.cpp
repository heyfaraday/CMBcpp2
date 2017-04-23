#include <random>
#include <iostream>

#include "aml.hpp"

#include "parameters.hpp"

void aml_gasdev(long double** cos_ml, long double** sin_ml, long double mean, long double std) {

    std::default_random_engine generator;
    //  std::default_random_engine generator(std::random_device{}());
    //  std::default_random_engine generator(1023515);
    //  or unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
    std::normal_distribution<long double> distribution(mean, std);

    for (unsigned int m = 0; m < nmod; ++m) {
        for (unsigned int l = m; l < nmod; ++l) {
            cos_ml[m][l] = distribution(generator);
        }
    }

    if (nmod > npix / 2 + 1) {
        std::cout << "Nyquist error!" << std::endl;
        exit (EXIT_FAILURE);
    }

    if ((nmod < 2) or (nback < 2)) {
        std::cout << "Not enough mods!" << std::endl;
        exit (EXIT_FAILURE);
    }

    if (nback < nmod) {
        std::cout << "NBACK < NMOD!" << std::endl;
        exit (EXIT_FAILURE);
    }

    for (unsigned int m = 1; m < nmod; ++m) {
        for (unsigned int l = m; l < nmod; ++l) {
            sin_ml[m][l] = distribution(generator);
        }
    }

    if (nmod == npix / 2 + 1) {
        std::cout << "Nyquist warning!" << std::endl;
        sin_ml[npix / 2][npix / 2] = 0.0L;
    }

    cos_ml[0][0] = 0.0L;
    cos_ml[0][1] = 0.0L;
    cos_ml[1][1] = 0.0L;

    sin_ml[0][0] = 0.0L;
    sin_ml[0][1] = 0.0L;
    sin_ml[1][1] = 0.0L;
}

void aml_from_cl(long double** cos_ml, long double** sin_ml, long double* cl) {

    std::default_random_engine generator;
    //  std::default_random_engine generator(std::random_device{}());
    //  std::default_random_engine generator(1023515);
    //  or unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
    std::normal_distribution<long double> distribution(0.0L, 1.0L);

    for (unsigned int m = 0; m < nmod; ++m) {
        for (unsigned int l = m; l < nmod; ++l) {
            cos_ml[m][l] = distribution(generator) * cl[l];
        }
    }

    if (nmod > npix / 2 + 1) {
        std::cout << "Nyquist error!" << std::endl;
        exit (EXIT_FAILURE);
    }

    if ((nmod < 2) or (nback < 2)) {
        std::cout << "Not enough mods!" << std::endl;
        exit (EXIT_FAILURE);
    }

    if (nback < nmod) {
        std::cout << "NBACK < NMOD!" << std::endl;
        exit (EXIT_FAILURE);
    }

    for (unsigned int m = 1; m < nmod; ++m) {
        for (unsigned int l = m; l < nmod; ++l) {
            sin_ml[m][l] = distribution(generator) * cl[l];
        }
    }

    if (nmod == npix / 2 + 1) {
        std::cout << "Nyquist warning!" << std::endl;
        sin_ml[npix / 2][npix / 2] = 0.0L;
    }

    cos_ml[0][0] = 0.0L;
    cos_ml[0][1] = 0.0L;
    cos_ml[1][1] = 0.0L;

    sin_ml[0][0] = 0.0L;
    sin_ml[0][1] = 0.0L;
    sin_ml[1][1] = 0.0L;
}

void healpix_transform(long double** aml) {

    for (unsigned int m = 0; m < nmod; ++m) {
        for (unsigned int l = 0; l < nmod; ++l) {
            if (l % 2 == 1) {
                aml[m][l] = - aml[m][l];
            }
        }
    }
}

void add_monopole(long double** cos_ml, long double** sin_ml, long double aml_0_0) {
    cos_ml[0][0] = aml_0_0;
    sin_ml[0][0] = 0.0L;
}

void remove_monopole(long double** cos_ml, long double** sin_ml) {
    cos_ml[0][0] = 0.0L;
    sin_ml[0][0] = 0.0L;
}

void add_dipole(long double** cos_ml, long double** sin_ml,
                long double aml_cos_0_1, long double aml_cos_1_1, long double aml_sin_1_1) {
    cos_ml[0][1] = aml_cos_0_1;
    sin_ml[0][1] = 0.0L;
    cos_ml[1][1] = aml_cos_1_1;
    sin_ml[1][1] = aml_sin_1_1;
}

void remove_dipole(long double** cos_ml, long double** sin_ml) {
    cos_ml[0][1] = 0.0L;
    sin_ml[0][1] = 0.0L;
    cos_ml[1][1] = 0.0L;
    sin_ml[1][1] = 0.0L;
}
