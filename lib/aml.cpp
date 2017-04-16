#include <random>
#include <iostream>

#include "aml.hpp"

#include "parameters.hpp"

void aml_gasdev(long double** cos_ml, long double** sin_ml, long double mean, long double std)
{


    std::default_random_engine generator;
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

    for (unsigned int m = 1; m < nmod; ++m) {
        for (unsigned int l = m; l < nmod; ++l) {
            sin_ml[m][l] = distribution(generator);
        }
    }

    if (nmod == npix / 2 + 1)
        std::cout << "Nyquist warning!" << std::endl;
        sin_ml[npix / 2][npix / 2] == 0.0;

    cos_ml[0][0] = 0.0;
    cos_ml[1][0] = 0.0;
    cos_ml[1][1] = 0.0;

    sin_ml[0][0] = 0.0;
    sin_ml[0][1] = 0.0;
    sin_ml[1][1] = 0.0;
}
