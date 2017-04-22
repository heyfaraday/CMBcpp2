#include <string>
#include <cmath>
#include <fstream>
#include <iostream>

#include "functionals.hpp"
#include "parameters.hpp"
#include "constants.hpp"

#include "monte.hpp"

void circles_area(long double** map, long double level, std::string name) {

    std::ofstream out_file;
    out_file.open(name);

    typedef std::numeric_limits<long double> dbl;

    out_file.precision(dbl::max_digits10);

    long double theta;
    unsigned int j_center = 0;
    long double v1, v2;

    for (unsigned int i_ring = 0; i_ring < nring / 2; ++i_ring) {

        theta = 2.0L * (long_npix / 4.0L + j_center) * PI / long_npix;

        if (j_center + static_cast<unsigned int>(sinl(theta)) * hring >= npix / 2) {
            std::cout << "Error" << std::endl;
            exit (EXIT_FAILURE);
        }

        v1 = area(map, level, npix / 4 + j_center, npix / 4 + j_center + static_cast<unsigned int>(hring * sinl(theta)) - 1);
        v2 = area(map, level, npix / 4 - j_center - static_cast<unsigned int>(hring * sinl(theta)), npix / 4 - j_center - 1);

        out_file << std::scientific << v1 << " "
                 << npix / 4 + j_center << " "
                 << npix / 4 + j_center + static_cast<unsigned int>(hring * sinl(theta)) - 1
                 << std::endl
                 << std::scientific << v2 << " "
                 << npix / 4 - j_center - static_cast<unsigned int>(hring * sinl(theta)) << " "
                 << npix / 4 - j_center - 1
                 << std::endl;

        j_center = j_center + static_cast<unsigned int>(hring * sinl(theta));

    }

    out_file.close();
}

void circles_length(long double** map, long double level, std::string name) {

    std::ofstream out_file;
    out_file.open(name);

    typedef std::numeric_limits<long double> dbl;

    out_file.precision(dbl::max_digits10);

    long double theta;
    unsigned int j_center = 0;
    long double v1, v2;

    for (unsigned int i_ring = 0; i_ring < nring / 2; ++i_ring) {

        theta = 2.0L * (long_npix / 4.0L + j_center) * PI / long_npix;

        if (j_center + static_cast<unsigned int>(sinl(theta)) * hring >= npix / 2) {
            std::cout << "Error" << std::endl;
            exit (EXIT_FAILURE);
        }

        v1 = area(map, level, npix / 4 + j_center, npix / 4 + j_center + static_cast<unsigned int>(hring * sinl(theta)) - 1);
        v2 = area(map, level, npix / 4 - j_center - static_cast<unsigned int>(hring * sinl(theta)), npix / 4 - j_center - 1);

        out_file << std::scientific << v1 << " "
                 << npix / 4 + j_center << " "
                 << npix / 4 + j_center + static_cast<unsigned int>(hring * sinl(theta)) - 1
                 << std::endl
                 << std::scientific << v2 << " "
                 << npix / 4 - j_center - static_cast<unsigned int>(hring * sinl(theta)) << " "
                 << npix / 4 - j_center - 1
                 << std::endl;

        j_center = j_center + static_cast<unsigned int>(hring * sinl(theta));

    }

    out_file.close();
}

void circles_points(long double** map, long double level, std::string name)