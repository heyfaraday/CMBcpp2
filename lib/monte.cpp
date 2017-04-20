#include <string>
#include "math.h"
#include <iostream>

#include "functionals.hpp"
#include "parameters.hpp"
#include "constants.hpp"

#include "monte.hpp"

void circles_area(long double** map, long double level, std::string name) {

    long double theta;
    unsigned int j1 = npix / 2;
    unsigned int j2 = npix / 2;
    unsigned int j_center = 0;
    long double v1, v2;

    for (unsigned int i_ring = 0; i_ring < nring / 2; ++i_ring) {

        theta = 2.0L * (long_npix / 4.0L + j_center) * PI / long_npix;

        if (j_center + static_cast<unsigned int>(sinl(theta)) * hring >= npix / 2) {
            std::cout << "Error" << std::endl;
            exit (EXIT_FAILURE);
        }

        v1 = area(map, level, npix / 4 + j_center, npix / 4 + j_center + static_cast<unsigned int>(hring * sinl(theta)));
        v2 = area(map, level, npix / 4 - j_center, npix / 4 - j_center - static_cast<unsigned int>(hring * sinl(theta)));

        j_center = j_center + static_cast<unsigned int>(hring * sinl(theta));

    }

}